#include "../../corelib/5mc_motif_finder.hpp"
#include "../../corelib/arg_parse.hpp"
#include "../../corelib/build_mod_bam.hpp"
#include "../../corelib/pdqsort.h"
#include "../../corelib/sam_batch.hpp"
#include "eval_kmer_features.hpp"
#include "mod_options.hpp"
#include "program_info.hpp"
#include <htslib/hts.h>

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <sys/stat.h>
#include <sys/utsname.h>
#include <dirent.h>

#include <ATen/Parallel.h>
#include <torch/csrc/api/include/torch/nn/functional/activation.h>
#include <torch/script.h>
#include <torch/cuda.h>

namespace ns_5mc_call {

static int
s_load_kmer_size_from_model_dir(const char* model_dir)
{
    int kmer_size = 0;
    char path[HBN_MAX_PATH_LEN];
    snprintf(path, HBN_MAX_PATH_LEN, "%s/kmer.txt", model_dir);
    hbn_dfopen(in, path, "r");
    HBN_SCANF(fscanf, in, 1, "%d", &kmer_size);
    hbn_fclose(in);
    return kmer_size;
}


static torch::jit::script::Module
s_load_model(const char* model_dir, const char* context)
{
    char path[HBN_MAX_PATH_LEN];
    snprintf(path, HBN_MAX_PATH_LEN, "%s/%s.pt", model_dir, context);
    if (access(path, F_OK)) HBN_LOG("Model file %s does not exist", path);

    if (!torch::cuda::is_available()) {
        throw std::runtime_error("CUDA device unavailable");
    }
    torch::jit::script::Module model = torch::jit::load(path, torch::kCUDA);
    model.eval();
    return model;
}

static void
add_cmd_to_sam_hdr(int argc, char* argv[], sam_hdr_t* hdr)
{
    std::ostringstream os;
    os << "@PG"
       << "\tID:" << HBN_PACKAGE_NAME
       << "\tPN:" << HBN_PACKAGE_NAME
       << "\tVN:" << HBN_PACKAGE_VERSION
       << "\tCL:";

    os << argv[0];
    for (int i = 1; i < argc; ++i) os << ' ' << argv[i];
    os << '\n';

    std::string cmd = os.str();
    sam_hdr_add_lines(hdr, cmd.c_str(), cmd.size());
}

class ThreadWorkData
{
public:
    ThreadWorkData(ModOptions* options, int argc, char* argv[])
    {
        M_options = options;
        M_sam = new SAM_Batch(options->bam_path.c_str(), options->read_batch_size);

        sam_hdr_t* hdr_in = M_sam->sam_hdr();
        sam_hdr_t* hdr_out = sam_hdr_dup(hdr_in);
        M_out = sam_open(options->mod_bam_path.c_str(), "wb");
		hts_set_threads(M_out, 8);
		hts_set_cache_size(M_out, 100000000);
        add_cmd_to_sam_hdr(argc, argv, hdr_out);
        if (sam_hdr_write(M_out, hdr_out)) HBN_ERR("Could not write BAM header to %s", options->mod_bam_path);
        sam_hdr_destroy(hdr_out);

        M_batch_cpg_samples = 0;
        M_batch_chg_samples = 0;
        M_batch_chh_samples = 0;

        A_cpg_samples = 0;
        A_chg_samples = 0;
        A_chh_samples = 0;
        A_num_reads = 0;
        A_num_bases = 0;

        M_thread_id = 0;

        if (!options->disable_cpg) M_cpg_model = s_load_model(options->model_dir.c_str(), "CpG");
        if (!options->disable_chg) M_chg_model = s_load_model(options->model_dir.c_str(), "CHG");
        if (!options->disable_chh) M_chh_model = s_load_model(options->model_dir.c_str(), "CHH");

	M_kmer_size = s_load_kmer_size_from_model_dir(options->model_dir.c_str());
    }

    ~ThreadWorkData()
    {
        delete M_sam;
        sam_close(M_out);

        HBN_LOG("***************************** Final sample stats:");
        std::string size = bytes_to_datasize(A_num_bases);        
        fprintf(stderr, "Reads: %zu\n", A_num_reads);
        fprintf(stderr, "Bases: %s\n", size.c_str());
        size = format_i64_with_commas(A_cpg_samples);
        fprintf(stderr, "CpG samples: %s\n", size.c_str());
        size = format_i64_with_commas(A_chg_samples);
        fprintf(stderr, "CHG samples: %s\n", size.c_str());
        size = format_i64_with_commas(A_chh_samples);
        fprintf(stderr, "CHH samples: %s\n", size.c_str());
    }

    bool get_next_sam(int& read_id, bam1_t* bam)
    {
        return M_sam->get_next_sam(read_id, bam);
    }

    bool reset_batch_read_idx()
    {
	    M_thread_id = 0;
        return M_sam->reset_batch_read_idx();
    }

    int thread_id()
    {
	    std::lock_guard<std::mutex> lg(M_model_mutex);
	    return M_thread_id++;
    }

    void add_mod_bams(std::pair<int, bam1_t*>* A, int C, size_t cpg_samples, size_t chg_samples, size_t chh_samples)
    {
        std::lock_guard<std::mutex> __(M_out_mutex);
        M_bam_list.insert(M_bam_list.end(), A, A + C);
        M_batch_cpg_samples += cpg_samples;
        M_batch_chg_samples += chg_samples;
        M_batch_chh_samples += chh_samples;
    }

    void save_mod_bam()
    {
        HBN_LOG("Save %zu mod bams", M_bam_list.size());
        std::sort(M_bam_list.begin(), M_bam_list.end(), 
            [](const std::pair<int, bam1_t*>& x, const std::pair<int, bam1_t*>& y) { return x.first < y.first; });
        sam_hdr_t* hdr = M_sam->sam_hdr();
        size_t batch_bases = 0;
        for (auto& rdmod : M_bam_list) {
            batch_bases += rdmod.second->core.l_qseq;
            if (sam_write1(M_out, hdr, rdmod.second) < 0) {
                HBN_ERR("Could not write BAM record to file");
            }
            bam_destroy1(rdmod.second);
        }

        std::string size = bytes_to_datasize(batch_bases);        
        fprintf(stderr, "===========> sample stats:\n");
        fprintf(stderr, "Reads: %zu\n", M_bam_list.size());
        fprintf(stderr, "Bases: %s\n", size.c_str());
        size = format_i64_with_commas(M_batch_cpg_samples);
        fprintf(stderr, "CpG samples: %s\n", size.c_str());
        size = format_i64_with_commas(M_batch_chg_samples);
        fprintf(stderr, "CHG samples: %s\n", size.c_str());
        size = format_i64_with_commas(M_batch_chh_samples);
        fprintf(stderr, "CHH samples: %s\n", size.c_str());

        A_num_reads += M_bam_list.size();
        A_num_bases += batch_bases;
        A_cpg_samples += M_batch_cpg_samples;
        A_chg_samples += M_batch_chg_samples;
        A_chh_samples += M_batch_chh_samples;

        M_batch_cpg_samples = 0;
        M_batch_chg_samples = 0;
        M_batch_chh_samples = 0;
        M_bam_list.clear();
    }

    void forward(std::vector<torch::jit::IValue>& inputs, torch::Tensor& output_cpu, EMethyContext ctx) {
        //std::lock_guard<std::mutex> __(M_model_mutex);
        try {
            //torch::NoGradGuard no_grad;
            if (ctx == EMethyContext::eCpG) {
                torch::Tensor output_gpu = M_cpg_model.forward(inputs).toTensor();
                output_cpu.copy_(output_gpu.to(torch::kCPU));
                //torch::cuda::synchronize();
            } else if (ctx == EMethyContext::eCHG) {
                torch::Tensor output_gpu = M_chg_model.forward(inputs).toTensor();
                output_cpu.copy_(output_gpu.to(torch::kCPU));
                //torch::cuda::synchronize();
            } else if (ctx == EMethyContext::eCHH) {
                torch::Tensor output_gpu = M_chh_model.forward(inputs).toTensor();
                output_cpu.copy_(output_gpu.to(torch::kCPU));
                //torch::cuda::synchronize();
            }
        } catch (const c10::Error& e) {
            std::cerr << "LibTorch error: " << e.what() << std::endl;
            abort();
        } catch(const std::exception& e) {
            std::cerr << "System error: " << e.what() << std::endl;
            abort();
        }
    }

public:
    ModOptions*        M_options;
    SAM_Batch*      M_sam;
    int             M_thread_id;

    size_t          M_batch_cpg_samples;
    size_t          M_batch_chg_samples;
    size_t          M_batch_chh_samples;

    size_t          A_cpg_samples;
    size_t          A_chg_samples;
    size_t          A_chh_samples;
    size_t          A_num_reads;
    size_t          A_num_bases;

    std::vector<std::pair<int, bam1_t*>> M_bam_list;
    std::mutex      M_model_mutex;

    std::mutex      M_out_mutex;
    samFile*        M_out;

    torch::jit::script::Module  M_cpg_model;
    torch::jit::script::Module  M_chg_model;
    torch::jit::script::Module  M_chh_model;
    int			M_kmer_size;
};

static void
fill_skipped_tags(std::set<int>& skipped_tags)
{
    int t0, t1, tx;

    t0 = 'f';
    t1 = 'p';
    tx = (t0 << 8) | t1;
    skipped_tags.insert(tx);

    t0 = 'r';
    t1 = 'p';
    tx = (t0 << 8) | t1;
    skipped_tags.insert(tx);

    t0 = 'f';
    t1 = 'i';
    tx = (t0 << 8) | t1;
    skipped_tags.insert(tx);

    t0 = 'r';
    t1 = 'i';
    tx = (t0 << 8) | t1;
    skipped_tags.insert(tx);
}

static void
s_call_one_batch_samples(ThreadWorkData* data,
    torch::TensorOptions& tensor_options,
    const EMethyContext ctx,
    std::vector<torch::jit::IValue>& inputs,
    torch::Tensor& prob_tensor_cpu,
    const int kmer_size,
    const int num_samples,
    std::vector<MolMethyCall>& samples)
{
    namespace F = torch::nn::functional;
    data->forward(inputs, prob_tensor_cpu, ctx);
    prob_tensor_cpu = prob_tensor_cpu.contiguous();
    at::Tensor prob = F::softmax(prob_tensor_cpu, F::SoftmaxFuncOptions(1));
    auto prob_a = prob.accessor<float, 2>();
    for (int i = 0; i < num_samples; ++i) {
        int v = prob_a[i][1] * 255;
        if (v > 255) v = 255;
        samples[i].scaled_prob = v;
    }
}

static void*
s_work_thread(void* params)
{
    ThreadWorkData* data = (ThreadWorkData*)(params);
    ModOptions* options = data->M_options;
    std::set<int> skipped_bam_tags; if (!options->keep_kinetics) fill_skipped_tags(skipped_bam_tags);
    int thread_id = data->thread_id();
    EvalKmerFeaturesGenerator ekfg(data->M_kmer_size);
    auto tensor_options = torch::TensorOptions().device(torch::kCPU).dtype(torch::kFloat32).requires_grad(false);
    std::vector<std::pair<int, bam1_t*>> processed_bams;
    std::vector<MolMethyCall> samples, batch_cpg_samples, batch_chg_samples, batch_chh_samples;
    size_t cpg_samples = 0, chg_samples = 0, chh_samples = 0;
    int added_cpg_samples = 0, added_chg_samples = 0, added_chh_samples = 0;

    torch::Tensor cpg_cpu = torch::empty(
        {options->sample_batch_size, data->M_kmer_size, EvalKmerFeaturesGenerator::EVAL_BASE_FEATURES},
        torch::dtype(torch::kFloat32).device(torch::kCPU).pinned_memory(true));
    float* Fcpg = cpg_cpu.data_ptr<float>();
    float* pFcpg = Fcpg;

    torch::Tensor chg_cpu = torch::empty(
        {options->sample_batch_size, data->M_kmer_size, EvalKmerFeaturesGenerator::EVAL_BASE_FEATURES},
        torch::dtype(torch::kFloat32).device(torch::kCPU).pinned_memory(true));
    float* Fchg = chg_cpu.data_ptr<float>();
    float* pFchg = Fchg;

    torch::Tensor chh_cpu = torch::empty(
        {options->sample_batch_size, data->M_kmer_size, EvalKmerFeaturesGenerator::EVAL_BASE_FEATURES},
        torch::dtype(torch::kFloat32).device(torch::kCPU).pinned_memory(true));
    float* Fchh = chh_cpu.data_ptr<float>();
    float* pFchh = Fchh;

    torch::Tensor F_gpu = torch::empty(
        {options->sample_batch_size, data->M_kmer_size, EvalKmerFeaturesGenerator::EVAL_BASE_FEATURES},
        torch::dtype(torch::kFloat32).device(torch::kCUDA));
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(F_gpu);

    torch::Tensor P_cpu = torch::empty({options->sample_batch_size, 2},
        torch::dtype(torch::kFloat32).device(torch::kCPU).pinned_memory(true));

    torch::NoGradGuard no_grad;

    while (1) {
        int read_id;
        bam1_t* bam = bam_init1();
        if (!data->get_next_sam(read_id, bam)) {
            bam_destroy1(bam);
            break;
        }
        //HBN_LOG("mod call bam %d:%s:%d", read_id, bam_get_qname(bam), bam->core.l_qseq);
        if ((read_id % 1000) == 0) HBN_LOG("%d\t%10d reads processed", thread_id, read_id);
        if (bam->core.l_qseq < options->min_read_size) {
            processed_bams.emplace_back(read_id, bam);
            continue;
        }
        if (!ekfg.init(bam)) {
            processed_bams.emplace_back(read_id, bam);
            continue;
        }

        if (!options->disable_cpg) {
            ekfg.extract_cpg_samples();
            int qoff, strand;
            MolMethyCall M; M.qid = read_id;
            while (ekfg.get_next_sample_features(pFcpg, qoff, strand)) {
                pFcpg += ekfg.M_features_per_sample;
                M.qoff = qoff;
                M.strand = strand;
                batch_cpg_samples.push_back(M);
                ++added_cpg_samples;
                if (added_cpg_samples == options->sample_batch_size) {
                    F_gpu.copy_(cpg_cpu, true);
                    s_call_one_batch_samples(data, tensor_options, EMethyContext::eCpG, 
                        inputs, P_cpu, data->M_kmer_size, added_cpg_samples, batch_cpg_samples);
                    samples.insert(samples.end(), batch_cpg_samples.begin(), batch_cpg_samples.end());
                    cpg_samples += added_cpg_samples;
                    batch_cpg_samples.clear();
                    added_cpg_samples = 0;
                    pFcpg = Fcpg;
                }
            }
        }

        if (!options->disable_chg) {
            ekfg.extract_chg_samples();
            int qoff, strand;
            MolMethyCall M; M.qid = read_id;
            while (ekfg.get_next_sample_features(pFchg, qoff, strand)) {
                pFchg += ekfg.M_features_per_sample;
                M.qoff = qoff;
                M.strand = strand;
                batch_chg_samples.push_back(M);
                ++added_chg_samples;
                if (added_chg_samples == options->sample_batch_size) {
                    F_gpu.copy_(chg_cpu, true);
                    s_call_one_batch_samples(data, tensor_options, EMethyContext::eCHG, 
                        inputs, P_cpu, data->M_kmer_size, added_chg_samples, batch_chg_samples);
                    samples.insert(samples.end(), batch_chg_samples.begin(), batch_chg_samples.end());
                    chg_samples += added_chg_samples;
                    batch_chg_samples.clear();
                    added_chg_samples = 0;
                    pFchg = Fchg;
                }
            }
        }

        if (!options->disable_chh) {
            ekfg.extract_chh_samples();
            int qoff, strand;
            MolMethyCall M; M.qid = read_id;
            while (ekfg.get_next_sample_features(pFchh, qoff, strand)) {
                pFchh += ekfg.M_features_per_sample;
                M.qoff = qoff;
                M.strand = strand;
                batch_chh_samples.push_back(M);
                ++added_chh_samples;
                if (added_chh_samples == options->sample_batch_size) {
                    F_gpu.copy_(chh_cpu, true);
                    s_call_one_batch_samples(data, tensor_options, EMethyContext::eCHH, 
                        inputs, P_cpu, data->M_kmer_size, added_chh_samples, batch_chh_samples);
                    samples.insert(samples.end(), batch_chh_samples.begin(), batch_chh_samples.end());
                    chh_samples += added_chh_samples;
                    batch_chh_samples.clear();
                    added_chh_samples = 0;
                    pFchh = Fchh;
                }
            }
        }

        processed_bams.emplace_back(read_id, bam);
    }
    if (added_cpg_samples) {
        F_gpu.copy_(cpg_cpu, true);
        s_call_one_batch_samples(data, tensor_options, EMethyContext::eCpG,
            inputs, P_cpu, data->M_kmer_size, added_cpg_samples, batch_cpg_samples);
        samples.insert(samples.end(), batch_cpg_samples.begin(), batch_cpg_samples.end());
        cpg_samples += added_cpg_samples;
        batch_cpg_samples.clear();
        added_cpg_samples = 0;
    }
    if (added_chg_samples) {
        F_gpu.copy_(chg_cpu, true);
        s_call_one_batch_samples(data, tensor_options, EMethyContext::eCHG,
            inputs, P_cpu, data->M_kmer_size, added_chg_samples, batch_chg_samples);
        samples.insert(samples.end(), batch_chg_samples.begin(), batch_chg_samples.end());
        chg_samples += added_chg_samples;
        batch_chg_samples.clear();
        added_chg_samples = 0;
    }
    if (added_chh_samples) {
        F_gpu.copy_(chh_cpu, true);
        s_call_one_batch_samples(data, tensor_options, EMethyContext::eCHH,
            inputs, P_cpu, data->M_kmer_size, added_chh_samples, batch_chh_samples);
        samples.insert(samples.end(), batch_chh_samples.begin(), batch_chh_samples.end());
        chh_samples += added_chh_samples;
        batch_chh_samples.clear();
        added_chh_samples = 0;
    }

    MolMethyCall* SA = samples.data();
    size_t SC = samples.size();
    pdqsort(SA, SA + SC, [](const MolMethyCall& x, const MolMethyCall& y) { return x.qid < y.qid; });
    std::map<int, std::pair<size_t, size_t>> bam_samples;
    size_t i = 0;
    while (i < SC) {
        size_t j = i + 1;
        while (j < SC && SA[i].qid == SA[j].qid) ++j;
        bam_samples[SA[i].qid] = std::pair<size_t, size_t>(i, j);
        i = j;
    }
    std::vector<MolMethyCall> fwd_calls, rev_calls;
    for (auto& bi : processed_bams) {
        int read_id = bi.first;
        bam1_t* bam = bi.second;
        auto iter = bam_samples.find(read_id);
        fwd_calls.clear();
        rev_calls.clear();
        if (iter != bam_samples.end()) {
            size_t s = iter->second.first;
            size_t e = iter->second.second;
            for (size_t k = s; k < e; ++k) {
                if (SA[k].strand == FWD) {
                    fwd_calls.push_back(SA[k]);
                } else {
                    rev_calls.push_back(SA[k]);
                }
            }
            pdqsort(fwd_calls.begin(), fwd_calls.end(), [](const MolMethyCall& x, const MolMethyCall& y) {
                return x.qoff < y.qoff;
            });
            pdqsort(rev_calls.begin(), rev_calls.end(), [](const MolMethyCall& x, const MolMethyCall& y) {
                return x.qoff < y.qoff;
            });            
        }
        build_one_mod_bam(bam, skipped_bam_tags, fwd_calls.data(), fwd_calls.size(), rev_calls.data(), rev_calls.size());
    }
    data->add_mod_bams(processed_bams.data(), processed_bams.size(), cpg_samples, chg_samples, chh_samples);

    return nullptr;
}

extern "C"
size_t getMemorySizeBytes();
static void dump_program_info(FILE* out)
{
    struct utsname _os_info_buf;
    struct utsname* os_info = nullptr;
    if (uname(&_os_info_buf) == 0) os_info = &_os_info_buf;

    size_t sys_mem_bytes = getMemorySizeBytes();
    std::string sys_mem = bytes_to_datasize(sys_mem_bytes);

    int cpu_threads = std::thread::hardware_concurrency();

    fprintf(out, "\n");
    fprintf(out, "PROGRAM:\n");
    fprintf(out, "  Name:           %s\n", HBN_PACKAGE_NAME);
    fprintf(out, "  Version:        %s\n", HBN_PACKAGE_VERSION);
    fprintf(out, "  htslib:         %s\n", hts_version());
    fprintf(out, "  Description:    Methylation states detection toolkit for HiFi reads\n");
    fprintf(out, "  Contact:        chenying2016@gmail.com\n");

    fprintf(out, "\n");
    fprintf(out, "SYSTEM:\n");
    if (os_info) {
    fprintf(out, "  Computer:       %s\n", os_info->nodename);
    fprintf(out, "  Name:           %s\n", os_info->sysname);
    fprintf(out, "  Release:        %s\n", os_info->release);
    fprintf(out, "  Version:        %s\n", os_info->version);
    fprintf(out, "  Machine:        %s\n", os_info->machine);
    }
    fprintf(out, "  CPU threads:    %d\n", cpu_threads);
    fprintf(out, "  RAM:            %s\n", sys_mem.c_str());
    fprintf(out, "\n");
}

int s_single_molecular_methy_call(int argc, char* argv[])
{
    ModOptions options;
    if (!options.parse(argc, argv)) {
        options.dump_usage(argc, argv);
        return EXIT_FAILURE;
    }
    HbnProgramInfo hpi(HBN_PACKAGE_NAME, dump_program_info);
    options.dump_parameters();
    fprintf(stderr, "\n\n");

    //at::set_num_threads(8);
    //at::set_num_interop_threads(8);

    ThreadWorkData data(&options, argc, argv);
    const int num_threads = options.num_threads;
    pthread_t jobs[num_threads];

    while (data.reset_batch_read_idx()) {
        for (int i = 0; i < num_threads; ++i) {
            pthread_create(jobs + i, NULL, s_work_thread, &data);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(jobs[i], NULL);
        }
        data.save_mod_bam();
    }

    return 0;
}

} // ns_5mc_call

int single_molecular_methy_call(int argc, char* argv[])
{
    return ns_5mc_call::s_single_molecular_methy_call(argc, argv);
}
