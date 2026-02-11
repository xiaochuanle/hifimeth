#include "mod_batch.hpp"

#include <corelib/build_mod_bam.hpp>
#include <corelib/get_core_count.hpp>
#include <corelib/pdqsort.h>
#include <corelib/sam_batch.hpp>
#include <htslib/hts.h>
#include "program_info.hpp"

#include <thread>

#include <sys/stat.h>
#include <sys/utsname.h>
#include <dirent.h>

namespace ns_mods {

struct ModModels {
    ModParams mod_params;

    ov::Core cpg_core;
    std::shared_ptr<ov::Model> cpg_model;
    ov::CompiledModel cpg_compiled_model;
    ov::Core chg_core;
    std::shared_ptr<ov::Model> chg_model;
    ov::CompiledModel chg_compiled_model;

    ov::Core chh_core;
    std::shared_ptr<ov::Model> chh_model;
    ov::CompiledModel chh_compiled_model;

    void s_load_one_model(const char* model_path,
        ov::Core& core,
        std::shared_ptr<ov::Model>& model,
        ov::CompiledModel& compiled_model,
        int& kmer_size,
        int& features_per_base,
        int& batch_size)
    {
        model = core.read_model(model_path);
        const auto& input_port = model->input();
        const ov::PartialShape& p_shape = input_port.get_partial_shape();
        if (p_shape.rank().get_length() != 3) {
            fprintf(stderr, "Ill-formated IR model: %s\n", model_path);
            fprintf(stderr, "  Model input rank must be 3 (Batch, Kmer, Features)\n");
            abort();
        }
        if (p_shape[1].is_dynamic() || p_shape[2].is_dynamic()) {
            fprintf(stderr, "Ill-formated IR model: %s\n", model_path);
            fprintf(stderr, "  Model Kmer or Feature dimension is dynamic, which is not supported.\n");
            abort();
        }
        if (p_shape[0].is_static()) {
            batch_size = p_shape[0].get_length();
        }
        kmer_size = p_shape[1].get_length();
        features_per_base = p_shape[2].get_length();

        fprintf(stderr, "######## PARAMS for model %s:\n", model_path);
        fprintf(stderr, "  ## KMER: %d\n", kmer_size);
        fprintf(stderr, "  ## BASE_FEATURES: %d\n", features_per_base);
        fprintf(stderr, "  ## BATCH: %d\n", batch_size);

        model->reshape({batch_size, kmer_size, features_per_base });
        compiled_model = core.compile_model(model, "CPU", 
            ov::hint::performance_mode(ov::hint::PerformanceMode::THROUGHPUT));
    }

    ModModels(const ModOptions* options) {
        char path[HBN_MAX_PATH_LEN];

        mod_params.disable_cpg = true; 
        if (!options->disable_cpg) {
            mod_params.disable_cpg = false;
            mod_params.cpg_batch_size = options->sample_batch_size;
            snprintf(path, HBN_MAX_PATH_LEN, "%s/CpG.onnx", options->model_dir.c_str());
            s_load_one_model(path, cpg_core, cpg_model, cpg_compiled_model, 
                mod_params.cpg_kmer_size, mod_params.cpg_features_per_base, mod_params.cpg_batch_size);
        }

        mod_params.disable_chg = true;
        if (!options->disable_chg) {
            mod_params.disable_chg = false;
            mod_params.chg_batch_size = options->sample_batch_size;
            snprintf(path, HBN_MAX_PATH_LEN, "%s/CHG.onnx", options->model_dir.c_str());
            s_load_one_model(path, chg_core, chg_model, chg_compiled_model,
                mod_params.chg_kmer_size, mod_params.chg_features_per_base, mod_params.chg_batch_size);
        }

        mod_params.disable_chh = true;
        if (!options->disable_chh) {
            mod_params.disable_chh = false;
            mod_params.chh_batch_size = options->sample_batch_size;
            snprintf(path, HBN_MAX_PATH_LEN, "%s/CHH.onnx", options->model_dir.c_str());
            s_load_one_model(path, chh_core, chh_model, chh_compiled_model,
                mod_params.chh_kmer_size, mod_params.chh_features_per_base, mod_params.chh_batch_size);
        }
    }
};

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
s_worker_thread(int thread_id,
    SAM_Batch* sam_reader,
    ModOptions* options,
    ModModels* models,
    std::vector<std::pair<size_t, bam1_t*>>* all_mod_bams,
    size_t* processed_cpg_samples,
    size_t* processed_chg_samples,
    size_t* processed_chh_samples,
    std::mutex* out_lock)
{
    EvalKmerFeaturesGenerator ekfg;
    std::vector<MolMethyCall> samples;
    ModBatch* cpg_caller = nullptr;
    if (!options->disable_cpg) cpg_caller = new ModBatch(models->mod_params.cpg_kmer_size,
            models->mod_params.cpg_features_per_base,
            options->sample_batch_size,
            &samples,
            models->cpg_compiled_model);
    ModBatch* chg_caller = nullptr;
    if (!options->disable_chg) chg_caller = new ModBatch(models->mod_params.chg_kmer_size,
        models->mod_params.chg_features_per_base,
        options->sample_batch_size,
        &samples,
        models->chg_compiled_model);
    ModBatch* chh_caller = nullptr;
    if (!options->disable_chh) chh_caller = new ModBatch(models->mod_params.chh_kmer_size,
        models->mod_params.chh_features_per_base,
        options->sample_batch_size,
        &samples,
        models->chh_compiled_model);

    std::set<int> skipped_bam_tags; if (!options->keep_kinetics) fill_skipped_tags(skipped_bam_tags);
    std::vector<std::pair<size_t, bam1_t*>> mod_bams;

    while (1) {
        int read_id;
        bam1_t* bam = bam_init1();
        if (!sam_reader->get_next_sam(read_id, bam)) {
            bam_destroy1(bam);
            break;
        }
        //HBN_LOG("mod call bam %d:%s:%d", read_id, bam_get_qname(bam), bam->core.l_qseq);
        if (((read_id+1) % 1000) == 0) HBN_LOG("%d\t%10d reads processed", thread_id, read_id+1);
        if (bam->core.l_qseq < options->min_read_size) {
            mod_bams.emplace_back(read_id, bam);
            continue;
        }
        if (!ekfg.init(bam)) {
            mod_bams.emplace_back(read_id, bam);
            continue;
        }

        if (cpg_caller) {
            ekfg.extract_cpg_samples();
            cpg_caller->call_mods_for_one_read(read_id, &ekfg);
        }
        if (chg_caller) {
            ekfg.extract_chg_samples();
            chg_caller->call_mods_for_one_read(read_id, &ekfg);
        }
        if (chh_caller) {
            ekfg.extract_chh_samples();
            chh_caller->call_mods_for_one_read(read_id, &ekfg);
        }

        mod_bams.emplace_back(read_id, bam);
    }
    if (cpg_caller && cpg_caller->num_samples) cpg_caller->call_current_batch();
    if (chg_caller && chg_caller->num_samples) chg_caller->call_current_batch();
    if (chh_caller && chh_caller->num_samples) chh_caller->call_current_batch();

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
    for (auto& bi : mod_bams) {
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

    {
        std::lock_guard<std::mutex> __(*out_lock);
        all_mod_bams->insert(all_mod_bams->end(), mod_bams.begin(), mod_bams.end());
        if (cpg_caller) *processed_cpg_samples += cpg_caller->total_samples;
        if (chg_caller) *processed_chg_samples += chg_caller->total_samples;
        if (chh_caller) *processed_chh_samples += chh_caller->total_samples;
    }
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

    int physical_cpu_threads = get_physical_core_count();
    int logical_cpu_threads = std::thread::hardware_concurrency();
    ov::Version ov_version = ov::get_openvino_version();

    fprintf(out, "\n");
    fprintf(out, "PROGRAM:\n");
    fprintf(out, "  Name:                   %s\n", HBN_PACKAGE_NAME);
    fprintf(out, "  Version:                %s\n", HBN_PACKAGE_VERSION);
    fprintf(out, "  htslib:                 %s\n", hts_version());
    fprintf(out, "  OpenVINO:               %s\n", ov_version.buildNumber);
    fprintf(out, "  Description:            Methylation states detection toolkit for HiFi reads\n");
    fprintf(out, "  Contact:                chenying2016@gmail.com\n");

    fprintf(out, "\n");
    fprintf(out, "SYSTEM:\n");
    if (os_info) {
    fprintf(out, "  Computer:                %s\n", os_info->nodename);
    fprintf(out, "  Name:                    %s\n", os_info->sysname);
    fprintf(out, "  Release:                 %s\n", os_info->release);
    fprintf(out, "  Version:                 %s\n", os_info->version);
    fprintf(out, "  Machine:                 %s\n", os_info->machine);
    }
    fprintf(out, "  Physical CPU threads:    %d\n", physical_cpu_threads);
    fprintf(out, "  Logical CPU threads:     %d\n", logical_cpu_threads);
    fprintf(out, "  RAM:                     %s\n", sys_mem.c_str());
    fprintf(out, "\n");
}

int mod_main(int argc, char* argv[])
{
    ModOptions options;
    if (!options.parse(argc, argv)) {
        options.dump_usage(argc, argv);
        return EXIT_FAILURE;
    }
    HbnProgramInfo hpi(HBN_PACKAGE_NAME, dump_program_info);
    options.dump_parameters();

    SAM_Batch* M_sam = new SAM_Batch(options.bam_path.c_str(), options.read_batch_size);
    sam_hdr_t* hdr_in = M_sam->sam_hdr();
    sam_hdr_t* hdr_out = sam_hdr_dup(hdr_in);
    samFile* M_out = sam_open(options.mod_bam_path.c_str(), "wb");
	hts_set_threads(M_out, 8);
	hts_set_cache_size(M_out, 100000000);
    add_cmd_to_sam_hdr(argc, argv, hdr_out);
    if (sam_hdr_write(M_out, hdr_out)) HBN_ERR("Could not write BAM header to %s", options.mod_bam_path);
    sam_hdr_destroy(hdr_out);
    size_t all_cpg = 0, all_chg = 0, all_chh = 0;
    size_t all_reads = 0, all_bases = 0;
    std::mutex out_lock;
    std::vector<std::pair<size_t, bam1_t*>> mod_bams;
    ModModels models(&options);

    fprintf(stderr, "\n\n");
    size_t batch_idx = 0;
    while (M_sam->reset_batch_read_idx()) {
        HBN_LOG("Batch %zu", batch_idx++);
        std::vector<std::thread> threads;
        size_t batch_cpg = 0, batch_chg = 0, batch_chh = 0;
        for (int i = 0; i < options.num_threads; ++i) {
                threads.emplace_back(
                    s_worker_thread,
                    i,
                    M_sam,
                    &options,
                    &models,
                    &mod_bams,
                    &batch_cpg,
                    &batch_chg,
                    &batch_chh,
                    &out_lock
                );
        }
        for (auto& t : threads) {
            if (t.joinable()) t.join();
        }

        HBN_LOG("Save %zu MOD BAM", mod_bams.size());
        pdqsort(mod_bams.begin(), mod_bams.end(), 
            [](const std::pair<int, bam1_t*>& x, const std::pair<int, bam1_t*>& y) { return x.first < y.first; });
        size_t batch_bases = 0;
        for (auto& rdmod : mod_bams) {
            batch_bases += rdmod.second->core.l_qseq;
            if (sam_write1(M_out, hdr_out, rdmod.second) < 0) {
                HBN_ERR("Could not write BAM record to file");
            }
            bam_destroy1(rdmod.second);
        }
        HBN_LOG("Done.");
        fprintf(stderr, "######## Batch stats:\n");
        fprintf(stderr, "  ## Reads: %zu\n", mod_bams.size());
        std::string size = bytes_to_datasize(batch_bases);
        fprintf(stderr, "  ## Bases: %s\n", size.c_str());
        if (batch_cpg) {
            size = format_i64_with_commas(batch_cpg);
            fprintf(stderr, "  ## CpG samples: %s\n", size.c_str());
        }
        if (batch_chg) {
            size = format_i64_with_commas(batch_chg);
            fprintf(stderr, "  ## CHG samples: %s\n", size.c_str());
        }
        if (batch_chh) {
            size = format_i64_with_commas(batch_chh);
            fprintf(stderr, "  ## CHH samples: %s\n", size.c_str());
        }
        
        all_reads += mod_bams.size();
        all_bases += batch_bases;
        all_cpg += batch_cpg;
        all_chg += batch_chg;
        all_chh += batch_chh;

        mod_bams.clear();
    }

    {
        fprintf(stderr, "******** Final stats:\n");
        fprintf(stderr, "  ## Reads: %zu\n", all_reads);
        std::string size = bytes_to_datasize(all_bases);
        fprintf(stderr, "  ## Bases: %s\n", size.c_str());
        if (all_cpg) {
            size = format_i64_with_commas(all_cpg);
            fprintf(stderr, "  ## CpG samples: %s\n", size.c_str());
        }
        if (all_chg) {
            size = format_i64_with_commas(all_chg);
            fprintf(stderr, "  ## CHG samples: %s\n", size.c_str());
        }
        if (all_chh) {
            size = format_i64_with_commas(all_chh);
            fprintf(stderr, "  ## CHH samples: %s\n", size.c_str());
        }
    }

    delete M_sam;
    sam_close(M_out);
    return 0;
}

} // ns_mods

int mod_main(int argc, char* argv[])
{
    return ns_mods::mod_main(argc, argv);
}
