#include "../../corelib/5mc_context.hpp"
#include "../../corelib/5mc_motif_finder.hpp"
#include "../../corelib/arg_parse.hpp"
#include "../../corelib/bam_info.hpp"
#include "../../corelib/bam_mod_parser.hpp"
#include "../../corelib/line_reader.hpp"
#include "../../corelib/pdqsort.h"
#include "../../corelib/sam_batch.hpp"

#include <algorithm>
#include <mutex>
#include <sstream>
#include <vector>

namespace ns_read_level_eval_samples {

static constexpr size_t kTargetSamples = 100000;

struct ReadLevelEvalBismarkSite
{
    size_t* offsets;
    char* labels;

    ReadLevelEvalBismarkSite(sam_hdr_t* hdr) {
        const int num_chr = sam_hdr_nref(hdr);
        offsets = new size_t[num_chr];
        size_t num_bases = 0;
        for (int i = 0; i < num_chr; ++i) {
            offsets[i] = num_bases;
            num_bases += sam_hdr_tid2len(hdr, i);
        }
        labels = new char[num_bases];
        std::fill(labels, labels + num_bases, -1);
    }

    ~ReadLevelEvalBismarkSite() {
        delete[] offsets;
        delete[] labels;
    }
};

static void
s_fill_chr_base_label_with_bismark(sam_hdr_t* hdr, const char* bismark_bed_path, ReadLevelEvalBismarkSite* bismark)
{
    std::string last_chr, chr, line;
    int last_sid = -1;
    size_t np = 0, nn = 0;
    HbnLineReader in(bismark_bed_path);
    while (!in.AtEOF()) {
        ++in;
        auto vline = *in;
        line.assign(vline.first, vline.second);
        const char* p = line.c_str();
        int pl = line.size();
        int i = 0;
        int j = 0;

        // chr
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        chr.assign(p+i, j - i);
        if (chr != last_chr) {
            last_chr = chr;
            last_sid = sam_hdr_name2tid(hdr, last_chr.c_str());
        }
        char* chr_labels = bismark->labels + bismark->offsets[last_sid];

        // soff
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int soff = atoi(p + i);

        /// send
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int send = atoi(p + i);
        hbn_assert(send - soff == 1, "%s", line.c_str());

        /// freq
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        //double freq = atof(p + i); 

        /// pos-cov
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int pcov = atoi(p + i);       

        /// neg-cov
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int ncov = atoi(p + i);  

        if (pcov + ncov < 10) continue;
        if (pcov == 0) {
            chr_labels[soff] = 0;
            ++nn;
        } else if (ncov == 0) {
            chr_labels[soff] = 1;
            ++np;
        }
    }
    HBN_LOG("Load %zu methylated sites and %zu unmethylated sites from %s", np, nn, bismark_bed_path);
}

/////////////////////////

class ProbThresholdThreadWorkData
{
public:
    ProbThresholdThreadWorkData(const char* mod_bam_path) {
        M_sam = new SAM_Batch(mod_bam_path, 10000);
        std::fill(M_cpg_bins, M_cpg_bins + 256, 0);
        std::fill(M_chg_bins, M_chg_bins + 256, 0);
        std::fill(M_chh_bins, M_chh_bins + 256, 0);
    }

    ~ProbThresholdThreadWorkData() {
        delete M_sam;
    }

    bool get_next_bam(int& read_id, bam1_t* bam) {
        return M_sam->get_next_sam(read_id, bam);
    }

    void add_counts(size_t* a, size_t* b, size_t* c) {
        std::lock_guard<std::mutex> __(M_mutex);
        for (int i = 0; i < 256; ++i) {
            M_cpg_bins[i] += a[i];
            M_chg_bins[i] += b[i];
            M_chh_bins[i] += c[i];
        }
    }

public:
    SAM_Batch*      M_sam;
    size_t          M_cpg_bins[256];
    size_t          M_chg_bins[256];
    size_t          M_chh_bins[256];
    std::mutex      M_mutex;
};

static void*
s_prob_bin_thread(void* params)
{
    ProbThresholdThreadWorkData* data = (ProbThresholdThreadWorkData*)(params);
    bam1_t* bam = bam_init1();
    BamQuerySequence query;
    MethylationContext methy_ctx;
    std::vector<BaseModInfo> base_mods;
    size_t cpg_bins[256];
    size_t chg_bins[256];
    size_t chh_bins[256];
    std::fill(cpg_bins, cpg_bins + 256, 0);
    std::fill(chg_bins, chg_bins + 256, 0);
    std::fill(chh_bins, chh_bins + 256, 0);
    int read_id;
    while (data->get_next_bam(read_id, bam)) {
        if (bam->core.flag & 0x900) continue;
        base_mods.clear();
        extract_bam_base_mods(bam, base_mods);
        if (base_mods.empty()) continue;
        if (!query.init(bam)) continue;
        for (auto& m : base_mods) {
            hbn_assert(query.fwd_rqs[m.qoff] == 'C' || query.fwd_rqs[m.qoff] == 'G');
            bool r;
            if (query.fwd_rqs[m.qoff] == 'C') {
                r = (m.qoff + 1 < query.size && query.fwd_rqs[m.qoff+1] == 'G');
                if (r) {
                    ++cpg_bins[m.scaled_prob];
                    continue;
                }
                r = (m.qoff + 2 < query.size && 
                    methy_ctx.get_fwd_chg_motif_idx(query.fwd_rqs + m.qoff) != methy_ctx.CHG_INVALID_MOTIF_IDX);
                if (r) {
                    ++chg_bins[m.scaled_prob];
                    continue;
                }
                r = (m.qoff + 2 < query.size &&
                    methy_ctx.get_fwd_chh_motif_idx(query.fwd_rqs + m.qoff) != methy_ctx.CHH_INVALID_MOTIF_IDX);
                if (r) {
                    ++chh_bins[m.scaled_prob];
                    continue;
                }
                continue;
            }
            {
                r = (m.qoff - 2 >= 0 &&
                    methy_ctx.get_rev_chh_motif_idx(query.fwd_rqs + m.qoff - 2) != methy_ctx.CHH_INVALID_MOTIF_IDX);
                if (r) {
                    ++chh_bins[m.scaled_prob];
                    continue;
                }
                continue;
            }
        }
    }
    data->add_counts(cpg_bins, chg_bins, chh_bins);

    return nullptr;
}

static void
s_resolve_scaled_prob_threshold(const char* mod_bam_path, const int num_threads,
    u8& cpg_threshold, u8& chg_threshold, u8& chh_threshold)
{
    ProbThresholdThreadWorkData data(mod_bam_path);
    pthread_t jobs[num_threads];
    while (data.M_sam->reset_batch_read_idx()) {
        for (int i = 0; i < num_threads; ++i) {
            pthread_create(jobs + i, nullptr, s_prob_bin_thread, &data);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(jobs[i], nullptr);
        }
    }

    const size_t* a = data.M_cpg_bins;
    size_t sum = 0;
    int min_i = -1;
    size_t min_cnt = std::numeric_limits<size_t>::max();
    int st = 20, en = 256-20;
    while (st < 256 && a[st] < 10) ++st;
    while (en && a[en-1] < 10) --en;
    if (en - st >= 50) {
        for (int i = st; i < en; ++i) {
            sum += a[i];
            if (min_cnt > a[i]) {
                min_cnt = a[i];
                min_i = i;
            }
        }
    }
    fprintf(stderr, "CpG samples: %zu\n", sum);
    if (sum < 10000 || min_i == -1) {
        fprintf(stderr, "Not enough samples for inferring scaled probability threshold, set it to 128\n");
        cpg_threshold = 128;
    } else {
        //fprintf(stderr, "st = %d, en = %d\n", st, en);
        fprintf(stderr, "CpG scaled probability threshold: %d\n", min_i);
        cpg_threshold = min_i;
    }

    a = data.M_chg_bins;
    sum = 0;
    min_i = -1;
    min_cnt = std::numeric_limits<size_t>::max();
    st = 20, en = 256-20;
    while (st < 256 && a[st] < 10) ++st;
    while (en && a[en-1] < 10) --en;
    if (en - st >= 50) {
        for (int i = st; i < en; ++i) {
            sum += a[i];
            if (min_cnt > a[i]) {
                min_cnt = a[i];
                min_i = i;
            }
        }
    }
    fprintf(stderr, "CHG samples: %zu\n", sum);
    if (sum < 10000 || min_i == -1) {
        fprintf(stderr, "Not enough samples for inferring scaled probability threshold, set it to 128\n");
        chg_threshold = 128;
    } else {
        //fprintf(stderr, "st = %d, en = %d\n", st, en);
        fprintf(stderr, "CHG scaled probability threshold: %d\n", min_i);
        chg_threshold = min_i;
    }

    a = data.M_chh_bins;
    sum = 0;
    min_i = -1;
    min_cnt = std::numeric_limits<size_t>::max();
    st = 20, en = 256-20;
    while (st < 256 && a[st] < 10) ++st;
    while (en && a[en-1] < 10) --en;
    if (en - st >= 50) {
        for (int i = st; i < en; ++i) {
            sum += a[i];
            if (min_cnt > a[i]) {
                min_cnt = a[i];
                min_i = i;
            }
        }
    }
    fprintf(stderr, "CHH samples: %zu\n", sum);
    if (sum < 10000 || min_i == -1) {
        fprintf(stderr, "Not enough samples for inferring scaled probability threshold, set it to 128\n");
        chh_threshold = 128;
    } else {
        //fprintf(stderr, "st = %d, en = %d\n", st, en);
        fprintf(stderr, "CHH scaled probability threshold: %d\n", min_i);
        chh_threshold = min_i;
    }
}

/////////////////////////////////

class ThreadWorkData
{
public:
    ThreadWorkData(const char* reference_path, const char* mod_bam_path,
        const char* bismark_bed_path) 
    {
        M_updb = new HbnDatabase(reference_path);
        M_sam = new SAM_Batch(mod_bam_path, 10000);
        M_bismark = new ReadLevelEvalBismarkSite(M_sam->sam_hdr());
        s_fill_chr_base_label_with_bismark(M_sam->sam_hdr(), bismark_bed_path, M_bismark);
    }

    ~ThreadWorkData() {
        delete M_updb;
        delete M_sam;
        delete M_bismark;
    }

    bool get_next_bam(int& read_id, bam1_t* bam) {
        return M_sam->get_next_sam(read_id, bam);
    }

    void add_predict_samples(u8* a_cpg_bp, size_t n_cpg_bp,
        u8* a_cpg_bn, size_t n_cpg_bn,
        u8* a_chg_bp, size_t n_chg_bp,
        u8* a_chg_bn, size_t n_chg_bn,
        u8* a_chh_bp, size_t n_chh_bp,
        u8* a_chh_bn, size_t n_chh_bn) 
    {
        std::lock_guard<std::mutex> __(M_mutex);
        M_cpg_bp.insert(M_cpg_bp.end(), a_cpg_bp, a_cpg_bp + n_cpg_bp);
        M_cpg_bn.insert(M_cpg_bn.end(), a_cpg_bn, a_cpg_bn + n_cpg_bn);

        M_chg_bp.insert(M_chg_bp.end(), a_chg_bp, a_chg_bp + n_chg_bp);
        M_chg_bn.insert(M_chg_bn.end(), a_chg_bn, a_chg_bn + n_chg_bn);

        M_chh_bp.insert(M_chh_bp.end(), a_chh_bp, a_chh_bp + n_chh_bp);
        M_chh_bn.insert(M_chh_bn.end(), a_chh_bn, a_chh_bn + n_chh_bn);
    }

    void over_sampling_eval_samples()
    {
        std::vector<u8> ss;

        std::vector<u8>* tt = &M_cpg_bp;
        const char* ctx = "CpG";
        const char* label = "positive";
        if (tt->size() > 0 && tt->size() < kTargetSamples) {
            fprintf(stderr, "Original %s %s samples: %zu\n", ctx, label, tt->size());
            ss = *tt;
            size_t x = 2 * kTargetSamples / tt->size();
            x *= 2;
            tt->clear();
            for (size_t i = 0; i < x; ++i) {
                tt->insert(tt->end(), ss.begin(), ss.end());
            }
            fprintf(stderr, "Over-sampled %s %s samples: %zu\n", ctx, label, tt->size());
        }

        tt = &M_cpg_bn;
        ctx = "CpG";
        label = "negative";
        if (tt->size() > 0 && tt->size() < kTargetSamples) {
            fprintf(stderr, "Original %s %s samples: %zu\n", ctx, label, tt->size());
            ss = *tt;
            size_t x = 2 * kTargetSamples / tt->size();
            x *= 2;
            tt->clear();
            for (size_t i = 0; i < x; ++i) {
                tt->insert(tt->end(), ss.begin(), ss.end());
            }
            fprintf(stderr, "Over-sampled %s %s samples: %zu\n", ctx, label, tt->size());
        }

        tt = &M_chg_bp;
        ctx = "CHG";
        label = "positive";
        if (tt->size() > 0 && tt->size() < kTargetSamples) {
            fprintf(stderr, "Original %s %s samples: %zu\n", ctx, label, tt->size());
            ss = *tt;
            size_t x = 2 * kTargetSamples / tt->size();
            x *= 2;
            tt->clear();
            for (size_t i = 0; i < x; ++i) {
                tt->insert(tt->end(), ss.begin(), ss.end());
            }
            fprintf(stderr, "Over-sampled %s %s samples: %zu\n", ctx, label, tt->size());
        }

        tt = &M_chg_bn;
        ctx = "CHG";
        label = "negative";
        if (tt->size() > 0 && tt->size() < kTargetSamples) {
            fprintf(stderr, "Original %s %s samples: %zu\n", ctx, label, tt->size());
            ss = *tt;
            size_t x = 2 * kTargetSamples / tt->size();
            x *= 2;
            tt->clear();
            for (size_t i = 0; i < x; ++i) {
                tt->insert(tt->end(), ss.begin(), ss.end());
            }
            fprintf(stderr, "Over-sampled %s %s samples: %zu\n", ctx, label, tt->size());
        }

        tt = &M_chh_bp;
        ctx = "CHH";
        label = "positive";
        if (tt->size() > 0 && tt->size() < kTargetSamples) {
            fprintf(stderr, "Original %s %s samples: %zu\n", ctx, label, tt->size());
            ss = *tt;
            size_t x = 2 * kTargetSamples / tt->size();
            x *= 2;
            tt->clear();
            for (size_t i = 0; i < x; ++i) {
                tt->insert(tt->end(), ss.begin(), ss.end());
            }
            fprintf(stderr, "Over-sampled %s %s samples: %zu\n", ctx, label, tt->size());
        }

        tt = &M_chh_bn;
        ctx = "CHH";
        label = "negative";
        if (tt->size() > 0 && tt->size() < kTargetSamples) {
            fprintf(stderr, "Original %s %s samples: %zu\n", ctx, label, tt->size());
            ss = *tt;
            size_t x = 2 * kTargetSamples / tt->size();
            x *= 2;
            tt->clear();
            for (size_t i = 0; i < x; ++i) {
                tt->insert(tt->end(), ss.begin(), ss.end());
            }
            fprintf(stderr, "Over-sampled %s %s samples: %zu\n", ctx, label, tt->size());
        }
    }

public:
    HbnDatabase*                M_updb;
    SAM_Batch*                  M_sam;
    ReadLevelEvalBismarkSite*   M_bismark;

    std::mutex                  M_mutex;
    std::vector<u8>           M_cpg_bp;
    std::vector<u8>           M_cpg_bn;
    std::vector<u8>           M_chg_bp;
    std::vector<u8>           M_chg_bn;
    std::vector<u8>           M_chh_bp;
    std::vector<u8>           M_chh_bn;
};

struct ReadBaseModInfo
{
    bool has_prob_value;
    uint8_t scaled_prob;

    void init() {
        has_prob_value = false;
        scaled_prob = 0;
    }
};

static void*
s_read_level_sample_thread(void* params)
{
    ThreadWorkData* data = (ThreadWorkData*)(params);
    BamQuerySequence query;
    BamMapInfo align;
    MethylationContext methy_ctx;
    std::vector<u8> cpg_bp, cpg_bn, chg_bp, chg_bn, chh_bp, chh_bn;
    std::vector<BaseModInfo> base_mods;
    std::vector<ReadBaseModInfo> read_mods;
    std::vector<MotifMappedHitInfo> mapped_samples;
    bam1_t* bam = bam_init1();
    int read_id;
    RandomFloatNumberGenerator rng;

    while (data->get_next_bam(read_id, bam)) {
        base_mods.clear();
        extract_bam_base_mods(bam, base_mods);
        if (base_mods.empty()) continue;
        //if (bam->core.qual < data->M_options->min_mapQ) continue;
        if (!query.init(bam)) continue;
        if (!align.init(data->M_sam->sam_hdr(), bam, data->M_updb, &query)) continue;
        //if (align.pi < data->M_options->min_pi) continue;

        read_mods.resize(query.size);
        for (auto& r : read_mods) r.init();
        for (auto& m : base_mods) {
            if (m.code != 'm') continue;
	        int qoff = m.qoff;
            read_mods[qoff].has_prob_value = true;
            read_mods[qoff].scaled_prob = m.scaled_prob;
        }
        const char* chr_bismark = data->M_bismark->labels + data->M_bismark->offsets[bam->core.tid];

        const char* qas = align.qas;
        const char* sas = align.sas;
        const int as_size = align.as_size;
        const int* qpos = align.qas_pos;
        const int* spos = align.sas_pos;

        for (int i = 0; i <= as_size - 2; ++i) {
            int soff = spos[i];
            if (strncmp(qas + i, "CG", 2) || strncmp(sas + i, "CG", 2)) continue;
            int qoff = (align.qdir == FWD) ? qpos[i] : query.size - 1 - (qpos[i] + 1);
            hbn_assert(query.fwd_rqs[qoff] == 'C');
            if (read_mods[qoff].has_prob_value) {
                if (chr_bismark[soff] == -1) continue;
                if (chr_bismark[soff] == 0) {
                    cpg_bn.push_back(read_mods[qoff].scaled_prob);
                } else {
                    cpg_bp.push_back(read_mods[qoff].scaled_prob);
                }
            }            
        }

        for (int i = 0; i <= as_size - 3; ++i) {
            int soff = spos[i];
            int qoff = -1;
            if (align.qdir == FWD) {
                if (strncmp(qas + i, "CCG", 3) == 0 && strncmp(sas + i, "CCG", 3) ==0) {
                    qoff = qpos[i];
                } else if (strncmp(qas + i, "CAG", 3) == 0 && strncmp(sas + i, "CAG", 3) == 0) {
                    qoff = qpos[i];
                } else if (strncmp(qas + i, "CTG", 3) == 0 && strncmp(sas + i, "CTG", 3) == 0) {
                    qoff = qpos[i];
                }
            } else {
                if (strncmp(qas + i, "CGG", 3) == 0 && strncmp(sas + i, "CGG", 3) ==0) {
                    qoff = query.size - 1 - (qpos[i] + 2);
                } else if (strncmp(qas + i, "CAG", 3) == 0 && strncmp(sas + i, "CAG", 3) == 0) {
                    qoff = query.size - 1 - (qpos[i] + 2);
                } else if (strncmp(qas + i, "CTG", 3) == 0 && strncmp(sas + i, "CTG", 3) == 0) {
                    qoff = query.size - 1 - (qpos[i] + 2);
                }                
            }
            if (qoff == -1) continue;
            hbn_assert(query.fwd_rqs[qoff] == 'C');
            if (read_mods[qoff].has_prob_value) {
                if (chr_bismark[soff] == -1) continue;
                if (chr_bismark[soff] == 0) {
                    chg_bn.push_back(read_mods[qoff].scaled_prob);
                } else {
                    chg_bp.push_back(read_mods[qoff].scaled_prob);
                }
            }  
        }

        mapped_samples.clear();
        extract_chh_mapped_samples(data->M_updb, methy_ctx, query, align, mapped_samples);
        for (auto& m : mapped_samples) {
            if (read_mods[m.qoff].has_prob_value) {
                if (chr_bismark[m.soff] == -1) continue;
                if (chr_bismark[m.soff] == 0) {
                    if (rng() <= 0.1) chh_bn.push_back(read_mods[m.qoff].scaled_prob);
                } else {
                    chh_bp.push_back(read_mods[m.qoff].scaled_prob);
                }
            }
        }
    }

    data->add_predict_samples(cpg_bp.data(), cpg_bp.size(),
        cpg_bn.data(), cpg_bn.size(),
        chg_bp.data(), chg_bp.size(),
        chg_bn.data(), chg_bn.size(),
        chh_bp.data(), chh_bp.size(),
        chh_bn.data(), chh_bn.size());

    return nullptr;
}

static void
s_dump_samples(u8* abp, size_t cbp, u8* abn, size_t cbn, u8 threshold, const char* output_prefix, const char* ctx)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    char path[HBN_MAX_PATH_LEN];

    for (int i = 0; i < 5; ++i) {
        snprintf(path, HBN_MAX_PATH_LEN, "%s.%s.%d", output_prefix, ctx, i);
        hbn_dfopen(out, path, "w");

        std::shuffle(abp, abp + cbp, gen);
        std::shuffle(abp, abp + cbp, gen);
        std::shuffle(abp, abp + cbp, gen);
        for (size_t k = 0; k < kTargetSamples; ++k) {
            int predict = (abp[k] >= threshold) ? 1 : 0;
            double prob = 1.0 * abp[k] / 255;
            fprintf(out, "1\t%d\t%g\n", predict, prob);
        }

        std::shuffle(abn, abn + cbn, gen);
        std::shuffle(abn, abn + cbn, gen);
        std::shuffle(abn, abn + cbn, gen);
        for (size_t k = 0; k < kTargetSamples; ++k) {
            int predict = (abn[k] >= threshold) ? 1 : 0;
            double prob = 1.0 * abn[k] / 255;
            fprintf(out, "0\t%d\t%g\n", predict, prob);
        }        

        hbn_fclose(out);
    }
}

static int
s_extract_read_leval_eval_samples(int argc, char* argv[])
{
    if (argc != 6) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s reference bismark mod-bam output-prefix\n", argv[0], argv[1]);
        return 1;
    }
    const char* reference_path = argv[2];
    const char* bismark_bed_path = argv[3];
    const char* mod_bam_path = argv[4];
    const char* output_prefix = argv[5];
    const int num_threads = 16;

    u8 cpg_threshold = 0, chg_threshold = 0, chh_threshold = 0;
    s_resolve_scaled_prob_threshold(mod_bam_path, num_threads, 
        cpg_threshold, chg_threshold, chh_threshold);

    ThreadWorkData data(reference_path, mod_bam_path, bismark_bed_path);
    pthread_t jobs[num_threads];
    while (data.M_sam->reset_batch_read_idx()) {
        for (int i = 0; i < num_threads; ++i) {
            pthread_create(jobs + i, nullptr, s_read_level_sample_thread, &data);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(jobs[i], nullptr);
        }
    }
    data.over_sampling_eval_samples();

    u8* abp = data.M_cpg_bp.data();
    size_t cbp = data.M_cpg_bp.size();
    u8* abn = data.M_cpg_bn.data();
    size_t cbn = data.M_cpg_bn.size();
    if (cbp > 0 && cbn > 0) {
        fprintf(stderr, "CpG positive samples: %zu, negative samples: %zu\n", cbp, cbn);
        s_dump_samples(abp, cbp, abn, cbn, cpg_threshold, output_prefix, "CpG");
    }

    abp = data.M_chg_bp.data();
    cbp = data.M_chg_bp.size();
    abn = data.M_chg_bn.data();
    cbn = data.M_chg_bn.size();
    if (cbp > 0 && cbn > 0) {
        fprintf(stderr, "CHG positive samples: %zu, negative samples: %zu\n", cbp, cbn);
        s_dump_samples(abp, cbp, abn, cbn, chg_threshold, output_prefix, "CHG");
    }

    abp = data.M_chh_bp.data();
    cbp = data.M_chh_bp.size();
    abn = data.M_chh_bn.data();
    cbn = data.M_chh_bn.size();
    if (cbp > 0 && cbn > 0) {
        fprintf(stderr, "CHH positive samples: %zu, negative samples: %zu\n", cbp, cbn);
        s_dump_samples(abp, cbp, abn, cbn, chh_threshold, output_prefix, "CHH");
    }    

    return 0;
}

} // ns_read_level_eval_samples

int eval_main(int argc, char* argv[])
{
    return ns_read_level_eval_samples::s_extract_read_leval_eval_samples(argc, argv);
}
