#include "../../corelib/5mc_context.hpp"
#include "../../corelib/5mc_motif_finder.hpp"
#include "../../corelib/arg_parse.hpp"
#include "../../corelib/bam_info.hpp"
#include "../../corelib/bam_mod_parser.hpp"
#include "../../corelib/pdqsort.h"
#include "../../corelib/sam_batch.hpp"
#include "program_info.hpp"

#include <mutex>
#include <sstream>
#include <vector>

#include <cstdio>

namespace ns_pileup {

static constexpr int kMinMapQ = 0;
static constexpr double kMinPi = 0.0;
static constexpr int kNumThreads = 8;

struct HbnAggrMethyCallOptions 
{
    int min_mapQ { kMinMapQ };
    double min_pi { kMinPi };
    int num_threads { kNumThreads };

    const char* reference_path { nullptr };
    const char* mod_bam_path { nullptr };
    const char* output_prefix { nullptr };

    bool parse(int argc, char* argv[]) {
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "-v") == 0) {
                fprintf(stderr, "%s\n", HBN_PACKAGE_VERSION);
                exit (0);
            }
            if (strcmp(argv[i], "-h") == 0) {
                dump_usage(argc, argv);
                exit (0);
            }
        }

        int i = 2;
        while (i < argc) {
            if (argv[i][0] != '-') break;
            if (argv[i][0] == '-' && strlen(argv[i]) == 1) break;

            if (parse_int_arg_value(argc, argv, i, "-q", min_mapQ)) continue;
            if (parse_real_arg_value(argc, argv, i, "-f", min_pi)) continue;
            if (parse_int_arg_value(argc, argv, i, "-t", num_threads)) continue;

            fprintf(stderr, "ERROR: unrecognised option %s", argv[i]);
            return false;
        }

        if (i >= argc) return false;
        reference_path = argv[i];
        ++i;

        if (i >= argc) return false;
        mod_bam_path = argv[i];
        ++i;

        if (i >= argc) return false;
        output_prefix = argv[i];
        ++i;

        if (i != argc) return false;
        return true;
    }

    void dump_usage(int argc, char* argv[]) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "  %s %s [OPTIONS] reference mod-bam output-prefix\n", argv[0], argv[1]);

        fprintf(stderr, "\n\n");
        fprintf(stderr, "DESCRIPTION:\n");
        fprintf(stderr, "  Compute aggregate cytosine methylation states on the genomic reference\n");

        fprintf(stderr, "\n\n");
        fprintf(stderr, "VERSION:\n");
        fprintf(stderr, "  %s\n", HBN_PACKAGE_VERSION);

        fprintf(stderr, "\n\n");
        fprintf(stderr, "OPTIONAL ARGUMENTS:\n");
        fprintf(stderr, "  -v\n");
        fprintf(stderr, "    Print version info and exit\n");
        fprintf(stderr, "  -h\n");
        fprintf(stderr, "    Print this help info and exit\n");
        fprintf(stderr, "  -q <mapQ>\n");
        fprintf(stderr, "    Minimum mapping quality score\n");
        fprintf(stderr, "    Default: %d\n", kMinMapQ);
        fprintf(stderr, "  -f <Alignment identity>\n");
        fprintf(stderr, "    Default: %g\n", kMinPi);
        fprintf(stderr, "  -t <CPU threads>\n");
        fprintf(stderr, "    Number of CPU threads\n");
        fprintf(stderr, "    Default: %d\n", kNumThreads);
    }

    void dump_parameters() {
        fprintf(stderr, "\n\n");
        fprintf(stderr, "====================> Parameters:\n");
        fprintf(stderr, "min-mapQ: %d\n", min_mapQ);
        fprintf(stderr, "min-identity: %g\n", min_pi);
        fprintf(stderr, "CPU threads: %d\n", num_threads);
        fprintf(stderr, "Genomic reference: %s\n", reference_path);
        fprintf(stderr, "mod-bam: %s\n", mod_bam_path);
        fprintf(stderr, "output prefix: %s\n", output_prefix);
        fprintf(stderr, "\n\n");
    }
};

class FreqThreadWorkData
{
public:
    struct MappedBaseModInfo {
        int sid;
        int soff;
        u8 scaled_prob;
        u8 motif;

        MappedBaseModInfo() {}

        MappedBaseModInfo(int _sid, int _soff, u8 _scaled_prob)
            : sid(_sid), soff(_soff), scaled_prob(_scaled_prob) {}

        bool operator<(const MappedBaseModInfo& rhs) const {
            return sid < rhs.sid;
        }
    };

    FreqThreadWorkData(samFile* sam_fp,
        sam_hdr_t* sam_hdr,
        HbnAggrMethyCallOptions* opts,
        HbnDatabase* db,
        size_t* chr_mod_counts,
        FILE* mod_out)
    {
        M_sam_fp = sam_fp;
        M_sam_hdr = sam_hdr;
        M_opts = opts;
        M_updb = db;
        M_chr_mod_counts = chr_mod_counts;
        for (int i = 0; i < 256; ++i) M_cpg_bins[i] = M_chg_bins[i] = M_chh_bins[i] = 0;
        M_mod_out = mod_out;
    }

    bool get_next_bam(bam1_t* bam) {
        std::lock_guard<std::mutex> __(M_mutex);
        int r = sam_read1(M_sam_fp, M_sam_hdr, bam);
        if (r == -1) return false;
        if (r < 0) HBN_ERR("Could not read BAM reacord");
        ++M_sam_cnt;
        return true;
    }

    void add_mods(std::vector<MappedBaseModInfo>& mods,
        size_t* a, size_t* b, size_t* c) {
        std::lock_guard<std::mutex> __(M_mutex);
        M_mods.insert(M_mods.end(), mods.begin(), mods.end());
        for (int i = 0; i < 256; ++i) {
            M_cpg_bins[i] += a[i];
            M_chg_bins[i] += b[i];
            M_chh_bins[i] += c[i];
        }
    }

    void dump_mods() {
        pdqsort(M_mods.begin(), M_mods.end());
        for (auto& M : M_mods) ++M_chr_mod_counts[M.sid];
        MappedBaseModInfo* A = M_mods.data();
        size_t C = M_mods.size();
        //HBN_LOG("Dump %zu mods for %zu reads", C, M_sam_cnt);
        hbn_fwrite(A, sizeof(MappedBaseModInfo), C, M_mod_out);
        //HBN_LOG("Done");
        M_mods.clear();
    }

public:
    samFile*    M_sam_fp;
    sam_hdr_t*  M_sam_hdr;
    size_t  M_sam_cnt;

    HbnAggrMethyCallOptions* M_opts;
    HbnDatabase*    M_updb;
    size_t*         M_chr_mod_counts;

    std::vector<MappedBaseModInfo>  M_mods;
    size_t          M_cpg_bins[256];
    size_t          M_chg_bins[256];
    size_t          M_chh_bins[256];
    std::mutex      M_mutex;
    FILE*           M_mod_out;
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
s_genomic_methy_freq_thread(void* params)
{
    FreqThreadWorkData* data = (FreqThreadWorkData*)(params);
    BamQuerySequence query;
    BamMapInfo align;
    MethylationContext methy_ctx;
    std::vector<FreqThreadWorkData::MappedBaseModInfo> chr_mods;
    std::vector<BaseModInfo> base_mods;
    std::vector<ReadBaseModInfo> read_mods;
    std::vector<MotifMappedHitInfo> mapped_samples;

    size_t cpg_bins[256];
    size_t chg_bins[256];
    size_t chh_bins[256];
    std::fill(cpg_bins, cpg_bins + 256, 0);
    std::fill(chg_bins, chg_bins + 256, 0);
    std::fill(chh_bins, chh_bins + 256, 0);

    bam1_t* bam = bam_init1();
    FreqThreadWorkData::MappedBaseModInfo M;
    int cnt = 0;
    while (data->get_next_bam(bam)) {
        base_mods.clear();
        extract_bam_base_mods(bam, base_mods);
        if (base_mods.empty()) continue;
        if (!query.init(bam)) continue;
        if (!align.init(data->M_sam_hdr, bam, data->M_updb, &query)) continue;

        if (!(bam->core.flag & 0x900)) {
            for (auto& m : base_mods) {
		    if (m.unmod_base != 'C' && m.unmod_base != 'G') continue;
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

        if (bam->core.qual < data->M_opts->min_mapQ) continue;
        if (align.pi < data->M_opts->min_pi) continue;

        read_mods.resize(query.size);
        for (auto& r : read_mods) r.init();
        for (auto& m : base_mods) {
            if (m.code != 'm') continue;
	        int qoff = m.qoff;
            read_mods[qoff].has_prob_value = true;
            read_mods[qoff].scaled_prob = m.scaled_prob;
        }

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
                M.scaled_prob = read_mods[qoff].scaled_prob;
                M.sid = align.sid;
                M.soff = soff;
                M.motif = 0;
                chr_mods.push_back(M);
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
                M.scaled_prob = read_mods[qoff].scaled_prob;
                M.sid = align.sid;
                M.soff = soff;
                M.motif = 1;
                chr_mods.push_back(M);
            }  
        }

        mapped_samples.clear();
        extract_chh_mapped_samples(data->M_updb, methy_ctx, query, align, mapped_samples);
        for (auto& m : mapped_samples) {
            if (read_mods[m.qoff].has_prob_value) {
                M.scaled_prob = read_mods[m.qoff].scaled_prob;
                M.sid = m.sid;
                M.soff = m.soff;
                M.motif = 2;
                chr_mods.push_back(M);
            }
        }
        if (++cnt == 100) break;
    }
    data->add_mods(chr_mods, cpg_bins, chg_bins, chh_bins);
    bam_destroy1(bam);
    return nullptr;
}

static void
s_resolve_scaled_prob_threshold(size_t* cpg_bins, size_t* chg_bins, size_t* chh_bins,
    u8& cpg_threshold, u8& chg_threshold, u8& chh_threshold)
{
    const size_t* a = cpg_bins;
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

    a = chg_bins;
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

    a = chh_bins;
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

static bool
s_bam_is_mapped_and_sorted(sam_hdr_t* hdr)
{
    bool bam_is_sorted = false;
    bool bam_is_mapped = false;

    kstring_t str = KS_INITIALIZE;
    if (sam_hdr_find_tag_hd(hdr, "SO", &str) == 0) {
        if (strcmp(str.s, "coordinate") == 0) bam_is_sorted = true;
    }
    ks_free(&str);

    if (hdr->n_targets > 0) bam_is_mapped = true;

    if ((!bam_is_mapped) || (!bam_is_sorted)) {
        fprintf(stderr, "ERROR: Methylation frequency could not be computed due to the following errors:\n");
        if (!bam_is_mapped) fprintf(stderr, "BAM is not mapped\n");
        if (!bam_is_sorted) fprintf(stderr, "BAM is not sorted\n");
        return false;
    }
    return true;
}

static int
s_compute_methy_freq(int argc, char* argv[])
{
    HbnAggrMethyCallOptions options;
    if (!options.parse(argc, argv)) {
        options.dump_usage(argc, argv);
        return EXIT_FAILURE;
    }
    options.dump_parameters();
    HbnProgramInfo hpi(HBN_PACKAGE_NAME, nullptr);

    samFile* sam_fp = sam_open(options.mod_bam_path, "rb");
    hts_set_threads(sam_fp, 8);
    sam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    if (!s_bam_is_mapped_and_sorted(sam_hdr)) return 1;

    char path[HBN_MAX_PATH_LEN];
    snprintf(path, HBN_MAX_PATH_LEN, "%s.CpG.cov.bed", options.output_prefix);
    hbn_dfopen(cpgout, path, "w");
    snprintf(path, HBN_MAX_PATH_LEN, "%s.CHG.cov.bed", options.output_prefix);
    hbn_dfopen(chgout, path, "w");
    snprintf(path, HBN_MAX_PATH_LEN, "%s.CHH.cov.bed", options.output_prefix);
    hbn_dfopen(chhout, path, "w");
    
    std::string base_mod_path("read_base_mods");
    hbn_dfopen(read_mod_out, base_mod_path.c_str(), "wb");

    HbnDatabase* updb = new HbnDatabase(options.reference_path);
    int num_chr = updb->num_seqs();
    size_t* chr_mod_counts = (size_t*)calloc(num_chr, sizeof(size_t));
    FreqThreadWorkData data(sam_fp, sam_hdr, &options, updb, chr_mod_counts, read_mod_out);
    pthread_t jobs[options.num_threads];

    while (1) {
        data.M_sam_cnt = 0;
        for (int i = 0; i < options.num_threads; ++i) {
            pthread_create(jobs + i, nullptr, s_genomic_methy_freq_thread, &data);
        }
        for (int i = 0; i < options.num_threads; ++i) {
            pthread_join(jobs[i], nullptr);
        }
        data.dump_mods();
        if (data.M_sam_cnt == 0) break;
    }
    hbn_fclose(read_mod_out);
    sam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);

    u8 cpg_threshold = 0, chg_threshold = 0, chh_threshold = 0;
    s_resolve_scaled_prob_threshold(data.M_cpg_bins, data.M_chg_bins, data.M_chh_bins,
        cpg_threshold, chg_threshold, chh_threshold);

    hbn_dfopen(read_mod_in, base_mod_path.c_str(), "rb");
    for (int i = 0; i < num_chr; ++i) {
        if (!chr_mod_counts[i]) continue;
        const char* chr_name = updb->seq_name(i);
        const int chr_size = updb->seq_length(i);
        //HBN_LOG("Dump frequency for %s:%d", chr_name, chr_size);
        u8* base_motifs = (u8*)calloc(chr_size, sizeof(u8));
        int* pcov = (int*)calloc(chr_size, sizeof(int));
        int* ncov = (int*)calloc(chr_size, sizeof(int));

        size_t left = chr_mod_counts[i];
        const size_t kBufSize = 1000000;
        FreqThreadWorkData::MappedBaseModInfo* buf = new FreqThreadWorkData::MappedBaseModInfo[kBufSize];
        while (left) {
            size_t N = std::min(left, kBufSize);
            int n = fread(buf, sizeof(FreqThreadWorkData::MappedBaseModInfo), N, read_mod_in);
            for (int k = 0; k < n; ++k) {
                hbn_assert(buf[k].sid == i);
                u8 prob = buf[k].scaled_prob;
                int soff = buf[k].soff;
                if (buf[k].motif == 0) {
                    base_motifs[soff] = 0;
                    if (prob >= cpg_threshold) {
                        ++pcov[soff];
                    } else {
                        ++ncov[soff];
                    }
                } else if (buf[k].motif == 1) {
                    base_motifs[soff] = 1;
                    if (prob >= chg_threshold) {
                        ++pcov[soff];
                    } else {
                        ++ncov[soff];
                    }
                } else if (buf[k].motif == 2) {
                    base_motifs[soff] = 2;
                    if (prob >= chh_threshold) {
                        ++pcov[soff];
                    } else {
                        ++ncov[soff];
                    }
                } else {
                    HBN_ERR("Invalid motif idx %d", buf[k].motif);
                }
            }
            left -= n;
        }
        delete[] buf;

        std::ostringstream cpg, chg, chh;
        for (int k = 0; k < chr_size; ++k) {
            int cov = pcov[k] + ncov[k];
            if (!cov) continue;
            double freq = 100.0 * pcov[k] / cov;
            if (base_motifs[k] == 0) {
                cpg << chr_name
                    << '\t' << k << '\t' << k+1
                    << '\t' << freq
                    << '\t' << pcov[k] << '\t' << ncov[k]
                    << '\n';
            } else if (base_motifs[k] == 1) {
                chg << chr_name
                    << '\t' << k << '\t' << k+1
                    << '\t' << freq
                    << '\t' << pcov[k] << '\t' << ncov[k]
                    << '\n';                
            } else if (base_motifs[k] == 2) {
                chh << chr_name
                    << '\t' << k << '\t' << k+1
                    << '\t' << freq
                    << '\t' << pcov[k] << '\t' << ncov[k]
                    << '\n';
            }
        }
        
        if (!cpg.str().empty()) hbn_fwrite(cpg.str().c_str(), 1, cpg.str().size(), cpgout);
        if (!chg.str().empty()) hbn_fwrite(chg.str().c_str(), 1, chg.str().size(), chgout);
        if (!chh.str().empty()) hbn_fwrite(chh.str().c_str(), 1, chh.str().size(), chhout);

        free(base_motifs);
        free(pcov);
        free(ncov);
    }

    delete updb;
    hbn_fclose(read_mod_in);
    hbn_fclose(cpgout);
    hbn_fclose(chgout);
    hbn_fclose(chhout);

    std::remove(base_mod_path.c_str());

    return 0;
}

} // ns_pileup

int pileup_main(int argc, char* argv[])
{
    return ns_pileup::s_compute_methy_freq(argc, argv);
}
