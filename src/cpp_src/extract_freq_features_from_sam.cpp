#include "5mc_aux.hpp"
#include "5mc_tags.hpp"
#include "line_reader.h"
#include "freq_features_aux.hpp"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <map>
#include <mutex>
#include <string>
#include <vector>

using namespace std;

#define HP0 0
#define HP1 1
#define HP2 2

static const int kMinMapQ = 10;

struct ModCallInfo
{
    int sid;
    int soff;
    double prob;
    int hp;
};

struct ModCallLoadThreadData
{
    const char* motif;
    int motif_size;
    int min_mapq;

    Reference* ref;
    HbnLineReader* in;
    mutex in_lock;
    vector<ModCallInfo>* call_list;
    mutex call_list_lock;

    bool get_next_sam(string& sam) {
        lock_guard<mutex> lg(in_lock);
        sam.clear();
        bool sam_is_loaded = false;
        while (!HbnLineReaderAtEof(in)) {
            HbnLineReaderReadOneLine(in);
            if (ks_empty(in->line)) continue;
            if (ks_front(in->line) == '@') continue;
            sam_is_loaded = true;
            sam.assign(ks_s(in->line), ks_size(in->line));
            break;
        }
        return sam_is_loaded;
    }

    void add_calls(vector<ModCallInfo>& read_call_list) {
        lock_guard<mutex> lg(call_list_lock);
        call_list->insert(call_list->end(), read_call_list.begin(), read_call_list.end());
    }
};

static void
s_extract_5mc_calls(const char* sam, const int sam_size, vector<pair<const char*, int>>& sam_cols, 
    const char* fwd_read, const int read_size,
    const char* motif, const int motif_size, 
    const int ref_id, map<int, int>& mapinfo,
    vector<ModCallInfo>& call_list)
{
    const char* mmtag = "MM:Z:C+m,";
    const int mmtag_size = strlen(mmtag);
    const char* mltag = "ML:B:C,";
    const int mltag_size = strlen(mltag);

    int hp = HP0;
    if (strstr(sam, "HP:i:1")) {
        hp = HP1;
    } else if (strstr(sam, "HP:i:2")) {
        hp = HP2;
    }

    call_list.clear();
    const char* mmstr = NULL;
    int mmstr_size = 0;
    const char* mlstr = NULL;
    int mlstr_size = 0;
    const int n_col = sam_cols.size();
    for (int i = 11; i < n_col && (mmstr == NULL || mlstr == NULL); ++i) {
        if (mmstr == NULL && sam_cols[i].second > mmtag_size && strncmp(mmtag, sam_cols[i].first, mmtag_size) == 0) {
            mmstr = sam_cols[i].first + mmtag_size;
            mmstr_size = sam_cols[i].second - mmtag_size;
            continue;
        }
        if (mlstr == NULL && sam_cols[i].second > mltag_size && strncmp(mltag, sam_cols[i].first, mltag_size) == 0) {
            mlstr = sam_cols[i].first + mltag_size;
            mlstr_size = sam_cols[i].second - mltag_size;
            continue;
        }
    }
    if (mmstr == NULL || mlstr == NULL) {
        //fprintf(stderr, "5mc tag is missing\n");
        return;
    }

    vector<pair<const char*, int>> mm_cols, ml_cols;
    split_string_by_char(mmstr, mmstr_size, ',', mm_cols);
    split_string_by_char(mlstr, mlstr_size, ',', ml_cols);
    if (mm_cols.size() != ml_cols.size()) {
        hbn_fwrite(sam, 1, sam_size, stderr); fprintf(stderr, "\n");
        HBN_LOG("motif counts in %s and in %s do not match (%zu v.s. %zu)", mmtag, mltag, mm_cols.size(), ml_cols.size());
        abort();
    }
    if (mm_cols.empty()) return;

    const int n_motif = mm_cols.size();
    int read_i = 0;
    for (int i = 0; i < n_motif; ++i) {
        const int total_skipped_c = atoi(mm_cols[i].first);
        int skipped_c = 0;
        while (skipped_c < total_skipped_c) {
            hbn_assert(read_i < read_size);
            if (fwd_read[read_i] == motif[0]) ++skipped_c;
            ++read_i;
        }
        hbn_assert(read_i < read_size);
        while (fwd_read[read_i] != motif[0]) ++read_i;
        hbn_assert(read_i < read_size);
        hbn_assert(fwd_read[read_i] == motif[0]);
        for (int k = 0; k < motif_size; ++k) {
            hbn_assert(fwd_read[read_i+k] == motif[k]);
        }
        auto pos = mapinfo.find(read_i);
        read_i += motif_size;
        if (pos == mapinfo.end()) continue;

        double p = atoi(ml_cols[i].first);
        p = p / 255.0;
        p = max(0.0, p);
        p = min(1.0, p);
        ModCallInfo mci;
        mci.sid = ref_id;
        mci.soff = pos->second;
        mci.prob = p;
        mci.hp = hp;
        call_list.push_back(mci);
    }
}

static void*
s_5mc_extract_thread(void* params)
{
    ModCallLoadThreadData* data = (ModCallLoadThreadData*)(params);
    Reference* ref = data->ref;
    const char* motif = data->motif;
    const int motif_size = data->motif_size;
    const int min_mapq = data->min_mapq;

    string line;
    vector<pair<const char*, int>> cols;
    vector<ModCallInfo> call_list;
    string fwd_read;
    string ref_name;
    map<int, int> mapinfo;

    while (data->get_next_sam(line)) {
        cols.clear();
        split_string_by_char(line.c_str(), line.size(), '\t', cols);

        const int read_size = cols[9].second;
        int sam_flag = atoi(cols[1].first);
        int strand = (sam_flag&16) ? REV : FWD;
        extract_forward_read_from_sam(cols[9].first, read_size, strand, fwd_read);
        const char* read = fwd_read.c_str();

        mapinfo.clear();
        int qdir = FWD;
        const char* sam = line.c_str();
        const int sam_size = line.size();
        extract_mapinfo(sam, sam_size, ref, motif, motif_size, min_mapq, ref_name, qdir, mapinfo); 
        if (mapinfo.empty()) {
            //fprintf(stderr, "empty map info\n");
            continue;
        }

        const int ref_id = ref->ref_id(ref_name);
        call_list.clear();
        s_extract_5mc_calls(sam, sam_size, cols, read, read_size, motif, motif_size, ref_id, mapinfo, call_list);
        if (call_list.empty()) {
            //fprintf(stderr, "no call is found\n");
            continue;
        }

        data->add_calls(call_list);
    }

    return NULL;
}

static void
s_dump_usage(int argc, char* argv[])
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s [OPTIONS] sam-path reference output\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS:\n");
    fprintf(stderr, "  -q <integer>\n");
    fprintf(stderr, "    Minimum mapping quality\n");
    fprintf(stderr, "    Default = '%d'\n", kMinMapQ);
    fprintf(stderr, "  -M <motif>\n");
    fprintf(stderr, "    The motif sequence\n");
    fprintf(stderr, "    Default = '%s'\n", default_motif());
    fprintf(stderr, "  -c <integer> Minimum coverage for building model features\n");
    fprintf(stderr, "    Default = '%d'\n", kMinCov);
    fprintf(stderr, "  -s <integer> Number of motif sites in one sample\n");
    fprintf(stderr, "    Default = '%d'\n", kNumSites);
    fprintf(stderr, "  -m <real> Methylation probablity cutoff\n");
    fprintf(stderr, "    Default = '%g'\n", kMethProb);
    fprintf(stderr, "  -u <real> Unmethlylation probability cutoff\n");
    fprintf(stderr, "    Default = '%g'\n", kUnmethProb);
    fprintf(stderr, "  -p Dump paternal and maternal freq info\n");
    fprintf(stderr, "  -t <integer>\n");
    fprintf(stderr, "    Number of CPU threads\n");
    fprintf(stderr, "\n");
}

static bool
s_parse_arguments(int argc, char* argv[],
    int* min_mapq,
    string& motif,
    int* min_cov,
    int* num_sites,
    double* meth_prob,
    double* unmeth_prob,
    bool* dump_phased_freq,
    int* num_threads,
    const char** sam_path,
    const char** ref_path,
    const char** output)
{
    *min_mapq = kMinMapQ;
    motif = default_motif();
    *min_cov = kMinCov;
    *num_sites = kNumSites;
    *meth_prob = kMethProb;
    *unmeth_prob = kUnmethProb;
    *dump_phased_freq = false;
    *num_threads = 1;
    *sam_path = NULL;
    *ref_path = NULL;
    *output = NULL;

    int i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') break;
        if (argv[i][0] == '-' && strlen(argv[i]) == 1) break;

        if (strcmp(argv[i], "-c") == 0) {
            *min_cov = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-q") == 0) {
            *min_mapq = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-M") == 0) {
            motif = argv[i+1];
            for (auto& c : motif) c = toupper(c);
            i += 2;
            continue;
        }   

        if (strcmp(argv[i], "-s") == 0) {
            *num_sites = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-m") == 0) {
            *meth_prob = atof(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-u") == 0) {
            *unmeth_prob = atof(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-p") == 0) {
            *dump_phased_freq = true;
            i += 1;
            continue;
        }    

        if (strcmp(argv[i], "-t") == 0) {
            *num_threads = atoi(argv[i+1]);
            i += 2;
            continue;
        } 

        fprintf(stderr, "unrecognised option '%s'\n", argv[i]);
        return false;
    }
    if (i + 3 != argc) return false;
    *sam_path = argv[i];
    *ref_path = argv[i+1];
    *output = argv[i+2];
    return true;
}

static void
s_dump_freq_call_for_one_chr(const char* target_chr_name,
    ModCallInfo* mia, int mic, double meth_prob, double unmeth_prob,
    FILE* freq_out)
{
    int i = 0;
    while (i < mic) {
        int j = i + 1;
        while (j < mic && mia[i].soff == mia[j].soff) ++j;
        int c1 = 0, c2 = 0;
        for (int k = i; k < j; ++k) {
            if (mia[k].prob >= meth_prob) {
                ++c1;
            } else if (mia[k].prob < unmeth_prob) {
                ++c2;
            }
        }
        double p = 1.0 * c1 / (c1 + c2);
        fprintf(freq_out, "%s\t%d\t%d\t%.4f\t%d\t%d\n", target_chr_name, mia[i].soff, mia[i].soff + 1, p, c1, c2);

        i = j;
    }
}

struct SiteModStat { int soff; int pos_cov; int neg_cov; };

static void
s_extract_call_features_for_one_chr(const char* target_chr_name,
    ModCallInfo* mia, int mic,
    const int bins, const int num_sites,
    double meth_prob, double unmeth_prob,
    FILE* features_out)
{
    ProbHist prob_hist(bins);
    vector<double> prob_list;
    vector<double> prob_features;
    vector<double> feature_list;
    vector<SiteModStat> modstat_list;
    SiteModStat modstat;

    int i = 0;
    int sample_cnt = 0;
    while (i < mic) {
        int j = i + 1;
        while (j < mic && mia[i].soff == mia[j].soff) ++j;
        prob_list.clear();
        for (size_t k = i; k < j; ++k) prob_list.push_back(mia[k].prob);
        prob_features.clear();
        prob_hist.process_prob_list(prob_list.data(), prob_list.size(), prob_features);    
        feature_list.insert(feature_list.end(), prob_features.begin(), prob_features.end());  
        int c1 = 0, c2 = 0;
        for (int k = i; k < j; ++k) {
            if (mia[k].prob >= meth_prob) ++c1; else if (mia[k].prob < unmeth_prob) ++c2;
        }
        modstat.soff = mia[i].soff;
        modstat.pos_cov = c2;
        modstat.neg_cov = c1;
        modstat_list.push_back(modstat);
        ++sample_cnt;
        i = j;
    }

    int flanking_sites = num_sites / 2;
    for (int p = 0; p < sample_cnt; ++p) {
        int from = p - flanking_sites, to = p + flanking_sites + 1;
        int soff = modstat_list[p].soff;
        fprintf(features_out, "%s\t%d\t%d\t%d", target_chr_name, soff, modstat_list[p].pos_cov, modstat_list[p].neg_cov);
        for (int i = from; i < to; ++i) {
            int ii = i;
            if (ii < 0) ii = 0;
            if (ii >= mic) ii = mic-1;
            size_t offset = ii * bins;
            double dist = fabs(soff - modstat_list[ii].soff);
            fprintf(features_out, "\t%lf", dist);
            for (int k = 0; k < bins; ++k) fprintf(features_out, "\t%lf", feature_list[offset+k]);
        }     
        fprintf(features_out, "\n");  
    }
    fprintf(stderr, "Build %d samples\n", sample_cnt);
}

static void
s_process_call_list(Reference* ref,
    ModCallInfo* mia,
    size_t mic,
    int bins, int num_sites,
    double meth_prob, double unmeth_prob,
    const char* output_dir,
    const char* output_prefix)
{
    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/%s.count.txt", output_dir, output_prefix);
    hbn_dfopen(freq_count_out, path, "w");
    sprintf(path, "%s/%s.model.features", output_dir, output_prefix);
    hbn_dfopen(freq_model_features_out, path, "w");

    sort(mia, mia + mic, [](const ModCallInfo& x, const ModCallInfo& y) { return (x.sid < y.sid) || (x.sid == y.sid && x.soff < y.soff); });
    size_t i = 0;
    while (i < mic) {
        size_t j = i + 1;
        while (j < mic && mia[i].sid == mia[j].sid) ++j;

        const char* chr_name = ref->ref_name(mia[i].sid);
        HBN_LOG("extract samples for %s", chr_name);
        s_dump_freq_call_for_one_chr(chr_name, mia + i, j - i, meth_prob, unmeth_prob, freq_count_out);
        s_extract_call_features_for_one_chr(chr_name, mia + i, j - i, bins, num_sites, meth_prob, unmeth_prob, freq_model_features_out);

        i = j;
    }

    hbn_fclose(freq_model_features_out);
    hbn_fclose(freq_count_out);
}

int main(int argc, char* argv[])
{
    int bins = 20;
    int min_mapq = 0;
    string motif;
    int min_cov = 0;
    int num_sites = 0;
    double meth_prob = 0;
    double unmeth_prob = 0;
    bool dump_phased_freq = false;
    int num_threads = 1;
    const char* sam_path = NULL;
    const char* ref_path = NULL;
    const char* output = NULL;
    if (!s_parse_arguments(argc, argv, &min_mapq, motif, &min_cov, &num_sites, &meth_prob, &unmeth_prob, &dump_phased_freq, &num_threads, &sam_path, &ref_path, &output)) {
        s_dump_usage(argc, argv);
        return 1;
    }

    fprintf(stderr, "Minimum mapping quality: %d\n", min_mapq);
    fprintf(stderr, "Motif: %s\n", motif.c_str());
    fprintf(stderr, "CPU threads: %d\n", num_threads);
    fprintf(stderr, "num sites: %d\n", num_sites);
    fprintf(stderr, "methylated cutoff: %g\n", meth_prob);
    fprintf(stderr, "unmethylated cutoff: %g\n", unmeth_prob);
    fprintf(stderr, "dump phased frequencies: %d\n", dump_phased_freq);
    fprintf(stderr, "SAM: %s\n", sam_path);
    fprintf(stderr, "reference: %s\n", ref_path);
    fprintf(stderr, "output: %s\n", output);

    Reference* ref = new Reference(ref_path);
    vector<ModCallInfo> call_list;
    HbnLineReader* in = HbnLineReaderNew(sam_path);

    ModCallLoadThreadData data;
    data.ref = ref;
    data.motif = motif.c_str();
    data.motif_size = motif.size();
    data.min_mapq = min_mapq;
    data.call_list = &call_list;
    data.in = in;
    
    pthread_t jobs[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(&jobs[i], NULL, s_5mc_extract_thread, &data);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobs[i], NULL);
    }
    in = HbnLineReaderFree(in);
    int hp1_calls = 0;
    int hp2_calls = 0;
    for (auto& x : call_list) {
        if (x.hp == HP1) ++hp1_calls;
        if (x.hp == HP2) ++hp2_calls;
    }
    HBN_LOG("Load %zu 5mc calls. HP1: %d, HP2: %d", call_list.size(), hp1_calls, hp2_calls);

    s_process_call_list(ref, call_list.data(), call_list.size(), bins, num_sites, meth_prob, unmeth_prob, output, "freq-call");

    if (dump_phased_freq) {
        vector<ModCallInfo> hp1_call_list, hp2_call_list;
        for (auto& x : call_list) if (x.hp == HP1) hp1_call_list.push_back(x); else if (x.hp == HP2) hp2_call_list.push_back(x);
        s_process_call_list(ref, hp1_call_list.data(), hp1_call_list.size(), bins, num_sites, meth_prob, unmeth_prob, output, "freq-call.hp1");
        s_process_call_list(ref, hp2_call_list.data(), hp2_call_list.size(), bins, num_sites, meth_prob, unmeth_prob, output, "freq-call.hp2");
    }

    delete ref;

    return 0;
}
