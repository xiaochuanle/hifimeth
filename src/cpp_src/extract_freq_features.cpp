#include "5mc_aux.hpp"
#include "5mc_tags.hpp"
#include "line_reader.h"
#include "freq_features_aux.hpp"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

using namespace std;

struct ChrModInfo 
{
    int sid;
    int soff;
    double prob;
};

static void
s_load_mod_prob(const char* mod_path, map<string, int>& chr_name2id, vector<string>& chr_name_list, vector<ChrModInfo>& mod_list)
{
    HbnLineReader* in = HbnLineReaderNew(mod_path);
    vector<pair<const char*, int>> cols;
    string chr_name, last_chr_name;
    int last_chr_id;
    ChrModInfo mod;
    int cnt = 0;
    while (!HbnLineReaderAtEof(in)) {
        HbnLineReaderReadOneLine(in);
        const char* line = ks_s(in->line);
        const int line_size = ks_size(in->line);
        cols.clear();
        split_string_by_char(line, line_size, '\t', cols);
        if (cols.size() == 5) continue;
        chr_name.assign(cols[4].first, cols[4].second);
        if (chr_name != last_chr_name) {
            last_chr_name = chr_name;
            auto pos = chr_name2id.find(last_chr_name);
            if (pos != chr_name2id.end()) {
                last_chr_id = pos->second;
            } else {
                last_chr_id = chr_name_list.size();
                chr_name_list.push_back(last_chr_name);
                chr_name2id.insert(pair<string, int>(last_chr_name, last_chr_id));
            }
        }
        mod.sid = last_chr_id;
        mod.soff = atoi(cols[5].first);
        mod.prob = atof(cols[7].first);
        mod_list.push_back(mod);
        if (++cnt < 10) fprintf(stderr, "%s\t%d\t%d\t%g\n", last_chr_name.c_str(), mod.sid, mod.soff, mod.prob);
    }
    HbnLineReaderFree(in);

    sort(mod_list.begin(), mod_list.end(), [](const ChrModInfo& x, const ChrModInfo& y) {
        return (x.sid < y.sid) || (x.sid == y.sid && x.soff < y.soff);
    });
    HBN_LOG("Load %zu 5mc call results\n", mod_list.size());
}

static void
s_load_label_from_bed(const char* bed_path, 
    const char* target_chr_name,
    const int min_cov, 
    vector<pair<int, double>>& label_list)
{
    label_list.clear();

    HbnLineReader* in = HbnLineReaderNew(bed_path);
    vector<pair<const char*, int>> cols;
    string chr_name;
    int cnt = 0;
    while (!HbnLineReaderAtEof(in)) {
        HbnLineReaderReadOneLine(in);
        const char* line = ks_s(in->line);
        const int line_size = ks_size(in->line);
        cols.clear();
        split_string_by_char(line, line_size, '\t', cols);
        chr_name.assign(cols[0].first, cols[0].second);
        if (chr_name != target_chr_name) continue;
        int soff = atoi(cols[1].first);
        int c1 = atoi(cols[4].first);
        int c2 = atoi(cols[5].first);
        if (c1 + c2 < min_cov) continue;
        double p = atof(cols[3].first);
	    if (c1 + c2 < 50 || c1 + c2 > 150) continue;
        label_list.push_back(pair<int, double>(soff, p));
	    if (++cnt < 10) {
            hbn_fwrite(line, 1, line_size, stderr);
            fprintf(stderr, "\n");
		    fprintf(stderr, "soff = %d, c1 = %d, c2 = %d, p = %g\n", soff, c1, c2, p);
	    }
    }
    HBN_LOG("Load %zu sites for %s\n", label_list.size(), target_chr_name);
}

struct SiteModStat { int soff; int pos_cov; int neg_cov; };

static void
s_extract_train_features_for_one_chr(const char* bed_path, 
    const char* target_chr_name,
    ChrModInfo* mia, int mic,
    const int min_cov, const int bins, const int num_sites,
    FILE* features_out)
{
    ProbHist prob_hist(bins);
    vector<double> prob_list;
    vector<double> prob_features;
    map<int, int> pos_ids;
    vector<double> feature_list;
    vector<SiteModStat> modstat_list;
    SiteModStat modstat;

    int pos_id = 0;
    int i = 0;
    while (i < mic) {
        int j = i + 1;
        while (j < mic && mia[i].soff == mia[j].soff) ++j;
        if (j - i >= min_cov) {
            prob_list.clear();
            for (size_t k = i; k < j; ++k) prob_list.push_back(mia[k].prob);
            prob_features.clear();
            prob_hist.process_prob_list(prob_list.data(), prob_list.size(), prob_features);    
            feature_list.insert(feature_list.end(), prob_features.begin(), prob_features.end());  
            int c1 = 0, c2 = 0;
            for (int k = i; k < j; ++k) if (mia[k].prob < 0.5) ++c1; else ++c2;
            modstat.soff = mia[i].soff;
            modstat.pos_cov = c2;
            modstat.neg_cov = c1;
            modstat_list.push_back(modstat);
            pos_ids[mia[i].soff] = pos_id++;      
        }
        i = j;
    }

    vector<pair<int, double>> label_list;
    s_load_label_from_bed(bed_path, target_chr_name, min_cov, label_list);
    int flanking_sites = num_sites / 2;
    int cnt = 0;
    vector<double> sample;
    for (auto& example : label_list) {
        int soff = example.first;
        auto iter = pos_ids.find(soff);
        if (iter == pos_ids.end()) continue;
        int pos_id = iter->second;
        hbn_assert(soff == modstat_list[pos_id].soff);
	//int cov = modstat_list[pos_id].pos_cov + modstat_list[pos_id].neg_cov;
	//if (cov < min_cov) continue;
        int from = pos_id - flanking_sites, to = pos_id + flanking_sites + 1;
        sample.clear();
        for (int i = from; i < to; ++i) {
            int ii = i;
            if (ii < 0) ii = 0;
            if (ii >= mic) ii = mic-1;
            size_t offset = ii * bins;
            double dist = fabs(soff - modstat_list[ii].soff);
            sample.push_back(dist);
            for (int k =0; k < bins; ++k) sample.push_back(feature_list[offset+k]);
        }
        for (auto f : sample) fprintf(features_out, "%lf\t", f);
        fprintf(features_out, "%lf\n", example.second);
        ++cnt;
    }
    fprintf(stderr, "Build %d samples\n", cnt);
}

static void
s_extract_call_features_for_one_chr(const char* target_chr_name,
    ChrModInfo* mia, int mic,
    const int min_cov, const int bins, const int num_sites,
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
        for (int k = i; k < j; ++k) if (mia[k].prob < 0.5) ++c1; else ++c2;
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
s_dump_freq_call_for_one_chr(const char* target_chr_name,
    ChrModInfo* mia, int mic, double meth_prob, double unmeth_prob,
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

static void
s_dump_usage(int argc, char* argv[])
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s [OPTIONS] mod-path output [bed-path]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS:\n");
    fprintf(stderr, "  -c Minimum coverage for building model features\n");
    fprintf(stderr, "    Default = '%d'\n", kMinCov);
    fprintf(stderr, "  -s Number of motif sites in one sample\n");
    fprintf(stderr, "    Default = '%d'\n", kNumSites);
    fprintf(stderr, "  -m Methylation probablity cutoff\n");
    fprintf(stderr, "    Default = '%g'\n", kMethProb);
    fprintf(stderr, "  -u Unmethlylation probability cutoff\n");
    fprintf(stderr, "    Default = '%g'\n", kUnmethProb);
    fprintf(stderr, "  -r Methylation frequency path\n");
    fprintf(stderr, "\n");
}

static bool
s_parse_arguments(int argc, char* argv[],
    int* min_cov,
    int* num_sites,
    double* meth_prob,
    double* unmeth_prob,
    const char** freq_path,
    const char** mod_path,
    const char** output,
    const char** bed_path)
{
    *min_cov = kMinCov;
    *num_sites = kNumSites;
    *meth_prob = kMethProb;
    *unmeth_prob = kUnmethProb;
    *freq_path = NULL;
    *mod_path = NULL;
    *output = NULL;
    *bed_path = NULL;

    int i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') break;
        if (argv[i][0] == '-' && strlen(argv[i]) == 1) break;

        if (strcmp(argv[i], "-c") == 0) {
            *min_cov = atoi(argv[i+1]);
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

        if (strcmp(argv[i], "-r") == 0) {
            *freq_path = argv[i+1];
            i += 2;
            continue;
        }     

        fprintf(stderr, "unrecognised option '%s'\n", argv[i]);
        return false;
    }
    if (i + 2 != argc && i + 3 != argc) return false;
    *mod_path = argv[i];
    *output = argv[i+1];
    if (i +3 == argc) *bed_path = argv[i+2];
    return true;
}

int main(int argc, char* argv[])
{
    int bins = 20;
    int min_cov = 0;
    int num_sites = 0;
    double meth_prob = 0;
    double unmeth_prob = 0;
    const char* freq_path = NULL;
    const char* mod_path = NULL;
    const char* output = NULL;
    const char* bed_path = NULL;
    if (!s_parse_arguments(argc, argv, &min_cov, &num_sites, &meth_prob, &unmeth_prob, &freq_path, &mod_path, &output, &bed_path)) {
        s_dump_usage(argc, argv);
        return 1;
    }

    fprintf(stderr, "Min cov: %d\n", min_cov);
    fprintf(stderr, "num sites: %d\n", num_sites);
    fprintf(stderr, "5mc-calls: %s\n", mod_path);
    fprintf(stderr, "output: %s\n", output);
    if (bed_path) fprintf(stderr, "bed path: %s\n", bed_path);
    if (freq_path) {
        fprintf(stderr, "freq-calls: %s\n", freq_path);
        fprintf(stderr, "meth prob: %g\n", meth_prob);
        fprintf(stderr, "unmeth prob: %g\n", unmeth_prob);
    }

    map<string, int> chr_name2id;
    vector<string> chr_name_list;
    vector<ChrModInfo> mod_list;
    s_load_mod_prob(mod_path, chr_name2id, chr_name_list, mod_list);
    hbn_dfopen(features_out, output, "w");
    FILE* freq_out = NULL;
    if (freq_path) hbn_fopen(freq_out, freq_path, "w");

    ChrModInfo* mia = mod_list.data();
    size_t mic = mod_list.size();
    size_t i = 0;
    while (i < mic) {
        size_t j = i + 1;
        while (j < mic && mia[i].sid == mia[j].sid) ++j;

        const char* chr_name = chr_name_list[mia[i].sid].c_str();
        HBN_LOG("extract samples for %s\n", chr_name);
        if (freq_out) s_dump_freq_call_for_one_chr(chr_name, mia + i, j - i, meth_prob, unmeth_prob, freq_out);
        if (bed_path) {
            s_extract_train_features_for_one_chr(bed_path, chr_name, mia + i, j -i, min_cov, bins, num_sites, features_out);
        } else {
            s_extract_call_features_for_one_chr(chr_name, mia + i, j -i, min_cov, bins, num_sites, features_out);
        }

        i = j;
    }
    hbn_fclose(features_out);
    if (freq_out) hbn_fclose(freq_out);

    return 0;
}
