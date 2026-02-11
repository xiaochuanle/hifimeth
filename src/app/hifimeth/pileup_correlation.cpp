#include "../../corelib/arg_parse.hpp"
#include "../../corelib/hbn_aux.hpp"
#include "../../corelib/line_reader.hpp"
#include "../../corelib/pdqsort.h"
#include "program_info.hpp"

#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace ns_pileup_correlation {

static constexpr int kMinCov = 5;

struct HbnAggrMethyCorrOptions 
{
    int min_cov { kMinCov };

    const char* bed1 { nullptr };
    const char* bed2 { nullptr };

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

            if (parse_int_arg_value(argc, argv, i, "-c", min_cov)) continue;

            fprintf(stderr, "ERROR: unrecognised option %s", argv[i]);
            return false;
        }

        if (i >= argc) return false;
        bed1 = argv[i];
        ++i;

        if (i >= argc) return false;
        bed2 = argv[i];
        ++i;

        if (i != argc) return false;
        return true;
    }

    void dump_usage(int argc, char* argv[]) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "  %s %s [OPTIONS] bed1 bed2\n", argv[0], argv[1]);

        fprintf(stderr, "\n\n");
        fprintf(stderr, "DESCRIPTION:\n");
        fprintf(stderr, "  Compute pearson correlation between two aggregate cytosine methylation bed files\n");

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
        fprintf(stderr, "    Only consider genomic loci covered by at leat this number of reads\n");
        fprintf(stderr, "    Default: %d\n", kMinCov);
    }

    void dump_parameters() {
        fprintf(stderr, "\n\n");
        fprintf(stderr, "====================> Parameters:\n");
        fprintf(stderr, "min-cov: %d\n", min_cov);
        fprintf(stderr, "bed1: %s\n", bed1);
        fprintf(stderr, "bed2: %s\n", bed2);
        fprintf(stderr, "\n\n");
    }
};

static void
s_load_bed_methy_call(const char* path, const int min_cov, 
    std::map<std::string, int>& chr_name2id, std::vector<std::pair<uint64_t, double>>& methy)
{
    HbnLineReader in(path);
    std::string line;
    std::string last_chr;
    int last_sid = -1;
    while (!in.AtEOF()) {
        ++in;
        auto vline = *in;
        line.assign(vline.first, vline.second);
        auto cols = split_cstr(line.c_str(), line.size(), '\t');

        int pcov = atoi(cols[4].first);
        int ncov = atoi(cols[5].first);
        if (pcov + ncov < min_cov) continue;
        
        if ((cols[0].second != last_chr.size()) || strncmp(cols[0].first, last_chr.c_str(), cols[0].second)) {
            last_chr.assign(cols[0].first, cols[0].second);
            auto iter = chr_name2id.find(last_chr);
            if (iter == chr_name2id.end()) {
                last_sid = chr_name2id.size();
                chr_name2id[last_chr] = last_sid;
            } else {
                last_sid = iter->second;
            }
        }

        u64 sid = last_sid;
        u64 soff = atoi(cols[1].first);
        u64 methy_id = (sid<<32) | soff;
        double freq = 1.0 * pcov / (pcov + ncov);
        methy.emplace_back(methy_id, freq);
    }
}

#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>

double s_pearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y) {
    size_t n = x.size();
    if (n != y.size() || n < 2) {
        throw std::invalid_argument("Arrays must be of equal length and size >= 2");
    }

    // 计算均值
    double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;

    // 计算协方差和方差
    double cov_xy = 0.0, var_x = 0.0, var_y = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - mean_x;
        double dy = y[i] - mean_y;
        cov_xy += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    // 处理方差为0的情况（避免除以0）
    if (var_x == 0 || var_y == 0) {
        return 0.0; // 常数组返回0（无相关性）
    }

    // 计算并返回相关系数
    return cov_xy / std::sqrt(var_x * var_y);
}

static int
s_aggregate_methy_correlation(int argc, char* argv[])
{
    HbnAggrMethyCorrOptions options;
    if (!options.parse(argc, argv)) {
        options.dump_usage(argc, argv);
        return 1;
    }
    options.dump_parameters();

    std::map<std::string, int> chr_name2id;
    std::vector<std::pair<uint64_t, double>> m1, m2;
    s_load_bed_methy_call(options.bed1, options.min_cov, chr_name2id, m1);
    s_load_bed_methy_call(options.bed2, options.min_cov, chr_name2id, m2);
    std::pair<uint64_t, double>* a1 = m1.data();
    size_t c1 = m1.size();
    std::pair<uint64_t, double>* a2 = m2.data();
    size_t c2 = m2.size();
    pdqsort(a1, a1 + c1, [](const std::pair<uint64_t, double>& x, const std::pair<uint64_t, double>& y) {
        return x.first < y.first;
    });
    pdqsort(a2, a2 + c2, [](const std::pair<uint64_t, double>& x, const std::pair<uint64_t, double>& y) {
        return x.first < y.first;
    });    
    std::vector<double> v1, v2;
    size_t i1 = 0, i2 = 0;
    while (i1 < c1 && i2 < c2) {
        if (a1[i1].first < a2[i2].first) {
            ++i1;
            continue;
        } else if (a1[i1].first > a2[i2].first) {
            ++i2;
            continue;
        }
        v1.push_back(a1[i1].second);
        v2.push_back(a2[i2].second);
        ++i1;
        ++i2;
    }
    size_t i0 = v1.size();
    if (i0 < 5) {
        HBN_LOG("Intersect genomic loci is less than 5. Skip computation");
        return 0;
    }
    double corr = s_pearsonCorrelation(v1, v2);
    fprintf(stdout, "Intersect loci: %zu\n", i0);
    fprintf(stderr, "correlation: %g\n", corr);

    return 0;
}

} // ns_pileup_correlation

int pileup_correlation_main(int argc, char* argv[])
{
    return ns_pileup_correlation::s_aggregate_methy_correlation(argc, argv);
}