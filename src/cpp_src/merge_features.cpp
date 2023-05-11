#include "kmer_signal.hpp"
#include "5mc_aux.hpp"

#include <pthread.h>

#include <cstring>
#include <limits>
#include <algorithm>

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <errno.h>
#include <pthread.h>
#include <sys/stat.h>

using namespace std;

void dump_merge_usage(int argc, char* argv[])
{
    fprintf(stderr, "\n");
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s [OPTIONS] pos-samples neg-samples merged-samples\n", argv[0]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");
    fprintf(stderr, "  -k <integer>\n");
    fprintf(stderr, "    Kmer size\n");
    fprintf(stderr, "    Default = '%d'\n", default_kmer_size());
    fprintf(stderr, "  -e\n");
    fprintf(stderr, "    Dump the same number of positive and negative samples\n");
}

bool parse_merge_arguments(int argc, char* argv[], 
    int* kmer_size,
    bool* equal_samples,
    const char** pos_samples_path,
    const char** neg_samples_path,
    const char** merged_samples)
{
    if (argc < 4) return false;
    *kmer_size = default_kmer_size();
    *equal_samples = false;
    *pos_samples_path = nullptr;
    *neg_samples_path = nullptr;
    *merged_samples = nullptr;

    int i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') break;
        if (argv[i][0] == '-' && strlen(argv[i]) == 1) break;    

	fprintf(stderr, "%d\t%s\n", i, argv[i]);

        if (strcmp(argv[i], "-k") == 0) {
            *kmer_size = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-e") == 0) {
	    *equal_samples = true;
            i += 1;
            continue;
        }

        fprintf(stderr, "unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i + 3 != argc) return false;
    *pos_samples_path = argv[i];
    *neg_samples_path = argv[i+1];
    *merged_samples = argv[i+2];
    return true;
}

static int
calc_sample_features(int flanking_bases)
{
    return flanking_bases * 2 * 6 + 1;
}

static int 
load_samples_count(const char* samples_path, const int flanking_bases)
{
    struct stat statbuf;
    int ret = stat(samples_path, &statbuf);
    if (ret != 0) {
        fprintf(stderr, "Fail to stat file '%s'\n", samples_path);
        abort();
    }
    int features = calc_sample_features(flanking_bases);
    int res = statbuf.st_size % features;
    if (res != 0) {
        fprintf(stderr, "File size does not match flanking bases (%d)\n", flanking_bases);
        abort();
    }
    return statbuf.st_size / features;
}

int main(int argc, char* argv[])
{
	for (int i = 0; i < argc; ++i) fprintf(stderr, "%d\t%s\n", i, argv[i]);
    int kmer_size = 0;
    bool equal_samples = false;
    const char* pos_samples_path = nullptr;
    const char* neg_samples_path = nullptr;
    const char* merged_samples_path = nullptr;
    if (!parse_merge_arguments(argc, argv, &kmer_size, &equal_samples, &pos_samples_path, &neg_samples_path, &merged_samples_path)) {
        dump_merge_usage(argc, argv);
        return 1;
    }
    int flanking_bases = kmer_size / 2;
    int n_pos = load_samples_count(pos_samples_path, flanking_bases);
    int n_neg = load_samples_count(neg_samples_path, flanking_bases);
    if (equal_samples) {
	    fprintf(stderr, "Dump equal number of positive and negative samples\n");
	    int n = min(n_pos, n_neg);
	    n_pos = n;
	    n_neg = n;
    }

    fprintf(stderr, "Flanking bases: %d\n", flanking_bases);
    fprintf(stderr, "Positive samples: %s\n", pos_samples_path);
    fprintf(stderr, "Negative samples: %s\n", neg_samples_path);
    fprintf(stderr, "Merged samples: %s\n", merged_samples_path);
    fprintf(stderr, "Positive samples: %d\n", n_pos);
    fprintf(stderr, "Negative samples: %d\n", n_neg);


    double n_to_p = (n_pos > 0) ? 1.0 * n_neg / n_pos : 1.0;
    hbn_dfopen(pos_in, pos_samples_path, "rb");
    hbn_dfopen(neg_in, neg_samples_path, "rb");
    hbn_dfopen(out, merged_samples_path, "wb");
    int pos_batch_samples = 1000;
    int neg_batch_samples = pos_batch_samples * n_to_p;
    int features = calc_sample_features(flanking_bases);
    size_t buf_size = (pos_batch_samples + neg_batch_samples) * features;
    char* buf = (char*)malloc(buf_size);

    while (n_pos || n_neg) {
        int xp = min<int>(pos_batch_samples, n_pos);
        int xn = min<int>(neg_batch_samples, n_neg);

        char* s = buf;
        hbn_fread(s, features, xp, pos_in);
        s += features * xp;
        hbn_fread(s, features, xn, neg_in);

        int n = xp + xn;
        hbn_fwrite(buf, features, n, out);

        n_pos -= xp;
        n_neg -= xn;
    }
    hbn_fclose(pos_in);
    hbn_fclose(neg_in);
    hbn_fclose(out);
    free(buf);

    return 0;
}
