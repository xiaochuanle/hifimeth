#include "kmer_signal.hpp"
#include "5mc_aux.hpp"
#include "5mc_tags.hpp"
#include "line_reader.h"

#include <pthread.h>

#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace std;

static const double kProbCutoff = 0.5;
static const int kMinReadSize = 500;

void dump_features_usage(int argc, char* argv[])
{
    fprintf(stderr, "\n");
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s [OPTIONS] input-sam output\n", argv[0]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");
    fprintf(stderr, "  -l <integer>\n");
    fprintf(stderr, "    Minimum read length\n");
    fprintf(stderr, "    Default = '%d'\n", kMinReadSize);
    fprintf(stderr, "  -k <integer>\n");
    fprintf(stderr, "    kmer size\n");
    fprintf(stderr, "    Default = '%d'\n", default_kmer_size());
    fprintf(stderr, "  -m <string>\n");
    fprintf(stderr, "    Motif\n");
    fprintf(stderr, "    Default = '%s'\n", default_motif());
    fprintf(stderr, "  -b <integer>\n");
    fprintf(stderr, "    Label of samples, only used for training purpose\n");
    fprintf(stderr, "    The default value (%d) means no label will be output\n", default_label());
    fprintf(stderr, "    Default = '%d'\n", default_label());
    fprintf(stderr, "  -t <integer>\n");
    fprintf(stderr, "    Number of CPU threads used\n");
    fprintf(stderr, "  -s <integer>\n");
    fprintf(stderr, "    Skip this number of leading sequences\n");
    fprintf(stderr, "    Default = '%d'\n", 0);
    fprintf(stderr, "  -n <integer>\n");
    fprintf(stderr, "    Extract features for this number of sequences\n");
    fprintf(stderr, "    The default value (-1) means extracting all sequences\n");
    fprintf(stderr, "    Default = '%d'\n", -1);
    fprintf(stderr, "  -x\n");
    fprintf(stderr, "    Dump read name and pos info for each sample\n");
    fprintf(stderr, "  -e <path>\n");
    fprintf(stderr, "    File containing previous 5mc calling results.\n");
    fprintf(stderr, "  -p <probability>\n");
    fprintf(stderr, "    Extracting >= p as label 1 example, and <= 1.0 - p as label 0 example\n");
    fprintf(stderr, "    Default = '%g'\n", kProbCutoff);
    fprintf(stderr, "\n");
}

bool parse_features_arguments(int argc, char* argv[], 
    int* min_read_size,
    int* kmer_size,
    string& motif,
    int* label,
    int* num_threads,
    int* skip_seqs,
    int* num_seqs,
    int* dump_info,
    const char** bed_path,
    double* prob_cutoff,
    const char** input_sam,
    const char** output)
{
    if (argc < 3) return false;
    *min_read_size = kMinReadSize;
    *kmer_size = default_kmer_size();
    *label = default_label();
    *num_threads = 1;
    motif = default_motif();
    *skip_seqs = 0;
    *num_seqs = -1;
    *dump_info = 0;
    *bed_path = nullptr;
    *prob_cutoff = kProbCutoff;
    *input_sam = nullptr;
    *output = nullptr;

    int i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') break;
        if (argv[i][0] == '-' && strlen(argv[i]) == 1) break;

        if (strcmp(argv[i], "-l") == 0) {
            *min_read_size = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-k") == 0) {
            *kmer_size = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-m") == 0) {
            motif = argv[i+1];
            i += 2;
            continue;
        }        

        if (strcmp(argv[i], "-b") == 0) {
            *label = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-t") == 0) {
            *num_threads = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-s") == 0) {
            *skip_seqs = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-n") == 0) {
            *num_seqs = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-x") == 0) {
            *dump_info = 1;
            i += 1;
            continue;
        }

        if (strcmp(argv[i], "-e") == 0) {
            *bed_path = argv[i+1];
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-p") == 0) {
            *prob_cutoff = atof(argv[i+1]);
            i += 2;
            continue;
        }

        fprintf(stderr, "unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i + 2 != argc) return false;
    *input_sam = argv[i];
    *output = argv[i+1];
    return true;
}

struct ThreadData
{
    int target_seqs;
    int seq_idx;
    pthread_mutex_t seq_idx_lock;
    FILE* out;
    FILE* info_out;
    pthread_mutex_t out_lock;
    int processed_seqs;
    int processed_samples;
    int padded_samples;
    int pos_samples;
    int neg_samples;
    const char* motif;
    int min_read_size;
    int flanking_bases;
    int label;
    double prob_cutoff;
    ModCallList* mod_list;
    HbnLineReader* sam_in;
};

static int
s_process_one_read(ReadInfo& readinfo, const char* motif, const int motif_size, const int flanking_size, const int label, 
    vector<uint8_t>& features, int& padded_samples,
    vector<string>* read_name_list = nullptr, vector<int>* fwd_pos_list = nullptr, vector<int>* rev_pos_list = nullptr)
{
    const char* seq = readinfo.fwd_seq();
    const int seq_size = readinfo.seq_size();
    int samples = 0;
    for (int i = 0; i <= seq_size - motif_size; ++i) {
        bool is_motif_loci = false;
        if (seq[i] == motif[0]) {
            is_motif_loci = true;
            for (int k = 1; k < motif_size; ++k) {
                if (seq[i+k] != motif[k]) {
                    is_motif_loci = false;
                    break;
                }
            }
        }
        if (!is_motif_loci) continue;
        bool is_padded = false;
        if (readinfo.extract_features(i, flanking_size, label, features, read_name_list, fwd_pos_list, rev_pos_list, is_padded)) ++samples;
        if (is_padded) ++padded_samples;
    }
    return samples;
}

static int
s_process_one_read_with_bed(ReadInfo& readinfo, const int read_id, const char* motif, const int motif_size, const int flanking_size, const int label,
    vector<uint8_t>& features, int& padded_samples, int& n_pos, int& n_neg,
    ModCallList* mod_list, const double prob_cutoff,
    vector<string>* read_name_list = nullptr, vector<int>* fwd_pos_list = nullptr, vector<int>* rev_pos_list = nullptr)
{
    const char* seq_name = readinfo.seq_name();
    //fprintf(stderr, "%s\n", seq_name);
    ModCall* mca = NULL;
    int mcc = 0;
    const char* name1 = NULL;
    mod_list->extract_mod_list(read_id, &mca, &mcc, &name1);
    if (mcc == 0) return 0;
    hbn_assert(strcmp(seq_name, name1) == 0);

    int samples = 0;
    for (int i = 0; i < mcc; ++i) {
        hbn_assert(mca[i].qid == read_id);
        bool is_padded = false;
        if (readinfo.extract_features(mca[i].fqoff, flanking_size, label, features, read_name_list, fwd_pos_list, rev_pos_list, is_padded)) ++samples;
        if (is_padded) ++padded_samples;           
    }

    return samples;
}

static int
s_read_next_sam_line(ThreadData* data, string& line)
{
    line.clear();
    int read_idx = -1;
    pthread_mutex_lock(&data->seq_idx_lock);
    if (data->seq_idx < data->target_seqs && (!HbnLineReaderAtEof(data->sam_in))) {
        HbnLineReaderReadOneLine(data->sam_in);
        line.assign(ks_s(data->sam_in->line), ks_size(data->sam_in->line));
        read_idx = data->seq_idx;
        ++data->seq_idx;
    }
    pthread_mutex_unlock(&data->seq_idx_lock);
    return read_idx;
}

void* feature_extract_thread(void* params)
{
    ThreadData* data = (ThreadData*)(params);
    ReadInfo readinfo(data->motif, strlen(data->motif));
    vector<uint8_t> features;
    vector<string> _read_name_list;
    vector<int> _fwd_pos_list, _rev_pos_list;
    vector<string>* read_name_list = NULL;
    vector<int>* fwd_pos_list = NULL;
    vector<int>* rev_pos_list = NULL;
    if (data->info_out) {
        read_name_list = &_read_name_list;
        fwd_pos_list = &_fwd_pos_list;
        rev_pos_list = &_rev_pos_list;
    }
    int n_samples = 0, pos_samples = 0, neg_samples = 0;
    int n_seqs = 0;
    int padded_samples = 0;
    string line;

    while (1) {
        int read_idx = s_read_next_sam_line(data, line);
        if (read_idx == -1) break;
        ++n_seqs;
        if (!readinfo.parse(line.c_str(), line.size(), data->min_read_size)) continue;
        if (data->mod_list) {
		    hbn_assert(data->label != -1);
            n_samples += s_process_one_read_with_bed(readinfo, read_idx, data->motif, strlen(data->motif),
                data->flanking_bases, data->label, features, padded_samples, pos_samples, neg_samples, data->mod_list, data->prob_cutoff,
                read_name_list, fwd_pos_list, rev_pos_list);
        } else {
            n_samples += s_process_one_read(readinfo, data->motif, strlen(data->motif),
                 data->flanking_bases, data->label, features, padded_samples,
                read_name_list, fwd_pos_list, rev_pos_list);
        }

        if (features.size() >= 10000) {
            const uint8_t* buf = features.data();
            const size_t buf_size = features.size();
            pthread_mutex_lock(&data->out_lock);
            fwrite(buf, sizeof(uint8_t), buf_size, data->out);
            if (data->info_out) {
                for (size_t i = 0; i < read_name_list->size(); ++i) {
                    fprintf(data->info_out, "%s\t%d\t%d\n", (*read_name_list)[i].c_str(), (*fwd_pos_list)[i], (*rev_pos_list)[i]);
                }
            }
            pthread_mutex_unlock(&data->out_lock);
            features.clear();
            if (read_name_list) read_name_list->clear();
            if (fwd_pos_list) fwd_pos_list->clear();
            if (rev_pos_list) rev_pos_list->clear();
        }
    }

    if (features.size() > 0) {
            const uint8_t* buf = features.data();
            const size_t buf_size = features.size();
            pthread_mutex_lock(&data->out_lock);
            fwrite(buf, sizeof(uint8_t), buf_size, data->out);
            if (data->info_out) {
                for (size_t i = 0; i < read_name_list->size(); ++i) {
                    fprintf(data->info_out, "%s\t%d\t%d\n", (*read_name_list)[i].c_str(), (*fwd_pos_list)[i], (*rev_pos_list)[i]);
                }
            }
            pthread_mutex_unlock(&data->out_lock);
            features.clear();
            if (read_name_list) read_name_list->clear();
            if (fwd_pos_list) fwd_pos_list->clear();
            if (rev_pos_list) rev_pos_list->clear();
    }

    pthread_mutex_lock(&data->out_lock);
    data->processed_samples += n_samples;
    data->processed_seqs += n_seqs;
    data->padded_samples += padded_samples;
    data->pos_samples += pos_samples;
    data->neg_samples += neg_samples;
    pthread_mutex_unlock(&data->out_lock);

    return NULL;
}

static void
s_skip_leading_lines(HbnLineReader* in, int skip_lines)
{
    int cnt = 0;
    while (cnt < skip_lines && (!HbnLineReaderAtEof(in))) {
        HbnLineReaderReadOneLine(in);
        ++cnt;
    }
}

static void
s_skip_sam_header(HbnLineReader* in)
{
    while (!HbnLineReaderAtEof(in)) {
        HbnLineReaderReadOneLine(in);
        if (ks_empty(in->line)) continue;
        if (ks_front(in->line) == '@') continue;
        HbnLineReaderUngetline(in);
        break;
    }
}

int main(int argc, char* argv[])
{
    int min_read_size = 0;
    int kmer_size = 0;
    int label = 0;
    string motif;
    int num_threads = 0;
    int skip_seqs = 0;
    int num_seqs = 0;
    int dump_info = 0;
    double prob_cutoff = 0.0;
    const char* bed_path = nullptr;
    const char* input = nullptr;
    const char* output = nullptr;
    if (!parse_features_arguments(argc, argv, &min_read_size, &kmer_size, motif, &label, &num_threads, &skip_seqs, &num_seqs, &dump_info, &bed_path, &prob_cutoff, &input, &output)) {
        dump_features_usage(argc, argv);
        return 1;
    }
    if (num_seqs < 0) num_seqs = numeric_limits<int>::max();

    int flanking_bases = kmer_size / 2;
    fprintf(stderr, "Minimum read length: %d\n", min_read_size);
    fprintf(stderr, "Kmer size: %d\n", kmer_size);
    fprintf(stderr, "Flanking bases: %d\n", flanking_bases);
    fprintf(stderr, "Motif: %s\n", motif.c_str());
    fprintf(stderr, "Label: %d\n", label);
    fprintf(stderr, "Skipping reads: %d\n", skip_seqs);
    fprintf(stderr, "Extracting reads: %d\n", num_seqs);
    fprintf(stderr, "CPU threads used: %d\n", num_threads);
    fprintf(stderr, "Input: %s\n", input);
    fprintf(stderr, "Output: %s\n", output);
    fprintf(stderr, "Dump sample info: %d\n", dump_info);
    if (bed_path) {
        fprintf(stderr, "BED file: %s\n", bed_path);
        fprintf(stderr, "Label 1 prob cutoff: %g\n", prob_cutoff);
    }
    fprintf(stderr, "\n");

    if (num_seqs != numeric_limits<int>::max()) num_seqs += skip_seqs;

    FILE* info_out = NULL;
    if (dump_info) {
        char path[2048];
        sprintf(path, "%s.sample_info", output);
        hbn_fopen(info_out, path, "w");
    }

    ModCallList* mod_list = NULL;
    if (bed_path) mod_list = new ModCallList(bed_path);

    HbnLineReader* sam_in = HbnLineReaderNew(input);
    s_skip_sam_header(sam_in);
    if (skip_seqs) s_skip_leading_lines(sam_in, skip_seqs);
    hbn_dfopen(out, output, "w");
    ThreadData data;
    data.min_read_size = min_read_size;
    data.flanking_bases = flanking_bases;
    data.label = label;
    data.motif = motif.c_str();
    data.out = out;
    data.info_out = info_out;
    pthread_mutex_init(&data.out_lock, NULL);
    data.processed_samples = 0;
    data.processed_seqs = 0;
    data.padded_samples = 0;
    data.sam_in = sam_in;
    data.seq_idx = skip_seqs;
    pthread_mutex_init(&data.seq_idx_lock, NULL);
    data.target_seqs = num_seqs;
    data.prob_cutoff = prob_cutoff;
    data.mod_list = mod_list;
    data.pos_samples = 0;
    data.neg_samples = 0;

    pthread_t jobids[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobids + i, NULL, feature_extract_thread, &data);
    }    
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobids[i], NULL);
    }
    fprintf(stderr, "Dump %d samples for %d reads\n", data.processed_samples, data.processed_seqs);
    double p = (data.processed_samples) ? 100.0 * data.padded_samples / data.processed_samples : 0.0;
    fprintf(stderr, "Padded samples: %d (%g%%)\n", data.padded_samples, p);
    fprintf(stderr, "Positive samples: %d, negative samples: %d\n", data.pos_samples, data.neg_samples);

    sam_in = HbnLineReaderFree(sam_in);
    fclose(out);
    if (info_out) fclose(info_out);
    fprintf(stderr, "Done\n");

    return 0;
}
