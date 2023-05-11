#include "kmer_signal.hpp"
#include "5mc_aux.hpp"
#include "5mc_tags.hpp"
#include "line_reader.h"

#include <pthread.h>

#include <cstring>
#include <fstream>
#include <limits>
#include <mutex>
#include <sstream>
#include <vector>

using namespace std;

#define BASE_FEATURES 6

static const int kMinReadSize = 500;
static const int kReadBatchSize = 10000;
static const int kMinMapQ = 10;

void dump_features_usage(int argc, char* argv[])
{
    fprintf(stderr, "\n");
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s [OPTIONS] input-sam output [reference]\n", argv[0]);

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
    fprintf(stderr, "  -c <integer>\n");
    fprintf(stderr, "    Number of reads in one chunk\n");
    fprintf(stderr, "    Default = '%d'\n", kReadBatchSize);
    fprintf(stderr, "  -t <integer>\n");
    fprintf(stderr, "    Number of CPU threads used\n");
    fprintf(stderr, "  -q <integer>\n");
    fprintf(stderr, "    Minimum maping quality\n");
    fprintf(stderr, "    Default = '%d'\n", kMinMapQ);
    fprintf(stderr, "\n");
}

bool parse_features_arguments(int argc, char* argv[], 
    int* min_read_size,
    int* kmer_size,
    int* read_batch_size,
    int* min_mapq,
    string& motif,
    int* num_threads,
    const char** input_sam,
    const char** output_dir,
    const char** ref_path)
{
    if (argc < 3) return false;
    *min_read_size = kMinReadSize;
    *kmer_size = default_kmer_size();
    *num_threads = 1;
    *read_batch_size = kReadBatchSize;
    *min_mapq = kMinMapQ;
    motif = default_motif();
    *input_sam = nullptr;
    *output_dir = nullptr;
    *ref_path = nullptr;

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

        if (strcmp(argv[i], "-t") == 0) {
            *num_threads = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-c") == 0) {
            *read_batch_size = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i], "-q") == 0) {
            *min_mapq = atoi(argv[i+1]);
            i += 2;
            continue;
        }

        fprintf(stderr, "unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i + 2 != argc && i + 3 != argc) return false;
    *input_sam = argv[i];
    *output_dir = argv[i+1];
    if (i + 3 == argc) *ref_path = argv[i+2];
    return true;
}

struct BaseFeatures
{
    uint8_t features[BASE_FEATURES];
};

struct SampleOffset
{
    int read_id;
    int fqoff, rqoff;
    int soff;
    int sample_offset;
};

class FeaturesIO
{
public:
    FeaturesIO(const int read_batch_size,
        const int min_read_size,
        const int kmer_size,
        const int min_mapq,
        const char* motif,
        const char* input,
        const char* output_dir,
        Reference* ref) {
        m_read_batch_size = read_batch_size;
        m_min_read_size = min_read_size;
        m_kmer_size = kmer_size;
        m_min_mapq = min_mapq;
        m_motif = motif; for (auto& c : m_motif) c = toupper(c);

        m_in = HbnLineReaderNew(input);
        m_read_id = 0;
        m_next_sam_idx = 0;

        m_output_dir = output_dir;
        m_features_out = NULL;
        m_sample_offset_out = NULL;
        m_sample_info_out = NULL;
        m_features_offset = 0;
        m_chunk_id = 0;

        m_samples = 0;
        m_ref = ref;
    }

    Reference* reference() {
        return m_ref;
    }

    const string& motif() const {
        return m_motif;
    }

    int min_read_size() const {
        return m_min_read_size;
    }

    int kmer_size() const {
        return m_kmer_size;
    }

    int map_quality() const {
        return m_min_mapq;
    }

    ~FeaturesIO() {
        if (m_in) HbnLineReaderFree(m_in);
        if (m_features_out) fclose(m_features_out);
        if (m_sample_offset_out) fclose(m_sample_offset_out);
        if (m_sample_info_out) fclose(m_sample_info_out);

        char path[2048];
        sprintf(path, "%s/num_chunks", m_output_dir);
        hbn_dfopen(out, path, "w");
        fprintf(out, "%d\n", m_chunk_id);
        hbn_fclose(out);

        HBN_LOG("Dump %d reads (%d samples)", m_read_id, m_samples);
    }

    int load_read_batch() {
        if (m_read_id) HBN_LOG("%d reads processed", m_read_id);

        int cnt = 0;
        m_sam_list.clear();
        m_sam_offset_list.clear();
        m_sam_size_list.clear();
        m_next_sam_idx = 0;

        while (!HbnLineReaderAtEof(m_in)) {
            HbnLineReaderReadOneLine(m_in);
            if (ks_empty(m_in->line)) continue;
            if (ks_front(m_in->line) == '@') continue;
            const char* line = ks_s(m_in->line);
            const int line_size = ks_size(m_in->line);
            m_sam_offset_list.push_back(m_sam_list.size());
            m_sam_size_list.push_back(line_size);
            m_sam_list.insert(m_sam_list.end(), line, line + line_size);
            m_sam_list.push_back('\0');
            if (++cnt == m_read_batch_size) break;
        }
        HBN_LOG("Load %d reads", cnt);

        return cnt;
    }

    char* get_next_read(int& read_id, int& sam_size) {
        lock_guard<mutex> lg(m_sam_idx_lock);

        int n = m_sam_offset_list.size();
        if (m_next_sam_idx >= n) return NULL;

	    size_t offset = m_sam_offset_list[m_next_sam_idx];
	    hbn_assert(offset < m_sam_list.size());
        char* s = m_sam_list.data() + m_sam_offset_list[m_next_sam_idx];
        read_id = m_read_id;
        sam_size = m_sam_size_list[m_next_sam_idx];
        ++m_next_sam_idx;
        ++m_read_id;

        return s;
    }

    void open_next_chunk_file() {
        if (m_features_out) fclose(m_features_out);
        if (m_sample_offset_out) fclose(m_sample_offset_out);
        if (m_sample_info_out) fclose(m_sample_info_out);

        char path[2048];
        sprintf(path, "%s/chunk_%d.bin", m_output_dir, m_chunk_id);
        hbn_fopen(m_features_out, path, "w");

        sprintf(path, "%s/chunk_%d.samples", m_output_dir, m_chunk_id);
        hbn_fopen(m_sample_offset_out, path, "w");

        sprintf(path, "%s/chunk_%d.info", m_output_dir, m_chunk_id);
        hbn_fopen(m_sample_info_out, path, "w");

        ++m_chunk_id;
        m_features_offset = 0;
    }

    void dump_features(const int read_idx, 
        const char* read_name,
        const char read_strand,
        const char* ref_name,
        vector<BaseFeatures>& feature_list, 
        vector<SampleOffset>& sample_offset_list) {
        lock_guard<mutex> lg(m_out_lock);

        fwrite(feature_list.data(), sizeof(BaseFeatures), feature_list.size(), m_features_out);
        for (auto& p : sample_offset_list) fprintf(m_sample_offset_out, "%d\t%d\n", read_idx, p.sample_offset + m_features_offset);

        for (auto& p : sample_offset_list) {
            if (p.soff != -1) {
                fprintf(m_sample_info_out, "%s\t%d\t%d\t%d", read_name, read_idx, p.fqoff, p.rqoff);
                fprintf(m_sample_info_out, "\t%s\t%d\t%c\n", ref_name, p.soff, read_strand);
            } else {
                fprintf(m_sample_info_out, "%s\t%d\t%d\t%d\n", read_name, read_idx, p.fqoff, p.rqoff);
            }
        }

        m_features_offset += feature_list.size();
        m_samples += sample_offset_list.size();
    }

private:
    int         m_read_batch_size;
    int         m_min_read_size;
    int         m_kmer_size;
    int         m_min_mapq;
    string      m_motif;

    HbnLineReader*       m_in;
    int         m_read_id;
    vector<char> m_sam_list;
    vector<size_t> m_sam_offset_list;
    vector<int> m_sam_size_list;
    int         m_next_sam_idx;
    mutex       m_sam_idx_lock;

    const char* m_output_dir;
    FILE*       m_features_out;
    FILE*       m_sample_info_out;
    FILE*       m_sample_offset_out;
    mutex       m_out_lock;
    int         m_features_offset;
    int         m_chunk_id;

    int         m_samples;
    Reference*  m_ref;
};

static void*
s_read_features_thread(void* params)
{
    FeaturesIO* data = (FeaturesIO*)(params);
    vector<BaseFeatures> feature_list;
    vector<SampleOffset> sample_offset_list;
    const char* motif = data->motif().c_str();
    const int motif_size = data->motif().size();
    ReadInfo readinfo(motif, motif_size);
    const int kmer_size = data->kmer_size();
    const int min_read_size = data->min_read_size();
    Reference* ref = data->reference();
    const int min_mapq = data->map_quality();
    string ref_name;
    map<int, int> mapinfo;
    char* sam = NULL;
    int read_id = 0, sam_size = 0;
    while ((sam = data->get_next_read(read_id, sam_size))) {
        if (!readinfo.parse(sam, sam_size, min_read_size)) continue;

        const int seq_size = readinfo.seq_size();
        const char* seq = readinfo.fwd_seq();
        const char* seq_name = readinfo.seq_name();
        const int* qual = readinfo.fwd_qual();
        const int* fipd = readinfo.fwd_ipd_values();
        const int* fpw = readinfo.fwd_pw_values();
        const int* ripd = readinfo.rev_ipd_values();
        const int* rpw = readinfo.rev_pw_values();

        feature_list.clear();

        BaseFeatures feature;
        feature.features[0] = nst_nt4_table[(int)seq[0]];
        feature.features[1] = qual[0];
        feature.features[2] = fipd[0];
        feature.features[3] = fpw[0];
        feature.features[4] = ripd[seq_size-1];
        feature.features[5] = rpw[seq_size-1];
        for (int i = 0; i < kmer_size/2; ++i) feature_list.push_back(feature);

        for (int i = 0; i < seq_size; ++i) {
            feature.features[0] = nst_nt4_table[(int)seq[i]];
            feature.features[1] = qual[i];
            feature.features[2] = fipd[i];
            feature.features[3] = fpw[i];
            feature.features[4] = ripd[seq_size-1-i];
            feature.features[5] = rpw[seq_size-1-i];
            feature_list.push_back(feature);            
        }

        feature.features[0] = nst_nt4_table[(int)seq[seq_size-1]];
        feature.features[1] = qual[seq_size-1];
        feature.features[2] = fipd[seq_size-1];
        feature.features[3] = fpw[seq_size-1];
        feature.features[4] = ripd[0];
        feature.features[5] = rpw[0];
        for (int i = 0; i < kmer_size/2; ++i) feature_list.push_back(feature); 

        mapinfo.clear();
        int qdir = FWD;
        extract_mapinfo(sam, sam_size, ref, motif, motif_size, min_mapq, ref_name, qdir, mapinfo); 

        sample_offset_list.clear();
        for (int i = 0; i <= seq_size - motif_size; ++i) {
            if (toupper(seq[i]) != motif[0]) continue;
            bool match_motif = true;
            for (int k = 0; k < motif_size; ++k) {
                if (toupper(seq[i+k]) == motif[k]) continue;
                match_motif = false;
                break;
            }
            if (!match_motif) continue;
            SampleOffset s;
            s.read_id = read_id;
            s.fqoff = i;
            s.rqoff = reverse_strand_motif_pos(i, motif_size, seq_size);
            s.soff = -1;
            auto pos = mapinfo.find(i);
            if (pos != mapinfo.end()) s.soff = pos->second;
            s.sample_offset = kmer_size/2 + i;
            sample_offset_list.push_back(s);
        }      
        if (sample_offset_list.empty()) continue;

        data->dump_features(read_id, seq_name, (qdir == FWD) ? '+' : '-', ref_name.c_str(), feature_list, sample_offset_list);
    }

    return NULL;
}

int main(int argc, char* argv[])
{
    int min_read_size = kMinReadSize;
    int kmer_size = 0;
    string motif;
    int num_threads = 1;
    int read_batch_size = 0;
    int min_mapq = 0;
    const char* input_sam = NULL;
    const char* output_dir = NULL;
    const char* ref_path = NULL;
    if (!parse_features_arguments(argc, argv, &min_read_size, &kmer_size, &read_batch_size, &min_mapq, motif, &num_threads, &input_sam, &output_dir, &ref_path)) {
        dump_features_usage(argc, argv);
        return 1;
    }

    int flanking_bases = kmer_size / 2;
    fprintf(stderr, "\n");
    fprintf(stderr, "Minimum read length: %d\n", min_read_size);
    fprintf(stderr, "Kmer size: %d\n", kmer_size);
    fprintf(stderr, "Flanking bases: %d\n", flanking_bases);
    fprintf(stderr, "Motif: %s\n", motif.c_str());
    fprintf(stderr, "CPU threads used: %d\n", num_threads);
    fprintf(stderr, "Read batch size: %d\n", read_batch_size);
    fprintf(stderr, "Minimum map score: %d\n", min_mapq);
    fprintf(stderr, "Input: %s\n", input_sam);
    fprintf(stderr, "Output dir: %s\n", output_dir);
    fprintf(stderr, "Reference: %s\n", ref_path);
    fprintf(stderr, "\n");

    Reference* ref = ref_path ? new Reference(ref_path) : NULL;
    FeaturesIO data(read_batch_size, min_read_size, kmer_size, min_mapq, motif.c_str(), input_sam, output_dir, ref);
    pthread_t jobs[num_threads];

    while(data.load_read_batch()) {
        data.open_next_chunk_file();
        for (int i = 0; i < num_threads; ++i) {
            pthread_create(jobs + i, NULL, s_read_features_thread, &data);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(jobs[i], NULL);
        }
    }

    return 0;
}
