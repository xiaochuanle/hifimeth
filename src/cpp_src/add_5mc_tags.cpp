#include "5mc_aux.hpp"
#include "kmer_signal.hpp"
#include "5mc_tags.hpp"
#include "line_reader.h"

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

static void
s_dump_one_sam(vector<pair<const char*, int>>& sam_cols, string& mm_tag, string& ml_tag, set<string>& skipped_tags, FILE* out)
{
    fwrite(sam_cols[0].first, 1, sam_cols[0].second, out);
    for (int i = 1; i < 11; ++i) {
        fprintf(out, "\t");
        fwrite(sam_cols[i].first, 1, sam_cols[i].second, out);
    }

    string tag;
    int n_col = sam_cols.size();
    for (int i = 11; i < n_col; ++i) {
        tag.assign(sam_cols[i].first, 2);
        if (skipped_tags.find(tag) != skipped_tags.end()) continue;
        fprintf(out, "\t");
        fwrite(sam_cols[i].first, 1, sam_cols[i].second, out);
    }

    fprintf(out, "\t");
    fprintf(out, "%s", mm_tag.c_str());
    fprintf(out, "\t");
    fprintf(out, "%s", ml_tag.c_str());
    fprintf(out, "\n");
}

int main(int argc, char* argv[])
{
    if (argc != 5) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s motif 5mc-path input-sam output-sam\n", argv[0]);
        return 1;
    }

    string motif = argv[1]; for (auto& c : motif) c = toupper(c);
    const char* mod_path = argv[2];
    const char* input_sam = argv[3];
    const char* output_sam = argv[4];

    fprintf(stderr, "Motif: %s\n", motif.c_str());
    fprintf(stderr, "5mc path: %s\n", mod_path);
    fprintf(stderr, "Input sam: %s\n", input_sam);
    fprintf(stderr, "Output sam: %s\n", output_sam);

    set<string> skipped_tags;
    skipped_tags.insert(string("fi"));
    skipped_tags.insert(string("ri"));
    skipped_tags.insert(string("fp"));
    skipped_tags.insert(string("rp"));
    skipped_tags.insert(string("MM"));
    skipped_tags.insert(string("ML"));

    ModCallList mod_list(mod_path);
    hbn_dfopen(out, output_sam, "w");
    HbnLineReader* in = HbnLineReaderNew(input_sam);
    string read_name, fwd_read;
    vector<pair<const char*, int>> cols;
    ostringstream mmos, mlos;
    string mmtag, mltag;
    int read_idx = -1;

    while (!HbnLineReaderAtEof(in)) {
        HbnLineReaderReadOneLine(in);
        const char* line = ks_s(in->line);
        const int line_size = ks_size(in->line);
        if (line_size == 0) continue;
        if (line[0] == '@') {
            hbn_fwrite(line, 1, line_size, out);
            fprintf(out, "\n");
            continue;
        }
        ++read_idx;
        mmos.str("");
        mmos << "MM:Z:" << motif[0] << "+m";
        mlos.str("");
        mlos << "ML:B:C";

        cols.clear();
        split_string_by_char(line, line_size, '\t', cols);
        ModCall* mca = NULL;
        int mcc = 0;
        const char* read_name1 = NULL;
        mod_list.extract_mod_list(read_idx, &mca, &mcc, &read_name1);
        if (mcc == 0) {
            mmos << ';';
            mmtag = mmos.str();
            mltag = mlos.str();
            s_dump_one_sam(cols, mmtag, mltag, skipped_tags, out);
            continue;
        }

        read_name.assign(cols[0].first, cols[0].second);
        hbn_assert(read_name == read_name1);
        const int read_size = cols[9].second;
        int sam_flag = atoi(cols[1].first);
        int strand = (sam_flag&16) ? REV : FWD;

        extract_forward_read_from_sam(cols[9].first, read_size, strand, fwd_read);
        const char* read = fwd_read.c_str();

        int last_offset = 0;
        for (int i = 0; i < mcc; ++i) {
            hbn_assert(mca[i].qid == read_idx);
	        for (int k = 0; k < motif.size(); ++k) {
		        int p = mca[i].fqoff + k;
		        int c = toupper(read[p]);
		        if (c != motif[k]) { fprintf(stderr, "fwd %s, %d, %c %c, fqoff = %d, rqoff = %d, i = %d\n", read_name.c_str(), p, c, motif[k], mca[i].fqoff, mca[i].rqoff, k); }
		        hbn_assert(c == motif[k]);
	        }

            int delta = 0;
            for (int p = last_offset; p < mca[i].fqoff; ++p) {
                hbn_assert(p < read_size);
                int c = toupper(read[p]);
                if (c == motif[0]) ++delta;
            }
            mmos << ',' << delta;
            int prob = mca[i].prob * 256;
            if (prob > 255) prob = 255;
            mlos << ',' << prob;
            last_offset = mca[i].fqoff + 1;
        }
        mmos << ';';
        mmtag = mmos.str();
        mltag = mlos.str();
        s_dump_one_sam(cols, mmtag, mltag, skipped_tags, out);
    }

    hbn_fclose(out);
    HbnLineReaderFree(in);

    return 0;
}
