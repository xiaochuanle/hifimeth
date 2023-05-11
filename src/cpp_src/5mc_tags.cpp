#include "5mc_tags.hpp"

#include "line_reader.h"

#include <algorithm>

using namespace std;

ModCallList::ModCallList(const char* path) {
        HBN_LOG("Load mod calls from %s\n", path);
        HbnLineReader* in = HbnLineReaderNew(path);
        std::vector<std::pair<const char*, int>> cols;
        ModCall mc;
        int last_qid = -1;
        std::string last_ref_name, ref_name;
        int last_ref_id = -1;
        while (!HbnLineReaderAtEof(in)) {
            HbnLineReaderReadOneLine(in);
            const char* line = ks_s(in->line);
            const int line_size = ks_size(in->line);
            cols.clear();
            split_string_by_char(line, line_size, '\t', cols);
            mc.qid = atoi(cols[1].first);
            mc.fqoff = atoi(cols[2].first);
            mc.rqoff = atoi(cols[3].first);
            mc.sid = -1;
            mc.soff = -1;
            mc.prob = atof(cols.back().first);

            if (cols.size() == 8) {
                ref_name.assign(cols[4].first, cols[4].second);
                if (ref_name != last_ref_name) {
                    last_ref_name = ref_name;
                    last_ref_id = add_ref_name(last_ref_name);
                }
                mc.sid = last_ref_id;
                mc.soff = atoi(cols[5].first);
                mc.qdir = cols[6].first[0] == '+' ? FWD : REV;
            }

            m_mod_list.push_back(mc);
            if (mc.qid != last_qid) {
                last_qid = mc.qid;
                m_read_name_list.insert(m_read_name_list.end(), cols[0].first, cols[0].first + cols[0].second);
                m_read_name_list.push_back('\0');
            }
        }
        in = HbnLineReaderFree(in);

        int new_qid = 0;
        size_t n_mod = m_mod_list.size();
        size_t i = 0;
        size_t name_offset = 0;
        while (i < n_mod) {
            size_t j = i + 1;
            while (j < n_mod && m_mod_list[i].qid == m_mod_list[j].qid) ++j;
            ReadModCall rmc;
            rmc.read_id = m_mod_list[i].qid;
            rmc.mod_offset = i;
            rmc.mod_cnt = j - i;
            rmc.name_offset = name_offset;
            while (name_offset < m_read_name_list.size() && m_read_name_list[name_offset] != '\0') ++name_offset;
            hbn_assert(name_offset < m_read_name_list.size());
            ++name_offset;
            m_read_mod_list.push_back(rmc);
	        auto pos = m_read_idmap.find(rmc.read_id);
	        hbn_assert(pos == m_read_idmap.end());
            m_read_idmap[rmc.read_id] = new_qid;
            ++new_qid;
            i = j;
        }

        fprintf(stderr, "Load %zu mod calls in %d reads\n", n_mod, new_qid);
}

int ModCallList::add_ref_name(const std::string& ref_name) {
        auto pos = m_ref_name2id.find(ref_name);
        if (pos != m_ref_name2id.end()) return pos->second;

        int ref_id = m_ref_name_offset_list.size();
        m_ref_name_offset_list.push_back(m_ref_name_list.size());
        m_ref_name_list.insert(m_ref_name_list.end(), ref_name.begin(), ref_name.end());
        m_ref_name_list.push_back('\0');
        m_ref_name2id.insert(std::pair<std::string, int>(ref_name, ref_id));
        return ref_id;
}

void ModCallList::extract_mod_list(int read_id, ModCall** mca, int* mcc, const char** name) {
        *mca = NULL;
        *mcc = 0;
        auto pos = m_read_idmap.find(read_id);
        if (pos == m_read_idmap.end()) return;

        int nid = pos->second;
        *mca = m_mod_list.data() + m_read_mod_list[nid].mod_offset;
        *mcc = m_read_mod_list[nid].mod_cnt;
        if (name) *name = m_read_name_list.data() + m_read_mod_list[nid].name_offset;
}

Reference::Reference(const char* ref_path) {
        size_t buf_size = 2048;
        char* buf = (char*)malloc(buf_size);
        sprintf(buf, "%s.fai", ref_path);
        HbnLineReader* fai = HbnLineReaderNew(buf);
        hbn_dfopen(fa, ref_path, "r");
        std::vector<std::pair<const char*, int>> cols;
        std::string ref_name;
        int ref_id = 0;
        size_t res_loaded = 0;
        while (!HbnLineReaderAtEof(fai)) {
            HbnLineReaderReadOneLine(fai);
            const char* line = ks_s(fai->line);
            const int line_size = ks_size(fai->line);
            cols.clear();
            split_string_by_char(line, line_size, '\t', cols);
            ref_name.assign(cols[0].first, cols[0].second);
            m_name2ids.insert(std::pair<std::string, int>(ref_name, ref_id));
            m_name_list.push_back(ref_name);
            ++ref_id;

            int ref_size = atoi(cols[1].first);
            size_t offset = atoll(cols[2].first);
            int res_per_line = atoi(cols[3].first);
            int char_per_line = atoi(cols[4].first);

            m_offset_list.push_back(res_loaded);
            m_size_list.push_back(ref_size);
            load_ref_seq(fa, offset, ref_size, res_per_line, char_per_line);
            res_loaded += ref_size;
        }
        fclose(fa);
        HbnLineReaderFree(fai);
        free(buf);
}

void Reference::load_ref_seq(FILE* fa, size_t offset, int ref_size, int res_per_line, int char_per_line) {
        fseek(fa, offset, SEEK_SET);
        char* line = new char[char_per_line];
        int left = ref_size;
        int eol = char_per_line - res_per_line;
        while (left > 0) {
            int n = std::min<int>(left, res_per_line);
            hbn_fread(line, 1, n+eol, fa);
            m_seq_list.insert(m_seq_list.end(), line, line + n);
            left -= n;
        }
        delete[] line;
}

bool
recover_align_string(std::vector<std::pair<const char*, int>>& cols, const char* subject, std::string& qas, std::string& sas, int& qb, int& qe, int& sb, int& se)
{
    const char* cigar = cols[5].first;
    const int cigar_size = cols[5].second;
	for (int i = 0; i < cigar_size; ++i) {
		if (cigar[i] == 'H') {
			fprintf(stderr, "'");
			fwrite(cigar, 1, cigar_size, stderr);
			fprintf(stderr, "' contains hard clip, skip this record\n");
			return false;
		}
		if (isdigit(cigar[i])) continue;
		if (cigar[i] != 'M' && cigar[i] != 'I' && cigar[i] != 'D' && cigar[i] != 'S' && cigar[i] != '=' && cigar[i] != 'X') {
			fprintf(stderr, "'");
			fwrite(cigar, 1, cigar_size, stderr);
			fprintf(stderr, "' contains unrecognised operation '%c', skip this record\n", cigar[i]);
			return false;
		}
	}

    qb = 0, qe = 0;
    int i = 0;
    while (i < cigar_size && isdigit(cigar[i])) ++i;
    hbn_assert(i < cigar_size);
    if (cigar[i] == 'S' || cigar[i] == 'H') {
        qb = atoi(cigar);
        ++i;
    } else {
        i = 0;
    }
    qe = qb;
    const char* query = cols[9].first;
    const int query_size = cols[9].second;

    sb = atoi(cols[3].first) - 1;
    se = sb;

    qas.clear();
    sas.clear();
    while (i < cigar_size) {
        hbn_assert(isdigit(cigar[i]));
        int num = atoi(cigar + i);
        while (i < cigar_size && isdigit(cigar[i])) ++i;
        hbn_assert(i < cigar_size);
        char op = cigar[i];
        ++i;
        if (op == '=' || op == 'X') {
            for (int k = 0; k < num; ++k) {
		        if (qe + k >= query_size) {
	                fwrite(cigar, 1, cigar_size, stderr);
	                fprintf(stderr, "\n");
			        fprintf(stderr, "i = %d, cigar_size = %d, num = %d, k = %d, qe = %d, query_size = %d\n", i, cigar_size, num, k, qe, query_size);
		        }
		        hbn_assert(qe + k < query_size);
                qas += query[qe + k];
                sas += subject[se+k];
            }
            qe += num;
            se += num;
        } else if (op == 'I') {
            for (int k = 0; k < num; ++k) {
		        hbn_assert(qe + k < query_size);
                qas += query[qe + k];
                sas += GAP_CHAR;
            }
            qe += num;
        } else if (op == 'D') {
            for (int k = 0; k < num; ++k) {
                qas += GAP_CHAR;
                sas += subject[se + k];
            }
            se += num;
        } else if (op == 'S' || op == 'H') {
            continue;
        } else {
            fprintf(stderr, "Illegal CIGAR operation '%c'\n", op);
	        fwrite(cigar, 1, cigar_size, stderr);
	        fprintf(stderr, "\n");
            abort();
        }
    }
    return true;
}

void
extract_mapinfo(const char* sam, const int sam_size, Reference* ref, const char* motif, const int motif_size, const int min_mapq, 
    std::string& ref_name, int& qdir, std::map<int, int>& mapinfo)
{
    if (!ref) return;
    vector<pair<const char*, int>> cols;
    split_string_by_char(sam, sam_size, '\t', cols);
    if (cols[5].first[0] == '*') return;
    if (cols[9].first[0] == '*') return;
    int mapq = atoi(cols[4].first);
    if (mapq < min_mapq) return;

    ref_name.assign(cols[2].first, cols[2].second);
    int flag = atoi(cols[1].first);
    qdir = (flag&16) ? REV : FWD;
    const char* subject = ref->ref_seq(ref_name);
    const char* query = cols[9].first;
    const int query_size = cols[9].second;

    int qb, qe, sb, se;
    string qas, sas;
    recover_align_string(cols, subject, qas, sas, qb, qe, sb, se);

    int qi = qb;
    int si = sb;
    int as_i = 0, as_size = qas.size();
    for (as_i = 0; as_i < as_size - motif_size; ++as_i) {
        if (qas[as_i] != '-') {
            if (query[qi] != qas[as_i]) {
                fprintf(stderr, "query bases do not match\n");
                abort();
            }
        }
        if (sas[as_i] != '-') {
            if (sas[as_i] != subject[si]) {
                fprintf(stderr, "ref bases do not match\n");
                abort();
            }
        }
                      
        if (toupper(qas[as_i]) == motif[0]) {
            bool is_match = true;
            for (int k = 0; k < motif_size; ++k) {
                if (toupper(qas[as_i+k]) != motif[k]) {
                    is_match = false;
                    break;
                }
                if (toupper(sas[as_i+k]) != motif[k]) {
                    is_match = false;
                    break;
                }
            }      
            if (is_match) {
                int qoff = qi;
                if (qdir == REV) qoff = reverse_strand_motif_pos(qoff, motif_size, query_size);
                mapinfo[qoff] = si;
            }
        }
        if (qas[as_i] != '-') ++qi;
        if (sas[as_i] != '-') ++si;
    }
}

void 
extract_forward_read_from_sam(const char* sam_read, const int read_size, const int strand, string& fwd_read)
{
    if (strand == FWD) {
        fwd_read.assign(sam_read, read_size);
        for (auto& c : fwd_read) c = toupper(c);
        return;
    }

    fwd_read.assign(sam_read, read_size);
    reverse(fwd_read.begin(), fwd_read.end());
    for (auto& c : fwd_read) {
        int xc = c;
        xc = nst_nt4_table[xc];
        xc = 3 - xc;
        c = "ACGT"[xc];
    }
}