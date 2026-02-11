#include "bam_info.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

#include <cmath>

using namespace std;

static double s_calc_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size)
{
	if (align_size == 0) return 0.0;
	
	int n = 0;
	for (int i = 0; i < align_size; ++i) {
		if (query_mapped_string[i] == target_mapped_string[i]) ++n;
	}

	return 100.0 * n / align_size;
}

static double s_calc_effective_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size)
{
	const int E = 8;
	int eff_len = 0;
	int eff_mat = 0;
	int i = 0;
	while (i < align_size) {
		int qc = query_mapped_string[i];
		int tc = target_mapped_string[i];

		if (qc != GAP_CHAR && tc != GAP_CHAR) {
			if (qc == tc) ++eff_mat;
			++eff_len;
			++i;
			continue;
		}

		if (qc == GAP_CHAR && tc == GAP_CHAR) {
			++i;
			continue;
		}

		if (qc == GAP_CHAR) {
			int j = i + 1;
			while (j < align_size) {
				if (query_mapped_string[j] == GAP_CHAR && target_mapped_string[j] == GAP_CHAR) {
					++j;
					continue;
				}
				if (query_mapped_string[j] != GAP_CHAR) break;
				++j;
			}
			if (j - i < E) {
				for (int k = i; k < j; ++k) {
					int qc = query_mapped_string[k];
					int tc = target_mapped_string[k];
					if (qc == GAP_CHAR && tc == GAP_CHAR) continue;
					if (qc == tc) ++eff_mat;
					++eff_len;
				}
			}
			i = j;
			continue;
		}

		hbn_assert(tc == GAP_CHAR);
		if (tc == GAP_CHAR) {
			int j = i + 1;
			while (j < align_size) {
				if (query_mapped_string[j] == GAP_CHAR && target_mapped_string[j] == GAP_CHAR) {
					++j;
					continue;
				}
				if (target_mapped_string[j] != GAP_CHAR) break;
				++j;
			}
			if (j - i < E) {
				for (int k = i; k < j; ++k) {
					qc = query_mapped_string[k];
					tc = target_mapped_string[k];
					if (qc == GAP_CHAR && tc == GAP_CHAR) continue;
					if (qc == tc) ++eff_mat;
					++eff_len;
				}
			}
			i = j;
			continue;
		}
	}
	if (eff_len == 0) return 0.0;
	return 100.0 * eff_mat / eff_len;
}

static inline char s_decode_bam_query_base(int c)
{
    // 1 -> A
    // 2 -> C
    // 4 -> G
    // 8 -> T
    // 15 -> N
    switch (c) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            HBN_ERR("Illegal BAM base encoded value %s", c);
    }
}

static inline char s_decode_bam_query_complement_base(int c)
{
    // 1 -> A
    // 2 -> C
    // 4 -> G
    // 8 -> T
    // 15 -> N
    switch (c) {
        case 1:
            return 'T'; // 'A';
        case 2:
            return 'G'; // 'C';
        case 4:
            return 'C'; // 'G';
        case 8:
            return 'A'; // 'T';
        case 15:
            return 'N';
        default:
            HBN_ERR("Illegal BAM base encoded value %s", c);
    }
}

static inline char completement_residue(const char c)
{
    switch (c) {
        case 'a':
        case 'A':
            return 'T';
        case 'c':
        case 'C':
            return 'G';
        case 'g':
        case 'G':
            return 'C';
        case 't':
        case 'T':
            return 'A';
        case 'n':
        case 'N':
            return 'N';
        default:
            return 'A';
    }
}

bool BamQuerySequence::init(bam1_t* bam)
{
    //if ((bam->core.flag&256) || (bam->core.flag&2048)) return false; 

    M_name = bam_get_qname(bam);
    size = bam->core.l_qseq;
        
    M_fwd_rqs.clear();
    M_rev_rqs.clear();
    M_fwd_qv.clear();
    M_rev_qv.clear();
    if (bam->core.flag&16) {
        uint8_t* seq = bam_get_seq(bam);
        uint8_t* qual = bam_get_qual(bam);
        for (int i = 0; i < size; ++i) {
            int rc = bam_seqi(seq, i);
            rc = s_decode_bam_query_base(rc);
            int c = completement_residue(rc);
            M_fwd_rqs.push_back(c);
            M_rev_rqs.push_back(rc);
            if (qual) M_rev_qv.push_back(qual[i]+33);
        }        
        reverse(M_fwd_rqs.begin(), M_fwd_rqs.end());
        if (qual) M_fwd_qv.assign(M_rev_qv.rbegin(), M_rev_qv.rend());
    } else {
        uint8_t* seq = bam_get_seq(bam);
        uint8_t* qual = bam_get_qual(bam);
        for (int i = 0; i < size; ++i) {
            int c = bam_seqi(seq, i);
            c = s_decode_bam_query_base(c);
            int rc = completement_residue(c);
            M_fwd_rqs.push_back(c);
            M_rev_rqs.push_back(rc);
            if (qual) M_fwd_qv.push_back(qual[i]+33);
        }
        reverse(M_rev_rqs.begin(), M_rev_rqs.end());
        if (qual) M_rev_qv.assign(M_fwd_qv.rbegin(), M_fwd_qv.rend());
    }

    M_fwd_qs.clear();
    for (auto c : M_fwd_rqs) M_fwd_qs.push_back(IUPACNA_TO_BLASTNA[(int)c]);
    M_rev_qs.clear();
    for (auto c : M_rev_rqs) M_rev_qs.push_back(IUPACNA_TO_BLASTNA[(int)c]);

    name = M_name.c_str();
    fwd_rqs = M_fwd_rqs.data();
    rev_rqs = M_rev_rqs.data();
    fwd_qs = M_fwd_qs.data();
    rev_qs = M_rev_qs.data();
    fwd_qv = M_fwd_qv.empty() ? nullptr : M_fwd_qv.data();
    rev_qv = M_rev_qv.empty() ? nullptr : M_rev_qv.data();

    return true;
}

char BamQuerySequence::get_bam_fwd_strand_base(bam1_t* bam, int fwd_strand_offset)
{
    uint8_t* seq = bam_get_seq(bam);
    return bam_is_rev(bam)
           ?
           s_decode_bam_query_complement_base(bam_seqi(seq, bam->core.l_qseq - 1 - fwd_strand_offset))
           :
           s_decode_bam_query_base(bam_seqi(seq, fwd_strand_offset));
}

char BamQuerySequence::get_bam_rev_strand_base(bam1_t* bam, int rev_strand_offset)
{
    uint8_t* seq = bam_get_seq(bam);
    return bam_is_rev(bam)
           ?
           s_decode_bam_query_base(bam_seqi(seq, rev_strand_offset))
           :
           s_decode_bam_query_complement_base(bam_seqi(seq, bam->core.l_qseq - 1 - rev_strand_offset));
}

//////////////////////

#if 0
CIGAR operation

OP      Description                                                 consume-query           consume-reference

M       alignment match (can be match or mismatch)                  yes                     yes
I       insertion to the reference                                  yes                     no 
D       deletion from the reference                                 no                      yes 
N       skipped region from the reference                           no                      yes
S       soft clipping (clipped sequences represent in SEQ)          yes                     no 
H       hard clipping (clipped sequences NOT represent in SEQ)      no                      no 
P       padding (silent deletion from padded reference)             no                      no 
=       sequence match                                              yes                     yes
X       sequence mismatch                                           yes                     yes
#endif 

static void cigar_to_alignment(const char* query,
    const int query_size,
    const char* subject,
    const int subject_size,
    uint32_t* cigar,
    int n_cigar,
    int& qb,
    int& qe,
    int& sb,
    int& se,
    std::string& qas,
    std::string& sas,
    std::vector<int>& q_pos_list,
    std::vector<int>& s_pos_list)
{
    qas.clear();
    sas.clear();
    q_pos_list.clear();
    s_pos_list.clear();

    char op = bam_cigar_opchr(cigar[0]);
    int num = bam_cigar_oplen(cigar[0]);
    int opi = 0;
    qb = 0;
    sb = 0;
    if (op == 'S') {
        qb = num;
        opi = 1;
    } else if (op == 'H') {
        opi = 1;
    }
    int qi = qb - 1;
    int si = sb - 1;

    qas.clear();
    sas.clear();
    q_pos_list.clear();
    s_pos_list.clear();
    for (; opi < n_cigar; ++opi) {
        op = bam_cigar_opchr(cigar[opi]);
        num = bam_cigar_oplen(cigar[opi]);

        if (op == 'M') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                ++si;
                qas += query[qi];
                sas += subject[si];
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else if (op == 'I') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                qas += query[qi];
                sas += GAP_CHAR;
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else if (op == 'D') {
            for (int i = 0; i < num; ++i) {
                ++si;
                qas += GAP_CHAR;
                sas += subject[si];
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else if (op == 'N') {
            for (int i = 0; i < num; ++i) {
                ++si;
                qas += GAP_CHAR;
                sas += subject[si];
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }            
        } else if (op == 'S') {
            continue;
        } else if (op == 'H') {
            continue;
        } else if (op == 'P') {
            continue;
        } else if (op == '=') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                ++si;
                qas += query[qi];
                sas += subject[si];
		        int qc = query[qi]; qc = toupper(qc);
		        int sc = subject[si]; sc = toupper(sc);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }            
        } else if (op == 'X') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                ++si;
                qas += query[qi];
                sas += subject[si];
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else {
            fprintf(stderr, "ERROR: Unrecognised CIGAR operation '%c' in\n", op);
            //cerr << NStr::CTempString(cigar, cigar_size) << '\n';
            exit(1);
        }
    }
    qe = qi;
    se = si;
}

bool BamMapInfo::init(sam_hdr_t* hdr, bam1_t* bam, HbnDatabase* updb, BamQuerySequence* query)
{
    if (bam->core.flag&0x4) return false;
    //if (bam->core.flag&0x100) return false;
    //if (bam->core.flag&0x800) return false;

    M_sname = hdr->target_name[bam->core.tid];
    M_sid = updb->seq_name2id(M_sname.c_str(), M_sname.size());
    const char* chr_seq = updb->seq_bases(M_sid);
    M_ssize = updb->seq_length(M_sid);
    M_qdir = bam_is_rev(bam);
    const char* query_seq = (M_qdir == FWD) ? query->fwd_rqs : query->rev_rqs;
    M_qsize = query->size;
    cigar_to_alignment(query_seq, M_qsize, chr_seq + bam->core.pos, M_ssize, 
        bam_get_cigar(bam), bam->core.n_cigar, M_qb, M_qe, M_sb, M_se, M_qas, M_sas, M_qas_pos_list, M_sas_pos_list);
    M_as_size = M_qas.size();
    M_sb += bam->core.pos;
    M_se += bam->core.pos;
    for (auto& p : M_sas_pos_list) p += bam->core.pos;
    ++M_qe;
    ++M_se;

    M_mapQ = bam->core.qual;
    M_pi = s_calc_ident_perc(M_qas.c_str(), M_sas.c_str(), M_as_size);
    M_epi = s_calc_effective_ident_perc(M_qas.c_str(), M_sas.c_str(), M_as_size);

    int qi = M_qb - 1, si = M_sb - 1;
    for (int i = 0; i < M_as_size; ++i) {
        if (M_qas[i] != GAP_CHAR) {
            int qc = M_qas[i];
            ++qi;
            hbn_assert(qi == M_qas_pos_list[i]);
            int qc1 = query_seq[qi];
            hbn_assert(qc == qc1);
        }
        if (M_sas[i] != GAP_CHAR) {
            int sc = M_sas[i];
            ++si;
            hbn_assert(si == M_sas_pos_list[i]);
            int sc1 = chr_seq[si];
            //sc1 = DECODE_RESIDUE(sc1);
            hbn_assert(sc == sc1, "sc = %c, sc1 = %c (%d)", sc, sc1, sc1);
        }
    }

    ///////////
    qdir = M_qdir;
    qb = M_qb;
    qe = M_qe;
    qsize = M_qsize;
    sname = M_sname.c_str();
    sid = M_sid;
    sb = M_sb;
    se = M_se;
    ssize = M_ssize;

    mapQ = M_mapQ;
    pi = M_pi;
    epi = M_epi;
    qas = M_qas.c_str();
    sas = M_sas.c_str();
    as_size = M_as_size;
    qas_pos = M_qas_pos_list.data();
    sas_pos = M_sas_pos_list.data();

    return true;
}

////////////////////

static bool
s_extract_kinetic_values(bam1_t* bam, char tagname[], uint8_t*& s, bool& kinetic_is_encoded)
{
    s = bam_aux_get(bam, tagname);
    if (!s) return false;
    if (bam_auxB_len(s) != (uint32_t)bam->core.l_qseq) return false;
    hbn_assert(s[0] == 'B');
    hbn_assert(s[1] == 'C' || s[1] == 'S');
    kinetic_is_encoded = (s[1] == 'C');
    return true;
}

static int s_encode_signal_value(int s)
{
    s = min(952, s);
    int t = 0;
    /// s = (t - 192) * 8 + 448
    /// t = (s - 448) / 8 + 192
    if (s >= 448) {
        t = (s - 448) / 8 + 192;
    }
    /// s = (t - 128) * 4 + 192
    /// t = (s - 192) / 4 + 128
    else if (s >= 192) {
        t = (s - 192) / 4 + 128;
    }
    /// s = (t - 64) * 2 + 64
    /// t = (s - 64) / 2 + 64
    else if (s >= 64) {
        t = (s - 64) / 2 + 64;
    }
    else {
        t = s;
    }
    return t;
}

/*

    CodeV1ToFrameTableSize = 256
    codev1_to_frame_table = [0] * CodeV1ToFrameTableSize
    idx = 0

    for i in range(64):
        codev1_to_frame_table[idx] = i
        idx += 1

    for i in range(64, 128):
        codev1_to_frame_table[idx] = (i - 64) * 2 + 64
        idx += 1

    for i in range(128, 192):
        codev1_to_frame_table[idx] = (i - 128) * 4 + 192
        idx += 1

    for i in range(192, 256):
        codev1_to_frame_table[idx] = (i - 192) * 8 + 448
        idx += 1

    return codev1_to_frame_table
*/

#if 0
static int s_decode_signal_value(int s)
{
    if (s < 64) {
        return s;
    } else if (s < 128) {
        return (s - 64) * 2 + 64;
    } else if (s < 192) {
        return (s - 128) * 4 + 192;
    } else {
        return (s - 192) * 8 + 448;
    }
}
#endif

int BamKinetics::ipd(const int strand, const int offset) const
{
    hbn_assert(strand == FWD || strand == REV);
    hbn_assert(offset >= 0 && offset < M_seq_size);
    int s = 0;
    if (strand == FWD) {
        s = bam_auxB2i(M_fwd_ipd, offset);
        if (!M_fwd_ipd_is_encoded) s = s_encode_signal_value(s);
    } else{
        s = bam_auxB2i(M_rev_ipd, offset);
        if (!M_rev_ipd_is_encoded) s = s_encode_signal_value(s);
    }
    return s;
}

int BamKinetics::pw(const int strand, const int offset) const
{
    hbn_assert(strand == FWD || strand == REV);
    hbn_assert(offset >= 0 && offset < M_seq_size);
    int s = 0;
    if (strand == FWD) {
        s = bam_auxB2i(M_fwd_pw, offset);
        if (!M_fwd_pw_is_encoded) s = s_encode_signal_value(s);
    } else {
        s = bam_auxB2i(M_rev_pw, offset);
        if (!M_rev_pw_is_encoded) s = s_encode_signal_value(s);
    }
    return s;
}

int BamKinetics::decoded_ipd(const int strand, const int offset) const
{
    int s = ipd(strand, offset);
    return M_codev1_table[s];
}

int BamKinetics::decoded_pw(const int strand, const int offset) const
{
    int s = pw(strand, offset);
    return M_codev1_table[s];
}

BamKinetics::BamKinetics()
{
    int p = 0;
    for (int i = 0; i < 64; ++i) M_codev1_table[p++] = i;
    for (int i = 64; i < 128; ++i) M_codev1_table[p++] = (i - 64) * 2 + 64;
    for (int i = 128; i < 192; ++i) M_codev1_table[p++] = (i - 128) * 4 + 192;
    for (int i = 192; i < 256; ++i) M_codev1_table[p++] = (i - 192) * 8 + 448;
    hbn_assert(p == 256);
}

bool BamKinetics::init(bam1_t* bam)
{
    x_clear();

    M_bam = bam;
    M_seq_size = bam->core.l_qseq;

    char kFwdIpdTag[2] = { 'f', 'i' };
    char kRevIpdTag[2] = { 'r', 'i' };
    char kFwdPwTag[2] = { 'f', 'p' };
    char kRevPwTag[2] = { 'r', 'p' };

    if (!s_extract_kinetic_values(bam, kFwdIpdTag, M_fwd_ipd, M_fwd_ipd_is_encoded)) return false;
    if (!s_extract_kinetic_values(bam, kRevIpdTag, M_rev_ipd, M_rev_ipd_is_encoded)) return false;
    if (!s_extract_kinetic_values(bam, kFwdPwTag, M_fwd_pw, M_fwd_pw_is_encoded)) return false;
    if (!s_extract_kinetic_values(bam, kRevPwTag, M_rev_pw, M_rev_pw_is_encoded)) return false;

    uint8_t* fn = bam_aux_get(bam, "fn");
    if (!fn) {
        M_fn = -1;
    } else {
        M_fn = bam_aux2i(fn);
    }
    uint8_t* rn = bam_aux_get(bam, "rn");
    if (!rn) {
        M_rn = -1;
    } else {
        M_rn = bam_aux2i(rn);
    }

    return true;
}
