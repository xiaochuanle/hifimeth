#include "bam_mod_parser.hpp"

#include <iostream>
#include <sstream>

#include <errno.h>

static void
s_extract_bam_mod_scaled_probs(bam1_t* bam, std::vector<uint8_t>& probs)
{
    const char tagname[] = {'M', 'L'};
    uint8_t* mltag = bam_aux_get(bam, tagname);
    if (!mltag) {
        if (errno == ENOENT) return;
        if (errno == EINVAL) {
            std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                      << "Could not extract ML aux tag from BAM record\n"
                      << "The bam record's aux data is corrupt (either a tag has an"
                      << "invalid type, or the last record is incomplete)\n";
            abort();
        }
    }
    uint32_t N = bam_auxB_len(mltag);
    for (uint32_t i = 0; i != N; ++i) {
        int64_t v = bam_auxB2i(mltag, i);
        if (v < 0 || v > 255) {
            std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                      << "Illegal scaled probability value " << v 
                      << ", which must be in range [0, 255]\n";
            abort();
        }
        probs.push_back(static_cast<uint8_t>(v));
    }
}

static char
s_chebi_to_iupac_code(bam1_t* bam, const char* s, int& si, const int sl)
{
    int c = 0;
    for (; si < sl; ++si) {
        if (!isdigit(s[si])) break;
        c = c * 10 + s[si] - '0';
    }
    if (s[si] >= sl || s[si] != ',') {
        std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                  << "Illegal edit series using EhEBI code: " << std::string(s, sl) << '\n';
        abort(); 
    }
    ++si;
    if (c == 27551) {
        return 'm';
    } else if (c == 76792) {
        return 'h';
    } else if (c == 76794) {
        return 'f';
    } else if (c == 76793) {
        return 'c';
    } else if (c == 16964) {
        return 'g';
    } else if (c == 80961) {
        return 'e';
    } else if (c ==  17477) {
        return 'b';
    } else if (c == 28871) {
        return 'a';
    } else if (c == 44605) {
        return 'o';
    } else if (c == 18107) {
        return 'n';
    }

    std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
              << "Unrecognised ChEBI code (" << c << ") in edit series " << std::string(s, sl) << '\n';
    abort();
    return '\0';
}

static void
s_extract_unmod_base_and_observed_strand(bam1_t* bam, const char* s, const int sl, char& unmod_base, int& strand)
{
    unmod_base = s[0];
    if (unmod_base != 'C' && unmod_base != 'G' && unmod_base != 'T' && unmod_base != 'A'
        && unmod_base != 'U' && unmod_base != 'N') {
        std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                  << "Unrecognosed unmodified base (" << unmod_base << ") in edit series " << std::string(s, sl) << '\n'
                  << "Plausible unmodified bases: {C, G, T, A, U, N}\n";
        abort();
    }
    if (s[1] != '+' && s[1] != '-') {
        std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                  << "Unrecognosed strand (" << s[1] << ") in edit series " << std::string(s, sl) << '\n'
                  << "Plausible strands: {+, -}\n";
        abort();        
    }
    strand = (s[1] == '+') ? FWD : REV;
}

static bool
s_is_valid_unmod_base_and_code(const char unmod_base, const std::string& codes)
{
    for (char c : codes) {
        if (c == 'm' || c == 'h' || c == 'f' || c == 'c' || c == 'C') {
            if (unmod_base != 'C' && unmod_base != 'G') {
                return false;
            }
        }
        if (c == 'g' || c == 'e' || c == 'b' || c == 'T') {
            if (unmod_base != 'T' && unmod_base != 'A') {
                return false;
            }
        }
        if (c == 'U') {
            if (unmod_base != 'U') {
                return false;
            }
        }
        if (c == 'a' || c == 'A') {
            if (unmod_base != 'A' && unmod_base != 'T') {
                return false;
            }
        }
        if (c == 'o' || c == 'G') {
            if (unmod_base != 'G' && unmod_base != 'C') {
                return false;
            }
        }
        if (c =='n' || c == 'N') {
            if (unmod_base != 'N') {
                return false;
            }
        }
    }
    return true;
}

static void
s_parse_one_mod_list(bam1_t* bam, const char* s, const int sl, const uint8_t* probs, const int num_probs, int& prob_idx,
    std::vector<BaseModInfo>& mods)
{
    if (sl < 4 || s[sl-1] != ';') {
        std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                  << "Corrupted edit series " << std::string(s, sl) << '\n';
        abort();
    }

    char unmod_base;
    int strand;
    s_extract_unmod_base_and_observed_strand(bam, s, sl, unmod_base, strand);
    std::string codes;
    int num_code = 0;
    int si = 2;
    if (isdigit(s[2])) {
        codes += s_chebi_to_iupac_code(bam, s, si, sl);
        num_code = 1;
    } else {
        for (; si < sl; ++si) {
            if (s[si] == ',' || s[si] == ';') break;
	    if (s[si] != '.' && s[si] != '?' ) {
            	codes += s[si];
            	++num_code;
	    }
        }
    }

    if (!s_is_valid_unmod_base_and_code(unmod_base, codes)) {
        std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
            << "In edit series " << std::string(s, sl) << '\n'
            << "Inconsistent combination of unmodified base (" << unmod_base << ") and modification code (" << codes << ")\n";
        abort();
    }

    for (int i = si; i < sl; ++i) {
        if (isdigit(s[si])) continue;
        if (s[si] != ',' && s[si] != ';') {
            std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                      << "Illegal character " << s[si] << " in edit series " << std::string(s, sl) << '\n';
            abort();
        }
    }

    std::vector<int> edit_list;
    hbn_assert(s[si] == ',' || s[si] == ';');
    ++si;
    while (si < sl) {
        hbn_assert(isdigit(s[si]));
        int delta = 0;
        while (si < sl) {
            if (s[si] == ',' || s[si] == ';') break;
            delta = delta * 10 + s[si] - '0';
            ++si;
        }
        edit_list.push_back(delta);
        hbn_assert(s[si] == ',' || s[si] == ';');
        ++si;
    }

    int qoff = 0;
    for (auto d : edit_list) {
	    //fprintf(stderr, "d = %d, strand = %d\n", d, strand);
        int cnt = 0;
        while (cnt < d) {
            hbn_assert(qoff < bam->core.l_qseq);
            int c = BamQuerySequence::get_bam_fwd_strand_base(bam, qoff);
            if (c == unmod_base) ++cnt;
            ++qoff;
        }
        hbn_assert(cnt == d);
        hbn_assert(qoff < bam->core.l_qseq, "qoff = %d, qlen = %d", qoff, bam->core.l_qseq);

        while (1) {
            hbn_assert(qoff < bam->core.l_qseq);
            int c = BamQuerySequence::get_bam_fwd_strand_base(bam, qoff);
            if (c == unmod_base) break;
	        ++qoff;
        }

        for (int j = 0; j < num_code; ++j) {
            BaseModInfo mod;
            mod.qoff = qoff;
            mod.observed_strand = strand;
            mod.unmod_base = unmod_base;
            mod.code = codes[j];
            hbn_assert(prob_idx < num_probs);
            mod.scaled_prob = probs[prob_idx++];
            mods.push_back(mod);
        }
        ++qoff;
    }
}

void extract_bam_base_mods(bam1_t* bam, std::vector<BaseModInfo>& mods)
{
    char tagname[] = { 'M', 'N' };
    uint8_t* tag;
#if 0
    uint8_t* tag = bam_aux_get(bam, tagname);
    if (tag) {
        int mn_qlen = bam_aux2i(tag);
        if (mn_qlen != bam->core.l_qseq) {
            std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                      << "Read length (" << mn_qlen << ", represented by the MN BAM tag) is when base modification tags"
                      << " (MM, ML, MN) are created is different from the current read length ("
                      << bam->core.l_qseq << ") of the SEQ field.\n";
            abort();
        }
    }
#endif

    std::vector<uint8_t> scaled_probs;
    s_extract_bam_mod_scaled_probs(bam, scaled_probs);
    if (scaled_probs.empty()) return;
    const int num_probs = scaled_probs.size();
    int prob_idx = 0;

    tagname[0] = 'M';
    tagname[1] = 'M';
    tag = bam_aux_get(bam, tagname);
    if (!tag) {
        if (errno == ENOENT) return;
        if (errno == EINVAL) {
            std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                      << "Could not extract MM aux tag from BAM record\n"
                      << "The bam record's aux data is corrupt (either a tag has an"
                      << "invalid type, or the last record is incomplete)\n";
            abort();
        }
    }
    const char* mms = bam_aux2Z(tag);
    const int mmsl = strlen(mms);
    if (mms[mmsl-1] != ';') {
        std::cerr << "ERROR at parsing read " << bam_get_qname(bam) << '\n'
                  << "The MM aux tag must end with ';'\n"
                  << std::string(mms, mmsl) << '\n';
        abort();
    }
    int i = 0;
    while (i < mmsl) {
        int j = i + 1;
        while (j < mmsl && mms[j] != ';') ++j;
	    hbn_assert (j < mmsl && mms[j] == ';');
	    ++j;
	    //std::cerr << "parse edit list " << std::string_view(mms + i, j - i) << '\n';
        s_parse_one_mod_list(bam, mms + i, j - i, scaled_probs.data(), num_probs, prob_idx, mods);
        i = j;
    }
}
