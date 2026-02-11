#include "build_mod_bam.hpp"

#include "../htslib/hts_endian.h"
#include "5mc_context.hpp"

#include <algorithm>
#include <sstream>
#include <string>

#include <errno.h>
#include <unistd.h>

static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end)
{
    int size;
    uint32_t n;
    if (s >= end) return end;
    size = aux_type2size(*s); ++s; // skip type
    switch (size) {
    case 'Z':
    case 'H':
        while (s < end && *s) ++s;
        return s < end ? s + 1 : end;
    case 'B':
        if (end - s < 5) return NULL;
        size = aux_type2size(*s); ++s;
        n = le_to_u32(s);
        s += 4;
        if (size == 0 || end - s < size * n) return NULL;
        return s + size * n;
    case 0:
        return NULL;
    default:
        if (end - s < size) return NULL;
        return s + size;
    }
}

#if 0
static void
iterate_bam_tags(bam1_t* bam, const std::set<int>& skipped_tags)
{
    uint8_t* end = bam->data + bam->l_data;
    uint8_t* s_from = bam_get_aux(bam);
    //uint8_t* s_to = s_from;
    kstring_t tag = KS_INITIALIZE;
    char tagname[2];

    while (s_from < end) {
        //int x = (int)s_from[0]<<8 | s_from[1];
        int t0 = s_from[0];
        int t1 = s_from[1];
        int tx = (t0 << 8) | t1;
        if (skipped_tags.find(tx) == skipped_tags.end()) {
            tagname[0] = t0;
            tagname[1] = t1;
            ks_clear(&tag);
            bam_aux_get_str(bam, tagname, &tag);
            fprintf(stderr, "%s\n", ks_c_str(&tag));
        }
        uint8_t* s = skip_aux(s_from + 2, end);
        s_from = s;
    }

    ks_free(&tag);
}
#endif

static void
s_remove_skipped_tags(bam1_t* bam, const std::set<int>& skipped_tags)
{
    char tagname[2];
    for (auto s : skipped_tags) {
        int t1 = s & 255;
        int t0 = s >> 8;
        tagname[0] = t0;
        tagname[1] = t1;
        uint8_t* aux = bam_aux_get(bam, tagname);
        if (aux) bam_aux_del(bam, aux);
    }

    tagname[0] = 'M';
    tagname[1] = 'L';
    uint8_t* aux = bam_aux_get(bam, tagname);
    if (aux) bam_aux_del(bam, aux);

    tagname[0] = 'M';
    tagname[1] = 'M';
    aux = bam_aux_get(bam, tagname);
    if (aux) bam_aux_del(bam, aux);
}

static void
s_add_one_mm_tag_delta(int delta, std::vector<char>& mmtags)
{
    char buf[64];
    int p = 0;
    do {
        buf[p] = (delta % 10) + '0';
        ++p;
        delta /= 10;
    } while (delta);
    std::reverse(buf, buf + p);
    for (int i = 0; i < p; ++i) mmtags.push_back(buf[i]);
}

void build_one_mod_bam(bam1_t* bam, const std::set<int>& skipped_tags, 
    MolMethyCall* fwd_strand_calls, const int num_fwd_strand_calls, 
    MolMethyCall* rev_strand_calls, const int num_rev_strand_calls)
{
    s_remove_skipped_tags(bam, skipped_tags);
    if (num_fwd_strand_calls == 0 && num_rev_strand_calls == 0) return;
    std::vector<char> mmtags;
    char tagname[2];

    mmtags.push_back('C');
    mmtags.push_back('+');
    mmtags.push_back('m');
    for (int last_qoff = 0, i = 0; i < num_fwd_strand_calls; ++i) {
        if (i < num_fwd_strand_calls - 1) hbn_assert(fwd_strand_calls[i].qoff < fwd_strand_calls[i+1].qoff);
        int c = BamQuerySequence::get_bam_fwd_strand_base(bam, fwd_strand_calls[i].qoff);
        hbn_assert(c == MethylationContext::FWD_MOD_BASE, "i = %c, qoff = %d, c = %c", i, fwd_strand_calls[i].qoff, c);
        int delta = 0;
        for (int k = last_qoff; k < fwd_strand_calls[i].qoff; ++k) {
            c = BamQuerySequence::get_bam_fwd_strand_base(bam, k);
            if (c == MethylationContext::FWD_MOD_BASE) ++delta;
        }
        mmtags.push_back(',');
        s_add_one_mm_tag_delta(delta, mmtags);
        last_qoff = fwd_strand_calls[i].qoff + 1;
    }
    mmtags.push_back(';');

    mmtags.push_back('G');
    mmtags.push_back('-');
    mmtags.push_back('m');
    for (int last_qoff = 0, i = 0; i < num_rev_strand_calls; ++i) {
        if (i < num_rev_strand_calls - 1) hbn_assert(rev_strand_calls[i].qoff < rev_strand_calls[i+1].qoff);
        int c = BamQuerySequence::get_bam_fwd_strand_base(bam, rev_strand_calls[i].qoff);
        hbn_assert(c == MethylationContext::REV_MOD_BASE);
        int delta = 0;
        for (int k = last_qoff; k < rev_strand_calls[i].qoff; ++k) {
            c = BamQuerySequence::get_bam_fwd_strand_base(bam, k);
            if (c == MethylationContext::REV_MOD_BASE) ++delta;
        }
        mmtags.push_back(',');
        s_add_one_mm_tag_delta(delta, mmtags);
        last_qoff = rev_strand_calls[i].qoff + 1;
    }
    mmtags.push_back(';');

    std::vector<uint8_t> mltags;
    for (int i = 0; i < num_fwd_strand_calls; ++i) {
        mltags.push_back(fwd_strand_calls[i].scaled_prob);
    }
    for (int i = 0; i < num_rev_strand_calls; ++i) {
        mltags.push_back(rev_strand_calls[i].scaled_prob);
    }

    tagname[0] = 'M';
    tagname[1] = 'M';
    if (bam_aux_update_str(bam, tagname, mmtags.size(), mmtags.data())) {
        const char* why;
        switch (errno) {
            case EINVAL:
                HBN_ERR("The bam record's aux data is corrupt or an existing tag with the given ID is not of type 'Z'.");
                break;
            case ENOMEM:
                why = "The bam data needs to be expanded and either the attempt to"
                    "reallocate the data buffer failed or the resulting buffer would be"
                    "longer than the maximum size allowed in a bam record (2Gbytes).";
                HBN_ERR("%s", why);
                break;
            default:
                HBN_ERR("Could not add MM tag to BAM record");
                break;
        }
    }

    tagname[0] = 'M';
    tagname[1] = 'L';
    if (bam_aux_update_array(bam, tagname, 'C', mltags.size(), mltags.data())) {
        const char* why;
        switch (errno)
        {
        case EINVAL:
            why = "The bam record's aux data is corrupt, an existing tag with the"
                " given ID is not of an array type or the type parameter is not one of"
                " the values listed above.";
            HBN_ERR("%s", why);
            break;
        case ENOMEM:
            why = "The bam data needs to be expanded and either the attempt to"
                " reallocate the data buffer failed or the resulting buffer would be"
                " longer than the maximum size allowed in a bam record (2Gbytes).";
            HBN_ERR("%s", why);
            break;
        default:
            HBN_ERR("Could not add ML tag to BAM record");
            break;
        }
    }

    tagname[0] = 'M';
    tagname[1] = 'N';
    if (bam_aux_update_int(bam, tagname, bam->core.l_qseq)) {
        const char* why;
        switch (errno) {
            case EINVAL:
                why = "The bam record's aux data is corrupt or an existing tag with the"
                    " given ID is not of an integer type (c, C, s, S, i or I).";
                HBN_ERR("%s", why);
                break;
            case EOVERFLOW:
            case ERANGE:
                why = "val is outside the range that can be stored in an integer bam tag (-2147483647 to 4294967295).";
                HBN_ERR("%s", why);
                break;
            case ENOMEM:
                why = "The bam data needs to be expanded and either the attempt to"
                    " reallocate the data buffer failed or the resulting buffer would be"
                    " longer than the maximum size allowed in a bam record (2Gbytes).";
                HBN_ERR("%s", why);
                break;
            default:
                HBN_ERR("Could not add MN tag to BAM record");
                break;
        }
    }
}
