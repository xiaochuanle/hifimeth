#ifndef __5MC_CONTEXT_HPP
#define __5MC_CONTEXT_HPP

#include "hbn_aux.hpp"

enum class EMethyContext
{
    eCpG,
    eCHG,
    eCHH,
    eUnknown
};

class MethylationContext
{
public:
    MethylationContext();
    
public:
    static constexpr char FWD_MOD_BASE = 'C';
    static constexpr u8 FWD_MOD_BASE_CODE = 1;
    static constexpr char REV_MOD_BASE = 'G';
    static constexpr u8 REV_MOD_BASE_CODE = 2;

    /// CpG
    static constexpr int CPG_MOTIF_SIZE = 2;
    static constexpr const char* CPG_MOTIF = "CG";

    /// CHG
    static constexpr int CHG_NUM_MOTIFS = 3;
    static constexpr int CHG_MOTIG_SIZE = 3;
    static constexpr int CHG_MAX_HASH = 64;
    static constexpr int CHG_HASH_MASK = 63;
    static constexpr u8  CHG_INVALID_MOTIF_IDX = 255;

    /// CHH
    static constexpr int CHH_NUM_MOTIFS = 9;
    static constexpr int CHH_MOTIF_SIZE = 3;
    static constexpr int CHH_MAX_HASH = 64;
    static constexpr int CHH_HASH_MASK = 63;
    static constexpr u8  CHH_INVALID_MOTIF_IDX = 255;

    static int cpg_opposite_strand_pos(int pos, int seq_size) {
        return seq_size - 1 - (pos + CPG_MOTIF_SIZE - 1);
    }

    static int chg_opposite_strand_pos(int pos, int seq_size) {
        return seq_size - 1 - (pos + CHG_MOTIG_SIZE - 1);
    }

    static int chh_opposite_strand_pos(int pos, int seq_size) {
        return seq_size - 1 - (pos + CHH_MOTIF_SIZE - 1);
    }

public:
    static bool is_cpg_motif_pos(const char* seq) {
        return toupper(seq[0]) == 'C' && toupper(seq[1]) == 'G';
    }

public:
    u8 get_fwd_chg_motif_idx(const char* seq) const {
        u8 hash = extract_motif_hash_value(seq, CHG_MOTIG_SIZE, CHG_MAX_HASH);
        if (hash == CHG_MAX_HASH) return CHG_INVALID_MOTIF_IDX;
        return M_fwd_chg_motif_idx[hash];
    }
    u8 get_rev_chg_motif_idx(const char* seq) const {
        u8 hash = extract_motif_hash_value(seq, CHG_MOTIG_SIZE, CHG_MAX_HASH);
        if (hash == CHG_MAX_HASH) return CHG_INVALID_MOTIF_IDX;
        return M_rev_chg_motif_idx[hash];
    }
    const char* get_fwd_chg_motif(const int motif_idx) const {
        return M_fwd_chg_motifs[motif_idx];
    }
    const char* get_rev_chg_motif(const int motif_idx) const {
        return M_rev_chg_motifs[motif_idx];
    }

    EMethyContext infer_ctx_for_C(const char* s) const {
        hbn_assert(s[0] == 'C');
        if (s[1] == 'G') return EMethyContext::eCpG;
        u8 mi = get_fwd_chg_motif_idx(s);
        if (mi != CHG_INVALID_MOTIF_IDX) return EMethyContext::eCHG;
        mi = get_fwd_chh_motif_idx(s);
        if (mi != CHH_INVALID_MOTIF_IDX) return EMethyContext::eCHH;
        return EMethyContext::eUnknown;
    }

private:
    u8 M_fwd_chg_motif_idx[CHG_MAX_HASH];
    u8 M_rev_chg_motif_idx[CHG_MAX_HASH];
    const char* M_fwd_chg_motifs[CHG_NUM_MOTIFS];
    const char* M_rev_chg_motifs[CHG_NUM_MOTIFS];    

public:
    u8 get_fwd_chh_motif_idx(const char* seq) const {
        u8 hash = extract_motif_hash_value(seq, CHH_MOTIF_SIZE, CHH_MAX_HASH);
        if (hash == CHH_MAX_HASH) return CHH_INVALID_MOTIF_IDX;
        return M_fwd_chh_motif_idx[hash];
    }
    u8 get_rev_chh_motif_idx(const char* seq) const {
        u8 hash = extract_motif_hash_value(seq, CHH_MOTIF_SIZE, CHH_MAX_HASH);
        if (hash == CHH_MAX_HASH) return CHH_INVALID_MOTIF_IDX;
        return M_rev_chh_motif_idx[hash];
    }
    const char* get_fwd_chh_motif(const int motif_idx) const {
        return M_fwd_chh_motifs[motif_idx];
    }
    const char* get_rev_chh_motif(const int motif_idx) const {
        return M_rev_chh_motifs[motif_idx];
    }
private:
    u8 M_fwd_chh_motif_idx[CHH_MAX_HASH];
    u8 M_rev_chh_motif_idx[CHH_MAX_HASH];
    const char* M_fwd_chh_motifs[CHH_NUM_MOTIFS];
    const char* M_rev_chh_motifs[CHH_NUM_MOTIFS];

private:
    static u8 extract_motif_hash_value(const char* seq, const int motif_size, const u8 invalid_hash) {
        u8 hash = 0;
        for (int i = 0; i < motif_size; ++i) {
            u8 c = IUPACNA_TO_BLASTNA[(int)seq[i]];
            if (c > 3) return invalid_hash;
            hash = (hash << 2) | c;
        }
        return hash;
    }
};

#endif // __5MC_CONTEXT_HPP