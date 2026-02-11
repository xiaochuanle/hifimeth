#include "5mc_context.hpp"

static const int kNumCHGMotifs = 3;
static const char* kFwdCHGMotifs[] = { "CAG", "CCG", "CTG" };
static const char* kRevCHGMotifs[] = { "CTG", "CGG", "CAG" };

#if 1
static const int kNumCHHMotifs = 9;
static const char* kFwdCHHMotifs[] = { "CAA", "CCA", "CTA", "CAC", "CCC", "CTC", "CAT", "CCT", "CTT" };
static const char* kRevCHHMotifs[] = { "TTG", "TGG", "TAG", "GTG", "GGG", "GAG", "ATG", "AGG", "AAG" };
#elif 0
static const int kNumCHHMotifs = 1;
static const char* kFwdCHHMotifs[] = { "CAT" };
static const char* kRevCHHMotifs[] = { "ATG" };
#elif 0
static const int kNumCHHMotifs = 1;
static const char* kFwdCHHMotifs[] = { "CCC" };
static const char* kRevCHHMotifs[] = { "GGG" };
#elif 1
static const int kNumCHHMotifs = 4;
static const char* kFwdCHHMotifs[] = { "CAA", "CTA", "CAT", "CTT" };
static const char* kRevCHHMotifs[] = { "TTG", "TAG", "ATG", "AAG" };
#else
static const int kNumCHHMotifs = 4;
static const char* kFwdCHHMotifs[] = { "CAC", "CCC", "CTC", "CCT" };
static const char* kRevCHHMotifs[] = { "GTG", "GGG", "GAG", "AGG" };
#endif

MethylationContext::MethylationContext()
{
    std::fill(M_fwd_chg_motif_idx, M_fwd_chg_motif_idx + CHG_MAX_HASH, CHG_INVALID_MOTIF_IDX);
    std::fill(M_rev_chg_motif_idx, M_rev_chg_motif_idx + CHG_MAX_HASH, CHG_INVALID_MOTIF_IDX);
    for (int i = 0; i < CHG_NUM_MOTIFS; ++i) {
        M_fwd_chg_motifs[i] = kFwdCHGMotifs[i];
        u8 hash = extract_motif_hash_value(M_fwd_chg_motifs[i], CHG_MOTIG_SIZE, CHG_MAX_HASH);
        M_fwd_chg_motif_idx[hash] = i;

        M_rev_chg_motifs[i] = kRevCHGMotifs[i];
        hash = extract_motif_hash_value(M_rev_chg_motifs[i], CHG_MOTIG_SIZE, CHG_MAX_HASH);
        M_rev_chg_motif_idx[hash] = i;
    }

    std::fill(M_fwd_chh_motif_idx, M_fwd_chh_motif_idx + CHH_MAX_HASH, CHH_INVALID_MOTIF_IDX);
    std::fill(M_rev_chh_motif_idx, M_rev_chh_motif_idx + CHH_MAX_HASH, CHH_INVALID_MOTIF_IDX);
    for (int i = 0; i < CHH_NUM_MOTIFS; ++i) {
        M_fwd_chh_motifs[i] = kFwdCHHMotifs[i];
        u8 hash = extract_motif_hash_value(M_fwd_chh_motifs[i], CHH_MOTIF_SIZE, CHH_MAX_HASH);
        M_fwd_chh_motif_idx[hash] = i;

        M_rev_chh_motifs[i] = kRevCHHMotifs[i];
        hash = extract_motif_hash_value(M_rev_chh_motifs[i], CHH_MOTIF_SIZE, CHH_MAX_HASH);
        M_rev_chh_motif_idx[hash] = i;
    }
}