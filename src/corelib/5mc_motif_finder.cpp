#include "5mc_motif_finder.hpp"

#define USE_CHG_CCG 1
#define USE_CHG_CAG 1
#define USE_CHG_CTG 1

static void
s_extract_ccg_mapped_samples(HbnDatabase* updb, MethylationContext& chg_motifs, 
    BamQuerySequence& query, BamMapInfo& align, std::vector<MotifMappedHitInfo>& samples)
{
    MotifMappedHitInfo msi;
    msi.sid = align.sid;
    const char* chr_seq = updb->seq_bases(align.sid);
    for (int i = 0; i <= align.as_size - MethylationContext::CHG_MOTIG_SIZE; ++i) {
        if (strncasecmp(align.qas + i, align.sas + i, 3) || strncasecmp(align.qas + i, "CCG", 3)) continue;
        if (align.qdir == FWD) {
            msi.qoff = align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::FWD_MOD_BASE_CODE);
        } else {
            msi.qoff = query.size - 1 - align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::REV_MOD_BASE_CODE);
        }
        msi.soff = align.sas_pos[i];
        hbn_assert(toupper(chr_seq[msi.soff]) == MethylationContext::FWD_MOD_BASE);
        samples.push_back(msi);
    }
    for (int i = 0; i <= align.as_size - MethylationContext::CHG_MOTIG_SIZE; ++i) {
        if (strncasecmp(align.qas + i, align.sas + i, 3) || strncasecmp(align.qas + i, "CGG", 3)) continue;
        if (align.qdir == FWD) {
            msi.qoff = align.qas_pos[i] + MethylationContext::CHG_MOTIG_SIZE - 1;
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::REV_MOD_BASE_CODE);
        } else {
            msi.qoff = query.size - 1 - (align.qas_pos[i] + MethylationContext::CHG_MOTIG_SIZE - 1);
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::FWD_MOD_BASE_CODE);
        }
        msi.soff = align.sas_pos[i] + MethylationContext::CHG_MOTIG_SIZE - 1;
        hbn_assert(toupper(chr_seq[msi.soff]) == MethylationContext::REV_MOD_BASE);
        samples.push_back(msi);
    } 
}

static void
s_extract_cag_mapped_samples(HbnDatabase* updb, MethylationContext& chg_motifs, 
    BamQuerySequence& query, BamMapInfo& align, std::vector<MotifMappedHitInfo>& samples)
{
    MotifMappedHitInfo msi;
    msi.sid = align.sid;
    const char* chr_seq = updb->seq_bases(align.sid);
    for (int i = 0; i <= align.as_size - MethylationContext::CHG_MOTIG_SIZE; ++i) {
        if (strncasecmp(align.qas + i, align.sas + i, 3) || strncasecmp(align.qas + i, "CAG", 3)) continue;
        if (align.qdir == FWD) {
            msi.qoff = align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::FWD_MOD_BASE_CODE);
        } else {
            msi.qoff = query.size - 1 - align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::REV_MOD_BASE_CODE);
        }
        msi.soff = align.sas_pos[i];
        hbn_assert(toupper(chr_seq[msi.soff]) == MethylationContext::FWD_MOD_BASE);
        samples.push_back(msi);
    }
}

static void
s_extract_ctg_mapped_samples(HbnDatabase* updb, MethylationContext& chg_motifs, 
    BamQuerySequence& query, BamMapInfo& align, std::vector<MotifMappedHitInfo>& samples)
{
    MotifMappedHitInfo msi;
    msi.sid = align.sid;
    const char* chr_seq = updb->seq_bases(align.sid);
    for (int i = 0; i <= align.as_size - MethylationContext::CHG_MOTIG_SIZE; ++i) {
        if (strncasecmp(align.qas + i, align.sas + i, 3) || strncasecmp(align.qas + i, "CTG", 3)) continue;
        if (align.qdir == FWD) {
            msi.qoff = align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::FWD_MOD_BASE_CODE);
        } else {
            msi.qoff = query.size - 1 - align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::REV_MOD_BASE_CODE);
        }
        msi.soff = align.sas_pos[i];
        hbn_assert(toupper(chr_seq[msi.soff]) == MethylationContext::FWD_MOD_BASE);
        samples.push_back(msi);
    }
}

void
extract_chg_mapped_samples(HbnDatabase* updb, MethylationContext& chg_motifs, 
    BamQuerySequence& query, BamMapInfo& align, std::vector<MotifMappedHitInfo>& samples)
{
    samples.clear();
#if USE_CHG_CCG
    s_extract_ccg_mapped_samples(updb, chg_motifs, query, align, samples);
#endif

#if USE_CHG_CAG
    s_extract_cag_mapped_samples(updb, chg_motifs, query, align, samples);
#endif

#if USE_CHG_CTG
    s_extract_ctg_mapped_samples(updb, chg_motifs, query, align, samples);
#endif
}

void
extract_chh_mapped_samples(HbnDatabase* updb, MethylationContext& chh_motifs, 
    BamQuerySequence& query, BamMapInfo& align, std::vector<MotifMappedHitInfo>& samples)
{
    samples.clear();
    MotifMappedHitInfo msi;
    msi.sid = align.sid;
    const char* chr_seq = updb->seq_bases(align.sid);
    for (int i = 0; i <= align.as_size - MethylationContext::CHH_MOTIF_SIZE; ++i) {
        u8 qmi = chh_motifs.get_fwd_chh_motif_idx(align.qas + i);
        if (qmi == MethylationContext::CHH_INVALID_MOTIF_IDX) continue;
        u8 smi = chh_motifs.get_fwd_chh_motif_idx(align.sas + i);
        if (qmi != smi) continue;
        if (align.qdir == FWD) {
            msi.qoff = align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::FWD_MOD_BASE_CODE);
        } else {
            msi.qoff = query.size - 1 - align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::REV_MOD_BASE_CODE);
        }
        msi.soff = align.sas_pos[i];
        hbn_assert(toupper(chr_seq[msi.soff]) == MethylationContext::FWD_MOD_BASE);
        samples.push_back(msi);
    }
    for (int i = 0; i <= align.as_size - MethylationContext::CHH_MOTIF_SIZE; ++i) {
        u8 qmi = chh_motifs.get_rev_chh_motif_idx(align.qas + i);
        if (qmi == MethylationContext::CHH_INVALID_MOTIF_IDX) continue;
        u8 smi = chh_motifs.get_rev_chh_motif_idx(align.sas + i);
        if (qmi != smi) continue;
        if (align.qdir == FWD) {
            msi.qoff = align.qas_pos[i] + MethylationContext::CHH_MOTIF_SIZE - 1;
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::REV_MOD_BASE_CODE);
        } else {
            msi.qoff = query.size - 1 - (align.qas_pos[i] + MethylationContext::CHH_MOTIF_SIZE - 1);
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::FWD_MOD_BASE_CODE);
        }
        msi.soff = align.sas_pos[i] + MethylationContext::CHH_MOTIF_SIZE - 1;
        hbn_assert(toupper(chr_seq[msi.soff]) == MethylationContext::REV_MOD_BASE);
        samples.push_back(msi);
    }
}

void
extract_cpg_mapped_samples(HbnDatabase* updb,
    BamQuerySequence& query, BamMapInfo& align, std::vector<MotifMappedHitInfo>& samples)
{
    samples.clear();
    MotifMappedHitInfo msi;
    msi.sid = align.sid;
    const char* chr_seq = updb->seq_bases(align.sid);
    for (int i = 0; i <= align.as_size - MethylationContext::CPG_MOTIF_SIZE; ++i) {
        if (strncasecmp(align.qas + i, "CG", 2) || strncasecmp(align.sas + i, "CG", 2)) continue;
        if (align.qdir == FWD) {
            msi.qoff = align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::FWD_MOD_BASE_CODE);
        } else {
            msi.qoff = query.size - 1 - align.qas_pos[i];
            hbn_assert(query.fwd_qs[msi.qoff] == MethylationContext::REV_MOD_BASE_CODE);
        }
        msi.soff = align.sas_pos[i];
        hbn_assert(toupper(chr_seq[msi.soff]) == MethylationContext::FWD_MOD_BASE);
        samples.push_back(msi);
    }
}
