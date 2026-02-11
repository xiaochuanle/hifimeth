#ifndef __5MC_MOTIF_FINDER_HPP
#define __5MC_MOTIF_FINDER_HPP

#include "5mc_context.hpp"
#include "bam_info.hpp"
#include "hbn_seqdb.hpp"

struct MolMethyCall
{
    int qid;
    int qoff;
    uint8_t strand;
    uint8_t scaled_prob;
};

struct MappedMolMethyCall
{
    int qid;
    int qoff;
    int sid;
    int soff;
    uint8_t scaled_prob;
};

struct MotifMappedHitInfo
{
    int qoff;
    int sid;
    int soff;
};

void
extract_chg_mapped_samples(HbnDatabase* updb, MethylationContext& chg_motifs, 
    BamQuerySequence& query, BamMapInfo& align, std::vector<MotifMappedHitInfo>& samples);

void
extract_chh_mapped_samples(HbnDatabase* updb, MethylationContext& chh_motifs, 
    BamQuerySequence& query, BamMapInfo& align, std::vector<MotifMappedHitInfo>& samples);

void
extract_cpg_mapped_samples(HbnDatabase* updb,
    BamQuerySequence& query, BamMapInfo& align, std::vector<MotifMappedHitInfo>& samples);

#endif // __5MC_MOTIF_FINDER_HPP