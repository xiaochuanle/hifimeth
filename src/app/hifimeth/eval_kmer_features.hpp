#ifndef __EVAL_KMER_FEATURES_HPP
#define __EVAL_KMER_FEATURES_HPP

#include <corelib/5mc_context.hpp>
#include <corelib/5mc_motif_finder.hpp>
#include <corelib/bam_info.hpp>
#include <corelib/hbn_seqdb.hpp>

#include <vector>

namespace ns_mods {

class EvalKmerFeaturesGenerator
{
public:

    bool init(bam1_t* bam) 
    {
        if (!M_query.init(bam)) return false;
        if (!M_kinetics.init(bam)) return false;
        return true;
    }

    void extract_chh_samples();
    void extract_cpg_samples();
    void extract_chg_samples();
    bool get_next_sample_features(const int kmer_size, const int features_per_base,
        float* features, int& offset, int& strand);

public:
    inline static constexpr int NA_SIZE = 4;
    inline static constexpr int ONE_HOT_ENCODING[NA_SIZE][NA_SIZE] = {
        { 1, 0, 0, 0 },
        { 0, 1, 0, 0 },
        { 0, 0, 1, 0 },
        { 0, 0, 0, 1 }
    };

public:
    MethylationContext          M_ctx;

    BamQuerySequence            M_query;
    BamKinetics                 M_kinetics;

    std::vector<int>            _sample_offsets;
    int*                        M_sample_offsets;
    int                         M_num_samples;
    int                         M_sample_idx;
};

} // ns_mods

#endif // __EVAL_KMER_FEATURES_HPP
