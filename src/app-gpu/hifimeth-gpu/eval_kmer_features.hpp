#ifndef __EVAL_KMER_FEATURES_HPP
#define __EVAL_KMER_FEATURES_HPP

#include "../../corelib/5mc_context.hpp"
#include "../../corelib/5mc_motif_finder.hpp"
#include "../../corelib/bam_info.hpp"
#include "../../corelib/hbn_seqdb.hpp"

#include <vector>

namespace ns_5mc_call {

class EvalKmerFeaturesGenerator
{
public:
    EvalKmerFeaturesGenerator(int kmer_size)
    {
        M_kmer_size = kmer_size;
        M_features_per_sample = kmer_size * EVAL_BASE_FEATURES;
    }

    ~EvalKmerFeaturesGenerator()
    {
        //delete[] M_one_sample_features;
    }

    bool init(bam1_t* bam) 
    {
        if (!M_query.init(bam)) return false;
        if (!M_kinetics.init(bam)) return false;
        return true;
    }

    void extract_chh_samples();
    void extract_cpg_samples();
    void extract_chg_samples();
    bool get_next_sample_features(float* features, int& offset, int& strand);

public:
    inline static constexpr int EVAL_BASE_FEATURES = 8;
    inline static constexpr int NA_SIZE = 4;
    inline static constexpr int ONE_HOT_ENCODING[NA_SIZE][NA_SIZE] = {
        { 1, 0, 0, 0 },
        { 0, 1, 0, 0 },
        { 0, 0, 1, 0 },
        { 0, 0, 0, 1 }
    };

public:
    int                         M_kmer_size;
    int                         M_features_per_sample;
    MethylationContext          M_ctx;

    BamQuerySequence            M_query;
    BamKinetics                 M_kinetics;

    std::vector<int>            _sample_offsets;
    int*                        M_sample_offsets;
    int                         M_num_samples;
    int                         M_sample_idx;
};

} // ns_5mc_call

#endif // __EVAL_KMER_FEATURES_HPP
