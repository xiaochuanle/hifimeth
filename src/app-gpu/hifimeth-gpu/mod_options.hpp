#ifndef __MOD_OPTIONS_HPP
#define __MOD_OPTIONS_HPP

#include <corelib/hbn_aux.hpp>

namespace ns_5mc_call {

struct ModOptions 
{
    int min_read_size;
    int sample_batch_size;
    int read_batch_size;
    bool keep_kinetics;
    bool disable_cpg;
    bool disable_chg;
    bool disable_chh;
    int num_threads;

    std::string model_dir;
    std::string bam_path;
    std::string mod_bam_path;

    ModOptions();
    bool parse(int argc, char* argv[]);
    void dump_usage(int argc, char* argv[]);
    void dump_parameters();
};

struct ModParams
{
    bool disable_cpg;
    int cpg_kmer_size;
    int cpg_features_per_base;
    int cpg_batch_size;

    bool disable_chg;
    int chg_kmer_size;
    int chg_features_per_base;
    int chg_batch_size;

    bool disable_chh;
    int chh_kmer_size;
    int chh_features_per_base;
    int chh_batch_size;
};

} // ns_5mc_call

#endif // __MOD_OPTIONS_HPP