#ifndef __BUILD_MOD_BAM_HPP
#define __BUILD_MOD_BAM_HPP

#include "5mc_motif_finder.hpp"

#include <set>

void build_one_mod_bam(bam1_t* bam, const std::set<int>& skipped_tags, 
    MolMethyCall* fwd_strand_calls, const int num_fwd_strand_calls, 
    MolMethyCall* rev_strand_calls, const int num_rev_strand_calls);

#endif // __BUILD_MOD_BAM_HPP