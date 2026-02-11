#ifndef __BAM_MOD_PARSER_HPP
#define __BAM_MOD_PARSER_HPP

#include "5mc_context.hpp"
#include "bam_info.hpp"

struct BaseModInfo
{
    int qoff;
    uint8_t observed_strand;
    char unmod_base;
    char code;
    uint8_t scaled_prob;
};

void extract_bam_base_mods(bam1_t* bam, std::vector<BaseModInfo>& mods);

#endif // __BAM_MOD_PARSER_HPP