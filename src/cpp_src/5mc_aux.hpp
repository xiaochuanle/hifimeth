#ifndef __5MC_AUX_HPP
#define __5MC_AUX_HPP

#include <fstream>
#include <string>
#include <vector>

#include "hbn_aux.h"

void
split_string_by_char(const char* s, const int sl, const char delim, std::vector<std::pair<const char*, int>>& tokens);

int default_kmer_size();

const char* default_motif();

int default_label();

int reverse_strand_motif_pos(int motif_pos, int motif_size, int seq_size);

char completement_residue(const char c);

#endif // __5MC_AUX_HPP