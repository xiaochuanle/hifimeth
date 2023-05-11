#include "5mc_aux.hpp"

#include <errno.h>

#include <cstdlib>
#include <cstring>

using namespace std;

static const int kKmerSize = 400;
static const char* kMotif = "CG";

int default_kmer_size()
{
    return kKmerSize;
}

const char* default_motif()
{
    return kMotif;
}

int default_label()
{
    return -1;
}

int reverse_strand_motif_pos(int motif_pos, int motif_size, int seq_size)
{
    return (seq_size - 1) - (motif_pos + motif_size - 1);
}

static inline int
s_string_find_first_not_of(const char* s, const int sl, const char delim, int pos)
{
    for (int i = pos; i < sl; ++i) {
        if (s[i] != delim) return i;
    }
    return sl;
}

static inline int
s_string_find_first_of(const char* s, const int sl, const char delim, int pos)
{
    for (int i = pos; i < sl; ++i) {
        if (s[i] == delim) return i;
    }
    return sl;
}

void
split_string_by_char(const char* s, const int sl, const char delim, vector<pair<const char*, int>>& tokens)
{
    tokens.clear();
    int last_pos = s_string_find_first_not_of(s, sl, delim, 0);
    int pos = s_string_find_first_of(s, sl, delim, last_pos);
    while (last_pos < sl) {
        tokens.push_back(pair<const char*, int>(s + last_pos, pos - last_pos));
        last_pos = s_string_find_first_not_of(s, sl, delim, pos);
        pos = s_string_find_first_of(s, sl, delim, last_pos);
    }
}

char completement_residue(const char c)
{
    switch (c)
    {
    case 'A':
    case 'a':
        return 'T';
        break;
    case 'C':
    case 'c':
        return 'G';
        break;
    case 'G':
    case 'g':
        return 'C';
        break;
    case 'T':
    case 't':
        return 'A';
    
    default:
        return 'A';
        break;
    }    
}