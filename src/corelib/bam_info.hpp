#ifndef __BAM_INFO_HPP
#define __BAM_INFO_HPP

#include "../htslib/sam.h"
#include "hbn_aux.hpp"
#include "hbn_seqdb.hpp"

#include <string>
#include <vector>

struct BamQuerySequence
{
public:
    bool init(bam1_t* bam);

    static char get_bam_fwd_strand_base(bam1_t* bam, int fwd_strand_offset);
    static char get_bam_rev_strand_base(bam1_t* bam, int rev_strand_offset);

public:
    int         id;
    const char* name;
    const char* fwd_rqs;
    const char* rev_rqs;
    const u8*   fwd_qs;
    const u8*   rev_qs;
    const char* fwd_qv;
    const char* rev_qv;
    int         size;

public:
    std::string         M_name;
    std::vector<char>   M_fwd_rqs;
    std::vector<char>   M_rev_rqs;
    std::vector<char>   M_fwd_qv;
    std::vector<char>   M_rev_qv;
    std::vector<u8>     M_fwd_qs;
    std::vector<u8>     M_rev_qs;
    int                 M_fn;
    int                 M_rn;
};

struct BamMapInfo
{
public:
    bool init(sam_hdr_t* hdr, bam1_t* bam, HbnDatabase* updb, BamQuerySequence* query);

public:
    int qdir, qb, qe, qsize;
    const char* sname;
    int sid, sb, se, ssize;
    int mapQ;
    double pi, epi;
    const char* qas;
    const char* sas;
    int as_size;
    const int* qas_pos;
    const int* sas_pos;

private:

    /// mapping info
    int                                         M_qdir;
    int                                         M_qb;
    int                                         M_qe;
    int                                         M_qsize;
    int                                         M_sid;
    std::string                                 M_sname;
    int                                         M_sb;
    int                                         M_se;
    int                                         M_ssize;
    int                                         M_mapQ;
    double                                      M_pi;
    double                                      M_epi;
    std::string                                 M_qas;
    std::string                                 M_sas;
    int                                         M_as_size;
    std::vector<int>                            M_qas_pos_list;
    std::vector<int>                            M_sas_pos_list;
};

///////////////

class BamKinetics
{
public:
    BamKinetics();

    bool init(bam1_t* bam);

    int ipd(const int strand, const int offset) const;
    int pw(const int strand, const int offset) const;

    int decoded_ipd(const int strand, const int offset) const;
    int decoded_pw(const int strand, const int offset) const;

    const int* codev1_table() const {
        return M_codev1_table;
    }

    int fn() const {
        return M_fn;
    }
    int rn() const {
        return M_rn;
    }

public:
    static constexpr int MAX_KINETIC_VALUE = 952;

private:
    void x_clear()
    {
        M_bam = nullptr;
        M_seq_size = 0;
        M_fwd_ipd = nullptr;
        M_rev_ipd = nullptr;
        M_fwd_pw = nullptr;
        M_rev_pw = nullptr;
    }

private:
    bam1_t*             M_bam;
    int                 M_seq_size;
    uint8_t*            M_fwd_ipd;
    bool                M_fwd_ipd_is_encoded;
    uint8_t*            M_rev_ipd;
    bool                M_rev_ipd_is_encoded;
    uint8_t*            M_fwd_pw;
    bool                M_fwd_pw_is_encoded;
    uint8_t*            M_rev_pw;
    bool                M_rev_pw_is_encoded;
    int                 M_fn;
    int                 M_rn;

    int                 M_codev1_table[256];
};

#endif // __BAM_INFO_HPP
