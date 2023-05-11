#ifndef __KMER_SIGNAL_HPP
#define __KMER_SIGNAL_HPP

#include <map>
#include <string>
#include <vector>

class ReadInfo
{
public:
    ReadInfo(const char* motif, int motif_size)
    : m_motif(motif),
      m_motif_size(motif_size) {

    }

    bool parse(const char* sam, const int sam_size, const int min_seq_size);

    void dump() const;

    const int* fwd_ipd_values() const {
        return m_fwd_ipd_values.data();
    }

    const int* rev_ipd_values() const {
        return m_rev_ipd_values.data();
    }

    const int* fwd_pw_values() const {
        return m_fwd_pw_values.data();
    }

    const int* rev_pw_values() const {
        return m_rev_pw_values.data();
    }

    const char* seq_name() const {
        return m_seq_name.c_str();
    }

    const char* fwd_seq() const {
        return m_fwd_seq.c_str();
    }

    const char* rev_seq() const {
        return m_rev_seq.c_str();
    }

    int seq_size() const {
        return m_fwd_seq.size();
    }

    const int* fwd_qual() const {
        return m_fwd_qual.data();
    }

    const int* rev_qual() const {
        return m_rev_qual.data();
    }

    bool extract_features(int offset, 
        int flanking_bases, 
        int label, 
        std::vector<uint8_t>& features, 
        std::vector<std::string>* read_name_list, 
        std::vector<int>* fwd_pos_list,
        std::vector<int>* rev_pos_list,
        bool& is_padded);

private:
    struct TagInfo {
        const char* tag;
        int tag_size;
    };

private:
    const char*                                 m_motif;
    const int                                   m_motif_size;
    const char*                                 m_sam;
    int                                         m_sam_size;
    std::vector<std::pair<const char*, int>>    m_cols;
    std::map<std::string, TagInfo>              m_tags;
    std::string                                 m_seq_name;
    std::string                                 m_fwd_seq;
    std::string                                 m_rev_seq;
    std::vector<int>                            m_fwd_qual;
    std::vector<int>                            m_rev_qual;
    std::vector<int>                            m_fwd_ipd_values;
    std::vector<int>                            m_fwd_pw_values;
    std::vector<int>                            m_rev_ipd_values;
    std::vector<int>                            m_rev_pw_values;
};

#endif // __KMER_SIGNAL_HPP