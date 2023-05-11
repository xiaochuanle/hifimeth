#include "kmer_signal.hpp"

#include "hbn_aux.h"
#include "5mc_aux.hpp"

using namespace std;

bool ReadInfo::parse(const char* sam, const int sam_size, const int min_seq_size)
{
    bool dump_fail_info = false;
    bool dump_sam_line_when_error = false;

    m_sam = sam;
    m_sam_size = sam_size;
    m_cols.clear();
    split_string_by_char(sam, sam_size, '\t', m_cols);

    m_seq_name.assign(m_cols[0].first, m_cols[0].second);
    int flag = atoi(m_cols[1].first);
    int strand = (flag&16) ? REV : FWD;
    
    if (strand == FWD) {
        m_fwd_seq.assign(m_cols[9].first, m_cols[9].second);
        m_rev_seq.assign(m_fwd_seq.rbegin(), m_fwd_seq.rend());
        for (auto& c : m_rev_seq) c = completement_residue(c);

        m_fwd_qual.clear();
        for (int i = 0; i < m_cols[10].second; ++i) m_fwd_qual.push_back(m_cols[10].first[i]);
        m_rev_qual.assign(m_fwd_qual.rbegin(), m_fwd_qual.rend());
    } else {
        m_rev_seq.assign(m_cols[9].first, m_cols[9].second);
        m_fwd_seq.assign(m_rev_seq.rbegin(), m_rev_seq.rend());
        for (auto& c : m_fwd_seq) c = completement_residue(c);

        m_rev_qual.clear();
        for (int i = 0; i < m_cols[10].second; ++i) m_rev_qual.push_back(m_cols[10].first[i]);
        m_fwd_qual.assign(m_rev_qual.rbegin(), m_rev_qual.rend());
    }

    if (m_fwd_seq.size() < min_seq_size) {
        if (dump_fail_info) fprintf(stderr, "[%s] is too short (%zu bp). Skip it.\n", m_seq_name.c_str(), m_fwd_seq.size());
        return false;
    }
    if (m_fwd_qual.size() == 1) {
        if (dump_fail_info) fprintf(stderr, "[%s] quality value is missing. Skip it.\n", m_seq_name.c_str());
        return false;
    }

    int n_cols = m_cols.size();
    string tag;
    TagInfo taginfo;
    m_tags.clear();
    for (int i = 11; i < n_cols; ++i) {
        tag.assign(m_cols[i].first, 5);
        taginfo.tag = m_cols[i].first + 5;
        taginfo.tag_size = m_cols[i].second - 5;
        m_tags.insert(pair<string, TagInfo>(tag, taginfo));
    }

    vector<pair<const char*, int>> signal_cols;
    tag = "fi:B:";
    auto pos = m_tags.find(tag);
    if (pos == m_tags.end()) {
        if (dump_fail_info) fprintf(stderr, "[%s] forward strand ipd value (tagged by %s) is not provided\n", m_seq_name.c_str(), tag.c_str());
        if (dump_sam_line_when_error) fprintf(stderr, "%s\n", sam);
        return false;
    }
    signal_cols.clear();
    split_string_by_char(pos->second.tag + 2, pos->second.tag_size - 2, ',',  signal_cols);
    m_fwd_ipd_values.clear();
    for (auto& x : signal_cols) {
        int s = atoi(x.first);
        m_fwd_ipd_values.push_back(s);
    }

    tag = "ri:B:";
    pos = m_tags.find(tag);
    if (pos == m_tags.end()) {
        if (dump_fail_info) fprintf(stderr, "[%s] reverse strand ipd value (tagged by %s) is not provided\n", m_seq_name.c_str(), tag.c_str());
        if (dump_sam_line_when_error) fprintf(stderr, "%s\n", sam);
        return false;
    }
    signal_cols.clear();
    split_string_by_char(pos->second.tag + 2, pos->second.tag_size - 2, ',',  signal_cols);
    m_rev_ipd_values.clear();
    for (auto& x : signal_cols) {
        int s = atoi(x.first);
        m_rev_ipd_values.push_back(s);
    }

    tag = "fp:B:";
    pos = m_tags.find(tag);
    if (pos == m_tags.end()) {
        if (dump_fail_info) fprintf(stderr, "[%s] forward pw value (tagged by %s) is not provided\n", m_seq_name.c_str(), tag.c_str());
        if (dump_sam_line_when_error) fprintf(stderr, "%s\n", sam);
        return false;
    }
    signal_cols.clear();
    split_string_by_char(pos->second.tag + 2, pos->second.tag_size - 2, ',', signal_cols);
    m_fwd_pw_values.clear();
    for (auto& x : signal_cols) {
        int s = atoi(x.first);
        m_fwd_pw_values.push_back(s);
    }

    tag = "rp:B:";
    pos = m_tags.find(tag);
    if (pos == m_tags.end()) {
        if (dump_fail_info) fprintf(stderr, "[%s] reverse pw value (tagged by %s) is not provided\n", m_seq_name.c_str(), tag.c_str());
        if (dump_sam_line_when_error) fprintf(stderr, "%s\n", sam);
        return false;
    }
    signal_cols.clear();
    split_string_by_char(pos->second.tag + 2, pos->second.tag_size - 2, ',', signal_cols);
    m_rev_pw_values.clear();
    for (auto& x : signal_cols) {
        int s = atoi(x.first);
        m_rev_pw_values.push_back(s);
    }

    if (m_fwd_seq.size() != m_fwd_qual.size()) {
        if (dump_fail_info) fprintf(stderr, "[%s] sequence length (%zu) and quality length (%zu) do not match\n", m_seq_name.c_str(), m_fwd_seq.size(), m_fwd_qual.size());
        if (dump_sam_line_when_error) fprintf(stderr, "%s\n", sam);
        return false;
    }

    int n_ipd = m_fwd_ipd_values.size();
    int m_seq_size = m_fwd_seq.size();
    if (n_ipd != m_seq_size) {
        if (dump_fail_info) fprintf(stderr, "[%s] forward ipd values (%d) and sequence length (%d) do not match\n", m_seq_name.c_str(), n_ipd, m_seq_size);
        if (dump_sam_line_when_error) fprintf(stderr, "%s\n", sam);
        return false;
    }

    n_ipd = m_rev_ipd_values.size();
    if (n_ipd != m_seq_size) {
        if (dump_fail_info) fprintf(stderr, "[%s] reverse ipd values (%d) and sequence length (%d) do not match\n", m_seq_name.c_str(), n_ipd, m_seq_size);
        if (dump_sam_line_when_error) fprintf(stderr, "%s\n", sam);
        return false;
    }

    int n_pw = m_fwd_pw_values.size();
    if (n_pw != m_seq_size) {
        if (dump_fail_info) fprintf(stderr, "[%s] forward pw values (%d) and sequence length (%d) do not match\n", m_seq_name.c_str(), n_pw, m_seq_size);
        if (dump_sam_line_when_error) fprintf(stderr, "%s\n", sam);
        return false;        
    }

    n_pw = m_rev_pw_values.size();
    if (n_pw != m_seq_size) {
        if (dump_fail_info) fprintf(stderr, "[%s] reverse pw values (%d) and sequence length (%d) do not match\n", m_seq_name.c_str(), n_pw, m_seq_size);
        if (dump_sam_line_when_error) fprintf(stderr, "%s\n", sam);
        return false;        
    }

    return true;
}

void ReadInfo::dump() const {
    fprintf(stderr, "Read name: %s\n", m_seq_name.c_str());
    int m_seq_size = m_fwd_seq.size();

    fprintf(stderr, "Sequence size: %d\n", m_seq_size);

    fprintf(stderr, "fwd ipd values:");
    for (int i = 0; i < 10 && i < m_seq_size; ++i) fprintf(stderr, " %d", m_fwd_ipd_values[i]);
    fprintf(stderr, "...\n");

    fprintf(stderr, "rev ipd values:");
    for (int i = 0; i < 10 && i < m_seq_size; ++i) fprintf(stderr, " %d", m_rev_ipd_values[i]);
    fprintf(stderr, "...\n");

    fprintf(stderr, "fwd pw values:");
    for (int i = 0; i < 10 && i < m_seq_size; ++i) fprintf(stderr, " %d", m_fwd_pw_values[i]);
    fprintf(stderr, "...\n");

    fprintf(stderr, "rev pw values:");
    for (int i = 0; i < 10 && i < m_seq_size; ++i) fprintf(stderr, " %d", m_rev_pw_values[i]);
    fprintf(stderr, "...\n");
}

bool ReadInfo::extract_features(int offset,
    int flanking_bases, 
    int label, 
    std::vector<uint8_t>& features, 
    std::vector<std::string>* read_name_list, 
    std::vector<int>* fwd_pos_list,
    std::vector<int>* rev_pos_list,
    bool& is_padded)
{
    int read_size = seq_size();
    const char* fwd_read = fwd_seq();
    //const char* rev_read = rev_seq();
    const int* fwd_qual_list = fwd_qual();
    //const int* rev_qual_list = rev_qual();
    const int* fwd_ipd_list = fwd_ipd_values();
    const int* rev_ipd_list = rev_ipd_values();
    const int* fwd_pw_list = fwd_pw_values();
    const int* rev_pwd_list = rev_pw_values();

    for (int i = 0; i < m_motif_size; ++i) {
        hbn_assert(fwd_read[offset+i] == m_motif[i]);
    }

    int left_padding = (offset < flanking_bases) ? flanking_bases - offset : 0;
    int right_padding = (offset + flanking_bases > read_size) ? offset + flanking_bases - read_size : 0;
    if (offset < flanking_bases) {
        if (offset < 5) return false;
    }
    if (offset + flanking_bases > read_size) {
	    if (read_size - offset < 5) return false;
    }

    is_padded = false;
    if (left_padding || right_padding) is_padded = true;

    if (read_name_list) read_name_list->push_back(m_seq_name);
    if (fwd_pos_list) fwd_pos_list->push_back(offset);
    if (rev_pos_list) rev_pos_list->push_back(reverse_strand_motif_pos(offset, m_motif_size, read_size));

    int from = max(0, offset - flanking_bases);
    int to = min(read_size, offset + flanking_bases);

    // qual
    int added = 0;
    for (int p = 0; p < left_padding; ++p, ++added) features.push_back(fwd_qual_list[0]);
    for (int p = from; p < to; ++p, ++added) features.push_back(fwd_qual_list[p]);
    for (int p = 0; p < right_padding; ++p, ++added) features.push_back(fwd_qual_list[read_size-1]);
    hbn_assert(added == flanking_bases * 2);

    // fwd seq
    added = 0;
    for (int p = 0; p < left_padding; ++p, ++added) {
        int c = fwd_read[0];
        //c = encode_residue(c);
        c = nst_nt4_table[c];
        features.push_back(c);
    }
    for (int p = from; p < to; ++p, ++added) {
        int c = fwd_read[p];
        //c = encode_residue(c);
        c = nst_nt4_table[c];
        features.push_back(c);
    }
    for (int p = 0; p < right_padding; ++p, ++added) {
        int c = fwd_read[read_size-1];
        //c = encode_residue(c);
        c = nst_nt4_table[c];
        features.push_back(c);
    }
    hbn_assert(added == flanking_bases * 2);

    // fwd ipd
    added = 0;
    for (int p = 0; p < left_padding; ++p, ++added) features.push_back(fwd_ipd_list[0]);
    for (int p = from; p < to; ++p, ++added) features.push_back(fwd_ipd_list[p]);
    for (int p = 0; p < right_padding; ++p, ++added) features.push_back(fwd_ipd_list[read_size-1]);
    hbn_assert(added == flanking_bases * 2);

    // fwd pw 
    added = 0;
    for (int p = 0; p < left_padding; ++p, ++added) features.push_back(fwd_pw_list[0]);
    for (int p = from; p < to; ++p, ++added) features.push_back(fwd_pw_list[p]);
    for (int p = 0; p < right_padding; ++p, ++added) features.push_back(fwd_pw_list[read_size-1]);
    hbn_assert(added == flanking_bases * 2);

    /// label
    if (label != -1) features.push_back(label);

    int x = from, y = to;
    from = read_size - y;
    to = read_size - x;

    // rev ipd
    added = 0;
    for (int p = 0; p < right_padding; ++p, ++added) features.push_back(rev_ipd_list[0]);
    for (int p = from; p < to; ++p, ++added) features.push_back(rev_ipd_list[p]);
    for (int p = 0; p < left_padding; ++p, ++added) features.push_back(rev_ipd_list[read_size-1]);
    hbn_assert(added == flanking_bases * 2);

    // rev pw
    added = 0;
    for (int p = 0; p < right_padding; ++p, ++added) features.push_back(rev_pwd_list[0]);
    for (int p = from; p < to; ++p, ++added) features.push_back(rev_pwd_list[p]);
    for (int p = 0; p < left_padding; ++p, ++added) features.push_back(rev_pwd_list[read_size-1]);
    hbn_assert(added == flanking_bases * 2);

    return true;
}
