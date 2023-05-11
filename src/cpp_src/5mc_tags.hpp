#ifndef __5MC_TAGS_HPP
#define __5MC_TAGS_HPP

#include "5mc_aux.hpp"
#include "kmer_signal.hpp"

#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

struct ModCall
{
    int qid;
    int qdir;
    int fqoff, rqoff;
    int sid;
    int soff;
    double prob;
};

struct ReadModCall
{
    int read_id;
    size_t mod_offset;
    int mod_cnt;
    size_t name_offset;
};

class ModCallList
{
public:

    ModCallList(const char* path);

    const char* ref_name(const int ref_id) const {
        return m_ref_name_list.data() + m_ref_name_offset_list[ref_id];
    }

    int add_ref_name(const std::string& ref_name);

    void extract_mod_list(int read_id, ModCall** mca, int* mcc, const char** name);

private:
    std::vector<ModCall>         m_mod_list;
    std::vector<ReadModCall>     m_read_mod_list;
    std::vector<char>            m_read_name_list;
    std::map<int, int>           m_read_idmap;

    std::vector<char>            m_ref_name_list;
    std::vector<size_t>          m_ref_name_offset_list;
    std::map<std::string, int>   m_ref_name2id;
};

class Reference
{
public:

    Reference(const char* ref_path);

    size_t db_size() const {
        return m_seq_list.size();
    }

    int ref_size(const int ref_id) const {
        return m_size_list[ref_id];
    }

    size_t ref_offset(const int ref_id) const {
        return m_offset_list[ref_id];
    }

    const char* ref_seq(const std::string& name) const {
        int id = ref_id(name);
        return m_seq_list.data() + m_offset_list[id];
    }

    int ref_id(const std::string& ref_name) const {
        auto pos = m_name2ids.find(ref_name);
        if (pos == m_name2ids.end()) {
            fprintf(stderr, "Ref name '%s' does not exist\n", ref_name.c_str());
            abort();
        }
        return pos->second;
    }

    const char* ref_name(const int ref_id) {
	    return m_name_list[ref_id].c_str();
    }

    void load_ref_seq(FILE* fa, size_t offset, int ref_size, int res_per_line, int char_per_line);

private:
    std::map<std::string, int> m_name2ids;
    std::vector<std::string> m_name_list;
    std::vector<size_t> m_offset_list;
    std::vector<int> m_size_list;
    std::vector<char> m_seq_list;
};

bool
recover_align_string(std::vector<std::pair<const char*, int>>& cols, const char* subject, std::string& qas, std::string& sas, int& qb, int& qe, int& sb, int& se);

void
extract_mapinfo(const char* sam, const int sam_size, Reference* ref, const char* motif, const int motif_size, const int min_mapq, 
    std::string& ref_name, int& qdir, std::map<int, int>& mapinfo);

void 
extract_forward_read_from_sam(const char* sam_read, const int read_size, const int strand, std::string& fwd_read);

#endif // __5MC_TAGS_HPP