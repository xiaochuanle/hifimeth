#ifndef __HBN_SEQDB_HPP
#define __HBN_SEQDB_HPP

#include "hbn_aux.hpp"
#include "spookyhash.h"

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <cstring>

class HbnDatabase
{
public:
    HbnDatabase(const char* path);

    int num_seqs() const {
        return M_num_seqs;
    }

    size_t num_bases() const {
        return M_num_bases;
    }

    const char* seq_name(int id) const {
        x_validate_target_id(id);
        return M_name_list.data() + M_seqinfo_list[id].name_offset;
    }

    const char* seq_bases(int id) const {
        x_validate_target_id(id);
        return M_base_list.data() + M_seqinfo_list[id].seq_offset;
    }

    int seq_length(int id) const {
        x_validate_target_id(id);
        return M_seqinfo_list[id].seq_length;
    }

    int seq_name2id(const char* name, const int name_length) const {
        auto iter = M_name2id.find(std::pair<const char*, size_t>(name, name_length));
        if (iter == M_name2id.end()) {
            std::cerr << "ERROR: Sequence name " << std::string(name, name_length) << " does not exist" << std::endl;
            abort();
        }
        return iter->second;
    }

private:

    void x_validate_target_id(const int id) const
    {
        if (id < 0 || id >= M_num_seqs) {
            HBN_ERR("Target id %d is out of plausible range [%d, %d)", id, 0, M_num_seqs);
        }
    }

private:
    struct HbnSeqInfo {
        size_t name_offset;
        size_t name_length;
        size_t seq_offset;
        size_t seq_length;
    };

    struct SpookyV2_cstr_to_u64
    {
        u64 operator()(const std::pair<const char*, size_t>& key) const {
            return spookyhash64((const void*)(key.first), key.second, 0);
        }
    };

    struct TempStrEqual 
    {
        bool operator()(const std::pair<const char*, size_t>& x, const std::pair<const char*, size_t>& y) const {
            return (x.second == y.second)
                   &&
                   (strncmp(x.first, y.first, x.second) == 0);
        }
    };

    int                     M_num_seqs;
    size_t                  M_num_bases;
    std::vector<char>       M_name_list;
    std::vector<char>       M_base_list;
    std::vector<HbnSeqInfo> M_seqinfo_list;
    std::unordered_map<std::pair<const char*, size_t>, int, SpookyV2_cstr_to_u64, TempStrEqual>
                            M_name2id;
};

#endif // __HBN_SEQDB_HPP
