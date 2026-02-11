#include "hbn_seqdb.hpp"

#include "line_reader.hpp"

#include <cstring>

bool s_IsSeqID(const std::string & line)
{
    static const int kMainAccSize = 32;
    size_t digit_pos = line.find_last_of("0123456789|", kMainAccSize);
    if (digit_pos != std::string::npos) {
    	return true;
    }

    return false;
}

static void
s_parse_chr_name_from_header(const char* s, const int sl, const char*& n, int& nl)
{
	int i = 0;
	if (s[0] == '>') {
		i = 1;
		n = s + 1;
	} else {
		i = 0;
		n = s;
	}
	nl = 0;
	for (; i < sl; ++i) {
		if (isspace(s[i])) break;
		++nl;
	}
}

HbnDatabase::HbnDatabase(const char* path)
{
	HbnSeqInfo seqinfo;
	seqinfo.name_offset = 0;
	seqinfo.name_length = 0;
	seqinfo.seq_offset = 0;
	seqinfo.seq_length = 0;
	HbnLineReader in(path);
	while (!in.AtEOF()) {
		++in;
		std::pair<const char*, size_t> full_line = *in;
		std::pair<const char*, size_t> line = truncate_spaces(full_line.first, full_line.second, eTrunc_Both);
		if (line.second == 0) continue;

		char c = line.first[0];
	    if (c == '!'  ||  c == '#' || c == ';') {
	    	continue;
	    }
		bool isId = s_IsSeqID(std::string(line.first, line.second));
	    if ( isId || ( c == '>' )) {
			if (seqinfo.name_length > 0) {
				M_seqinfo_list.push_back(seqinfo);
			}
			seqinfo.name_offset = M_name_list.size();
			seqinfo.name_length = 0;
			seqinfo.seq_offset = M_base_list.size();
			seqinfo.seq_length = 0;

			const char* n = nullptr;
			int nl = 0;
			s_parse_chr_name_from_header(line.first, line.second, n, nl);

			seqinfo.name_length = nl;
			M_name_list.insert(M_name_list.end(), n, n + nl);
			M_name_list.push_back('\0');
	    } else {
			seqinfo.seq_length += line.second;
			M_base_list.insert(M_base_list.end(), line.first, line.first + line.second);
		}
	}
	if (seqinfo.name_length > 0) M_seqinfo_list.push_back(seqinfo);
	for (auto& c : M_base_list) c = toupper(c);

	M_num_seqs = M_seqinfo_list.size();
	M_num_bases = M_base_list.size();

	for (int i = 0; i < M_num_seqs; ++i) {
		const char* name = this->seq_name(i);
		std::pair<const char*, size_t> key(name, strlen(name));
		auto pos = M_name2id.find(key);
		if (pos != M_name2id.end()) {
			std::cerr << "ERROR: Duplicate sequence name " << name << '\n';
			abort();
		}
		M_name2id[key] = i;
	}

	std::string size = bytes_to_datasize(M_num_bases);
	HBN_LOG("Load %d sequences (%s) from %s", M_num_seqs, size.c_str(), path);
}
