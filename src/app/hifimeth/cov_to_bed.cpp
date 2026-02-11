#include <corelib/5mc_context.hpp>
#include <corelib/hbn_seqdb.hpp>
#include <corelib/line_reader.hpp>

#include <iostream>
#include <string>

#include <cstring>

namespace ns_bismark_call_to_bed {

struct BismarkModInfo
{
    int pcov;
    int ncov;
    const char* motif;
};

static void
s_dump_one_chr_mods(const char* chr_name, const int chr_size, BismarkModInfo* chr_mods, FILE* out)
{
    for (int i = 0; i < chr_size; ++i) {
        if (!chr_mods[i].motif) continue;
        fprintf(out, "%s", chr_name);
        fprintf(out, "\t%d\t%d", i, i + 1);
        int cov = chr_mods[i].pcov + chr_mods[i].ncov;
        hbn_assert(cov > 0);
        double freq = 100.0 * chr_mods[i].pcov / cov;
        fprintf(out, "\t%g", freq);
        fprintf(out, "\t%d\t%d", chr_mods[i].pcov, chr_mods[i].ncov);
        fprintf(out, "\t%s", chr_mods[i].motif);
        fprintf(out, "\n");
    }
}

static int 
s_normalize_cpg_bismark(const char* reference_path, const char* bismark_path, const char* bed_path)
{
    HbnDatabase updb(reference_path);
    hbn_dfopen(out, bed_path, "w");
    BismarkModInfo* chr_mods = nullptr;
    std::string chr, line;
    int last_sid = -1;
    size_t fs = 0, rs = 0;
    HbnLineReader in(bismark_path);
    while (!in.AtEOF()) {
        ++in;
        auto vline = *in;
        line.assign(vline.first, vline.second);
        const char* p = line.c_str();
        int pl = line.size();
        int i = 0;
        int j = 0;

        // chr
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        chr.assign(p+i, j - i);
        int sid = updb.seq_name2id(chr.c_str(), chr.size());
        
        if (sid != last_sid) {
            if (last_sid != -1) {
                s_dump_one_chr_mods(updb.seq_name(last_sid), updb.seq_length(last_sid), chr_mods, out);
                free(chr_mods);
            }
            last_sid = sid;
            int chr_size = updb.seq_length(last_sid);
            chr_mods = (BismarkModInfo*)calloc(chr_size, sizeof(BismarkModInfo));
        }

        // soff
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int soff = atoi(p + i);

        /// send
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int send = atoi(p + i);
        hbn_assert(send - soff == 0, "%s", line.c_str());

        /// freq
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        //double freq = atof(p + i); 

        /// pos-cov
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int pcov = atoi(p + i);       

        /// neg-cov
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int ncov = atoi(p + i);  

        soff -= 1;
        const char* chr_seq = updb.seq_bases(sid);
        int cc0 = chr_seq[soff];

        if (cc0 == MethylationContext::FWD_MOD_BASE) {
            int cc1 = toupper(chr_seq[soff+1]);
            if (cc1 == 'G') {
                BismarkModInfo& mod = chr_mods[soff];
                mod.pcov = pcov;
                mod.ncov = ncov;
                mod.motif = "CG";
                ++fs;
            }
        }
        if (cc0 == MethylationContext::REV_MOD_BASE) {
            int cc1 = chr_seq[soff-1];
            if (cc1 == 'C') {
                BismarkModInfo& mod = chr_mods[soff-1];
                mod.pcov += pcov;
                mod.ncov += ncov;
                mod.motif = "CG";
                ++rs;
            }
        }
    }
    if (chr_mods) {
        s_dump_one_chr_mods(updb.seq_name(last_sid), updb.seq_length(last_sid), chr_mods, out);
        free(chr_mods);
    }
    hbn_fclose(out);
    fprintf(stderr, "forward-strand-sites: %zu, reverse-strand-sites: %zu\n", fs, rs);

    return 0;
}

static int 
s_normalize_chg_bismark(const char* reference_path, const char* bismark_path, const char* bed_path)
{
    HbnDatabase updb(reference_path);
    hbn_dfopen(out, bed_path, "w");
    BismarkModInfo* chr_mods = nullptr;
    std::string chr, line;
    int last_sid = -1;
    size_t fs = 0, rs = 0;
    HbnLineReader in(bismark_path);

    while (!in.AtEOF()) {
	    ++in;
        auto vline = *in;
        line.assign(vline.first, vline.second);
        const char* p = line.c_str();
        int pl = line.size();
        int i = 0;
        int j = 0;

        // chr
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        chr.assign(p+i, j - i);
        int sid = updb.seq_name2id(chr.c_str(), chr.size());
        
        if (sid != last_sid) {
            if (last_sid != -1) {
                s_dump_one_chr_mods(updb.seq_name(last_sid), updb.seq_length(last_sid), chr_mods, out);
                free(chr_mods);
            }
            last_sid = sid;
            int chr_size = updb.seq_length(last_sid);
            chr_mods = (BismarkModInfo*)calloc(chr_size, sizeof(BismarkModInfo));
        }

        // soff
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int soff = atoi(p + i);

        /// send
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int send = atoi(p + i);
        hbn_assert(send - soff == 0, "%s", line.c_str());

        /// freq
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        //double freq = atof(p + i); 

        /// pos-cov
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int pcov = atoi(p + i);       

        /// neg-cov
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int ncov = atoi(p + i);  

        soff -= 1;
        const char* chr_seq = updb.seq_bases(sid);
        int cc0 = chr_seq[soff];

        if (cc0 == MethylationContext::FWD_MOD_BASE) {
            int cc1 = chr_seq[soff+1];
            int cc2 = chr_seq[soff+2];
            if (cc1 == 'C' && cc2 == 'G') {
                BismarkModInfo& mod = chr_mods[soff];
                mod.motif = "CCG";
                mod.pcov = pcov;
                mod.ncov = ncov;
                ++fs;
            }
        }
        if (cc0 == MethylationContext::REV_MOD_BASE) {
            int cc1 = chr_seq[soff-1];
            int cc2 = chr_seq[soff-2];
            if (cc2 == 'C' && cc1 == 'G') {
                BismarkModInfo& mod = chr_mods[soff];
                mod.motif = "CCG";
                mod.pcov = pcov;
                mod.ncov = ncov;
                ++rs;
            }
        }

        if (cc0 == MethylationContext::FWD_MOD_BASE) {
            int cc1 = chr_seq[soff+1];
            int cc2 = chr_seq[soff+2];
            if (cc1 == 'A' && cc2 == 'G') {
                BismarkModInfo& mod = chr_mods[soff];
                mod.motif = "CAG";
                mod.pcov = pcov;
                mod.ncov = ncov;
                ++fs;
            }
        }
        if (cc0 == MethylationContext::REV_MOD_BASE) {
            int cc1 = chr_seq[soff-1];
            int cc2 = chr_seq[soff-2];
            if (cc2 == 'C' && cc1 == 'A') {
                BismarkModInfo& mod = chr_mods[soff-2];
                if (!mod.motif) mod.motif = "CAG";
                mod.pcov += pcov;
                mod.ncov += ncov;
                ++rs;
            }
        }

        if (cc0 == MethylationContext::FWD_MOD_BASE) {
            int cc1 = chr_seq[soff+1];
            int cc2 = chr_seq[soff+2];
            if (cc1 == 'T' && cc2 == 'G') {
                BismarkModInfo& mod = chr_mods[soff];
                mod.motif = "CTG";
                mod.pcov = pcov;
                mod.ncov = ncov;
                ++fs;
            }
        }
        if (cc0 == MethylationContext::REV_MOD_BASE) {
            int cc1 = chr_seq[soff-1];
            int cc2 = chr_seq[soff-2];
            if (cc2 == 'C' && cc1 == 'T') {
                BismarkModInfo& mod = chr_mods[soff-2];
                if (!mod.motif) mod.motif = "CTG";
                mod.pcov += pcov;
                mod.ncov += ncov;
                ++rs;
            }
        }
    }
    if (chr_mods) {
        s_dump_one_chr_mods(updb.seq_name(last_sid), updb.seq_length(last_sid), chr_mods, out);
        free(chr_mods);
    }
    hbn_fclose(out);
    fprintf(stderr, "forward-strand-sites: %zu, reverse-strand-sites: %zu\n", fs, rs);

    return 0;
}

static int 
s_normalize_chh_bismark(const char* reference_path, const char* bismark_path, const char* bed_path)
{
    HbnDatabase updb(reference_path);
    hbn_dfopen(out, bed_path, "w");
    BismarkModInfo* chr_mods = nullptr;
    std::string chr, line;
    int last_sid = -1;
    size_t fs = 0, rs = 0;
    HbnLineReader in(bismark_path);
    MethylationContext ctx;

    while (!in.AtEOF()) {
	    ++in;
        auto vline = *in;
        line.assign(vline.first, vline.second);
        const char* p = line.c_str();
        int pl = line.size();
        int i = 0;
        int j = 0;

        // chr
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        chr.assign(p+i, j - i);
        int sid = updb.seq_name2id(chr.c_str(), chr.size());
        
        if (sid != last_sid) {
            if (last_sid != -1) {
                s_dump_one_chr_mods(updb.seq_name(last_sid), updb.seq_length(last_sid), chr_mods, out);
                free(chr_mods);
            }
            last_sid = sid;
            int chr_size = updb.seq_length(last_sid);
            chr_mods = (BismarkModInfo*)calloc(chr_size, sizeof(BismarkModInfo));
        }

        // soff
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int soff = atoi(p + i);

        /// send
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int send = atoi(p + i);
        hbn_assert(send - soff == 0, "%s", line.c_str());

        /// freq
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        //double freq = atof(p + i); 

        /// pos-cov
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int pcov = atoi(p + i);       

        /// neg-cov
        i = j + 1;
        j = i;
        if (j >= pl) HBN_ERR("corrupted bismark record %s", p);
        while (j < pl && p[j] != '\t') ++j;
        int ncov = atoi(p + i);  

        soff -= 1;
        const char* chr_seq = updb.seq_bases(sid);
        int cc0 = chr_seq[soff];
        if (cc0 == MethylationContext::FWD_MOD_BASE) {
            u8 mi = ctx.get_fwd_chh_motif_idx(chr_seq+soff);
            if (mi != ctx.CHH_INVALID_MOTIF_IDX) {
                BismarkModInfo& mod = chr_mods[soff];
                mod.pcov = pcov;
                mod.ncov = ncov;
                mod.motif = ctx.get_fwd_chh_motif(mi);
                ++fs;
            }
        } else if (cc0 == MethylationContext::REV_MOD_BASE) {
            u8 mi = ctx.get_rev_chh_motif_idx(chr_seq + soff + 1 - ctx.CHH_MOTIF_SIZE);
            if (mi != ctx.CHH_INVALID_MOTIF_IDX) {
                BismarkModInfo& mod = chr_mods[soff];
                mod.pcov = pcov;
                mod.ncov = ncov;
                mod.motif = ctx.get_fwd_chh_motif(mi);
                ++rs;
            }
        }
    }
    if (chr_mods) {
        s_dump_one_chr_mods(updb.seq_name(last_sid), updb.seq_length(last_sid), chr_mods, out);
        free(chr_mods);
    }
    hbn_fclose(out);
    fprintf(stderr, "forward-strand-sites: %zu, reverse-strand-sites: %zu\n", fs, rs);

    return 0;
}

static int
s_bismark_call_to_bed(int argc, char* argv[])
{
    if (argc != 6) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s reference context bismark-call bed\n", argv[0], argv[1]);
        return 1;
    }
    const char* reference_path = argv[2];
    const char* context = argv[3];
    const char* bismark_call_path = argv[4];
    const char* bed_path = argv[5];

    const int ctxlen = strlen(context);
    if (ctxlen == 3) {
        if (strncasecmp(context, "CpG", 3) == 0) {
            return s_normalize_cpg_bismark(reference_path, bismark_call_path, bed_path);
        } else if (strncasecmp(context, "CHG", 3) == 0) {
            return s_normalize_chg_bismark(reference_path, bismark_call_path, bed_path);
        } else if (strncasecmp(context, "CHH", 3) == 0) {
            return s_normalize_chh_bismark(reference_path, bismark_call_path, bed_path);
        }
    }

    std::cerr << "Illegal 5mc context: " << context << '\n'
              << "Plausible contexts: CpG, CHG, CHH\n";
    return 1;
}

} // ns_bismark_call_to_bed

int cov_to_bed_main(int argc, char* argv[])
{
    return ns_bismark_call_to_bed::s_bismark_call_to_bed(argc, argv);
}