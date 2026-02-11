#include "../../corelib/bam_info.hpp"
#include "../../corelib/hbn_aux.hpp"

#include <algorithm>

namespace ns_unmmapped_bam_extraction_random {

struct UnmappedBamInfo
{
    int id;
    bool is_valid;
    bool is_selected;
    int length;
};

static void
s_set_bam_info(bam1_t* bam, BamQuerySequence& query, BamKinetics& kinetics, UnmappedBamInfo& ubi)
{
    ubi.is_valid = false;
    ubi.is_selected = false;
    if (bam->core.l_qseq < 5000) return;
    if (!query.init(bam)) return;
    if (!kinetics.init(bam)) return;

    ubi.is_valid = true;
    ubi.length = bam->core.l_qseq;
}

static int
s_extract_unmapped_bam_random(int argc, char* argv[])
{
    if (argc != 6) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "  %s %s reference input-bam coverage output-bam\n", argv[0], argv[1]);
        return 1;
    }

    const char* reference_path = argv[2];
    const char* input_bam_path = argv[3];
    const int cov = atoi(argv[4]);
    const char* output_bam_path = argv[5];

    HbnDatabase updb(reference_path);
    const size_t dbsize = updb.num_bases();
    const size_t target = dbsize * cov;

    samFile* in = sam_open(input_bam_path, "rb");
    hts_set_threads(in, 8);
    sam_hdr_t* hdr = sam_hdr_read(in);
    bam1_t* bam = bam_init1();

    samFile* out = sam_open(output_bam_path, "wb");
    hts_set_threads(out, 8);
    if (sam_hdr_write(out, hdr)) HBN_ERR("Could not write BAM header");

    BamQuerySequence query;
    BamKinetics kinetics;
    std::vector<UnmappedBamInfo> bam_list;
    size_t id = 0;
    while (1) {
        int r = sam_read1(in, hdr, bam);
        if (r < 0) break;
        UnmappedBamInfo ubi;
        s_set_bam_info(bam, query, kinetics, ubi);
        ubi.id = id++;
        if ((id % 100000) == 0) HBN_LOG("%10zu reads processed", id);
        bam_list.push_back(ubi);
    }
    sam_hdr_destroy(hdr);
    sam_close(in);
    in = sam_open(input_bam_path, "rb");
    hts_set_threads(in, 8);
    hdr = sam_hdr_read(in);

    size_t total_bases = 0;
    for (auto& ubi : bam_list) if (ubi.is_valid) total_bases += ubi.length;
    std::string size = bytes_to_datasize(dbsize);
    fprintf(stderr, "DB size: %s\n", size.c_str());
    size = bytes_to_datasize(target);
    fprintf(stderr, "coverage: %d, target size: %s\n", cov, size.c_str());
    size = bytes_to_datasize(total_bases);
    fprintf(stderr, "BAM size: %s\n", size.c_str());

    RandomFloatNumberGenerator rng;
    int extracted_reads = 0;
    size_t extracted_bases = 0;
    {
        std::shuffle(bam_list.begin(), bam_list.end(), rng.gen());
        std::shuffle(bam_list.begin(), bam_list.end(), rng.gen());
        size_t num_bases = 0;
        for (auto& ubi : bam_list) {
            if (!ubi.is_valid) continue;
            num_bases += ubi.length;
            ubi.is_selected = true;
            if (num_bases >= target) break;
        }
        std::sort(bam_list.begin(), bam_list.end(), [](const UnmappedBamInfo& x, const UnmappedBamInfo& y) { return x.id < y.id; });
        id = 0;
        while (1) {
            int r = sam_read1(in, hdr, bam);
            if (r < 0) break;
            UnmappedBamInfo& ubi = bam_list[id++];
            if ((id % 100000) == 0) HBN_LOG("%10zu reads processed", id);
            if (!ubi.is_valid) continue;
            if (!ubi.is_selected) continue;
            if (sam_write1(out, hdr, bam) < 0) HBN_ERR("Could not write BAM record to file");
            ++extracted_reads;
            extracted_bases += ubi.length;
        }
    }

    sam_hdr_destroy(hdr);
    sam_close(in);
    bam_destroy1(bam);
    sam_close(out);

    size = bytes_to_datasize(target);
    fprintf(stderr, "Target: %s\n", size.c_str());
    size = bytes_to_datasize(extracted_bases);
    fprintf(stderr, "Extracted reads: %d (%s)\n", extracted_reads, size.c_str());

    return 0;
}

} // ns_unmmapped_bam_extraction_random

int subsample_bam_main(int argc, char* argv[])
{
    return ns_unmmapped_bam_extraction_random::s_extract_unmapped_bam_random(argc, argv);
}