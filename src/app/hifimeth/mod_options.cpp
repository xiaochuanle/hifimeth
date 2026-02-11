#include "mod_options.hpp"

#include <corelib/arg_parse.hpp>
#include <corelib/get_core_count.hpp>

#include "program_info.hpp"

namespace ns_mods {

static constexpr int kMinReadSize = 1000;
static constexpr int kSampleBatchSize = 32;
static constexpr int kReadBatchSize = 10000;
static constexpr bool kKeepKinetics = false;
static constexpr bool kDisableCpG = false;
static constexpr bool kDisableCHG = false;
static constexpr bool kDisableCHH = false;
static constexpr int kNumThreads = 0;

ModOptions::ModOptions()
{
    min_read_size = kMinReadSize;
    sample_batch_size = kSampleBatchSize;
    read_batch_size = kReadBatchSize;
    keep_kinetics = kKeepKinetics;
    disable_cpg = kDisableCpG;
    disable_chg = kDisableCHG;
    disable_chh = kDisableCHH;
    num_threads = kNumThreads;
}

static bool
s_parse_context_string(const char* option, const char* argument,
    bool& skip_cpg, bool& skip_chg, bool& skip_chh)
{
    skip_cpg = skip_chg = skip_chh = true;
    std::string ctx;
    size_t N = strlen(argument);
    size_t i = 0;
    while (i < N) {
        ctx.clear();
        while (i < N) {
            if (argument[i] == ',') break;
            ctx += toupper(argument[i]);
            ++i;
        }
        if (ctx == "CPG") {
            skip_cpg = false;
        } else if (ctx == "CHG") {
            skip_chg = false;
        } else if (ctx == "CHH") {
            skip_chh = false;
        } else {
            fprintf(stderr, "Illegal argument to option '%s': %s\n", option, argument);
            abort();
        }
        ++i;
    }
    return true;
}

bool ModOptions::parse(int argc, char* argv[]) {
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-v") == 0) {
            fprintf(stderr, "%s\n", HBN_PACKAGE_VERSION);
            exit (0);
        }
        if (strcmp(argv[i], "-h") == 0) {
            dump_usage(argc, argv);
            exit (0);
        }
    }

    int default_cpu_threads = get_physical_core_count();
    std::string default_model_dir = get_executable_dir();
    if (default_model_dir.size() && default_model_dir.back() != '/') {
        default_model_dir += '/';
    }
    default_model_dir += "models";

    int i = 2;
    const char* cstr = nullptr;
    while (i < argc) {
        if (argv[i][0] != '-') break;
        if (argv[i][0] == '-' && strlen(argv[i]) == 1) break;

        if (parse_string_arg_value(argc, argv, i, "-m", cstr)) {
            model_dir = cstr;
            continue;
        }
        if (parse_int_arg_value(argc, argv, i, "-l", min_read_size)) continue;
        if (parse_int_arg_value(argc, argv, i, "-s", sample_batch_size)) continue;
        if (parse_int_arg_value(argc, argv, i, "-b", read_batch_size)) continue;
        if (parse_bool_arg_value_x(argc, argv, i, "-k", keep_kinetics, true)) continue;
        if (parse_int_arg_value(argc, argv, i, "-t", num_threads)) continue;
        
        if (strcmp(argv[i], "-c") == 0) {
            if (i + 1 == argc) {
                fprintf(stderr, "ERROR: Argument to option '-c' is missing\n");
                abort();
            }
            s_parse_context_string(argv[i], argv[i+1], disable_cpg, disable_chg, disable_chh);
            i += 2;
            continue;
        }

        fprintf(stderr, "ERROR: unrecognised option %s", argv[i]);
        return false;
    }

    if (i >= argc) return false;
    bam_path = argv[i];
    ++i;

    if (i >= argc) return false;
    mod_bam_path = argv[i];
    ++i;

    if (i != argc) return false;

    if (model_dir.empty()) {
        if (!path_exists(default_model_dir)) {
            fprintf(stderr, "ERROR: the default model directory does not exist: %s\n", default_model_dir.c_str());
            fprintf(stderr, "Please specified the model directory via '-m' option\n");
            return false;
        }
        model_dir = default_model_dir;
    }

    if (num_threads == 0) {
        num_threads = default_cpu_threads;
    }

    return true;
}

void ModOptions::dump_usage(int argc, char* argv[]) {
    int default_cpu_threads = get_physical_core_count();
    std::string default_model_dir = get_executable_dir();
    if (default_model_dir.size() && default_model_dir.back() != '/') {
        default_model_dir += '/';
    }
    default_model_dir += "models";

    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "  %s %s [OPTIONS] BAM MOD-BAM\n", argv[0], argv[1]);

    fprintf(stderr, "\n\n");
    fprintf(stderr, "DESCRIPTION:\n");
    fprintf(stderr, "  Compute single molecular cytosine methylation states in BAM file reads\n");

    fprintf(stderr, "\n\n");
    fprintf(stderr, "VERSION:\n");
    fprintf(stderr, "  %s\n", HBN_PACKAGE_VERSION);

    fprintf(stderr, "\n\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS:\n");
    fprintf(stderr, "  -v\n");
    fprintf(stderr, "    Print version info and exit\n");
    fprintf(stderr, "  -h\n");
    fprintf(stderr, "    Print this help info and exit\n");
    fprintf(stderr, "  -m <Model directory>\n");
    fprintf(stderr, "    Path to the model directory\n");
    fprintf(stderr, "    Default: %s\n", default_model_dir.c_str());
    fprintf(stderr, "  -l <Integer>\n");
    fprintf(stderr, "    Minumum read length consider for 5mC calling\n");
    fprintf(stderr, "    Default = %d\n", kMinReadSize);
    fprintf(stderr, "  -s <Integer>\n");
    fprintf(stderr, "    Number of samples in one batch for 5mC calling\n");
    fprintf(stderr, "    Default = %d\n", kSampleBatchSize);
    fprintf(stderr, "  -b <Integer>\n");
    fprintf(stderr, "    Number of reads in one batch\n");
    fprintf(stderr, "    Default = %d\n", kReadBatchSize);
    fprintf(stderr, "  -k\n");
    fprintf(stderr, "    Keep kinetic values (fi, ri, fp, rp) in modified BAM output\n");
    fprintf(stderr, "  -c <string>\n");
    fprintf(stderr, "    5mC contexts to detect; comma separated\n");
    fprintf(stderr, "    Default = cpg,chg,chh\n");
    fprintf(stderr, "  -t <Integer>\n");
    fprintf(stderr, "    Number of CPU threads used\n");
    fprintf(stderr, "    Default = %d\n", default_cpu_threads);
}

void ModOptions::dump_parameters() 
{
        fprintf(stderr, "\n\n");
        fprintf(stderr, "######## Parameters:\n");
        fprintf(stderr, "  ## min-read-length: %d\n", min_read_size);
        fprintf(stderr, "  ## sample-batch-size: %d\n", sample_batch_size);
        fprintf(stderr, "  ## read-batch-size: %d\n", read_batch_size);
        fprintf(stderr, "  ## keep-kinetics: %s\n", keep_kinetics ? "yes" : "no");
        fprintf(stderr, "  ## skip-CpG: %s\n", disable_cpg ? "yes" : "no");
        fprintf(stderr, "  ## skip-CHG: %s\n", disable_chg ? "yes" : "no");
        fprintf(stderr, "  ## skip-CHH: %s\n", disable_chh ? "yes" : "no");
        fprintf(stderr, "  ## CPU-threads: %d\n", num_threads);
        fprintf(stderr, "  ## model-dir: %s\n", model_dir.c_str());
        fprintf(stderr, "  ## BAM-path: %s\n", bam_path.c_str());
        fprintf(stderr, "  ## MOD-bam-path: %s\n", mod_bam_path.c_str());
    }

} // ns_mods
