#include "../../corelib/hbn_aux.hpp"

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <thread>

typedef int main_fun_type(int argc, char* argv[]);

int mod_main(int argc, char* argv[]);
int pileup_main(int argc, char* argv[]);
int pileup_correlation_main(int argc, char* argv[]);
int cov_to_bed_main(int argc, char* argv[]);
int subsample_bam_main(int argc, char* argv[]);
int eval_main(int argc, char* argv[]);

struct HbnCommand
{
    std::string app;
    std::string description;
    main_fun_type* app_main;

    HbnCommand(const char* _app, const char* _description, main_fun_type* _app_main)
        : app(_app), description(_description), app_main(_app_main) {}

    void dump_usage() const {
        //std::cerr << std::left << std::setw(10) << "  " << app << description << std::endl;
        fprintf(stderr, "  ");
        fprintf(stderr, "%-15s", app.c_str());
        fprintf(stderr, "%s\n", description.c_str());
    }    
};

static const HbnCommand
AppCall("call", "Detect per-read 5mC methylation from PacBio HiFi reads", mod_main);
static const HbnCommand
AppPileup("pileup", "Compute methylation frequency at genomic 5mC site from aligned BAM", pileup_main);
static const HbnCommand
AppCorrelation("corr", "Calculate Pearson correlation coefficient between two pileup files.", pileup_correlation_main);
static const HbnCommand
AppCovToBed("cov2bed", "Convert 1-based Bismark .cov file to 0-based BED format", cov_to_bed_main);
static const HbnCommand
AppSubsample("sample", "Randomly subsample reads from a BAM file", subsample_bam_main);
static const HbnCommand
AppEval("eval", "Extract read 5mC sites for benchmark", eval_main);

class HifiMethyCommandList
{
public:
    HifiMethyCommandList() {
        M_apps.push_back(&AppCall);
        M_apps.push_back(&AppPileup);
        M_apps.push_back(&AppCorrelation);
        M_apps.push_back(&AppCovToBed);
        M_apps.push_back(&AppSubsample);
        M_apps.push_back(&AppEval);
    }

    void dump_cmds(int argc, char* argv[]) {
        fprintf(stderr, "\n");
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "  %s <commands> [options]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "COMMANDS:\n");
        for (auto app : M_apps) app->dump_usage();
        fprintf(stderr, "\n");
    }

    int run_cmd(int argc, char* argv[]) {
        if (argc >= 2) {
            for (auto app : M_apps) {
                if (app->app == argv[1]) {
                    return app->app_main(argc, argv);
                }
            }
        }
        if (argc >= 2) {
            fprintf(stderr, "\n");
            fprintf(stderr, "ERROR: Unrecognised command '%s'\n", argv[1]);
        }
        dump_cmds(argc, argv);
        return EXIT_FAILURE;   
    }

private:
    std::vector<const HbnCommand*>     M_apps;
};

int main(int argc, char* argv[])
{
    return HifiMethyCommandList().run_cmd(argc, argv);
}
