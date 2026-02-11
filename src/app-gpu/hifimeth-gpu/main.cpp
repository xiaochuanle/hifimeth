#include "../../corelib/hbn_aux.hpp"

#include <map>
#include <string>
#include <thread>

typedef int main_fun_type(int argc, char* argv[]);

int single_molecular_methy_call(int argc, char* argv[]);

class HifiMethyCommandList
{
public:
    HifiMethyCommandList() {
        add_one_cmd("call", single_molecular_methy_call);
    }

    void add_one_cmd(const char* cmdname, main_fun_type* fun) {
        M_cmds[std::string(cmdname)] = fun;
    }

    void dump_cmds(FILE* out = stderr) {
        fprintf(stderr, "\n");
        fprintf(stderr, "COMMANDS:\n");
        for (auto& cmd : M_cmds) {
            fprintf(stderr, "  %s\n", cmd.first.c_str());
        }
    }

    int run_cmd(int argc, char* argv[]) {
        if (argc < 2) return 1;
        auto pos = M_cmds.find(std::string(argv[1]));
        if (pos == M_cmds.end()) {
            fprintf(stderr, "ERROR: Unrecognised command '%s'\n", argv[1]);
            return 1;
        }
        return (*pos->second)(argc, argv);
    }

private:
    std::map<std::string, main_fun_type*>    M_cmds;
};

int main(int argc, char* argv[])
{
    HifiMethyCommandList cmds;
    int r = cmds.run_cmd(argc, argv);
    if (r) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s COMMANDS\n", argv[0]);
        cmds.dump_cmds();
        return 1;
    }

    return 0;
}
