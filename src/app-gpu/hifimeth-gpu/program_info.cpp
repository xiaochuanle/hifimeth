#include "program_info.hpp"

using namespace std;

extern "C"
size_t getMemorySizeBytes();
extern "C" size_t getPeakRSS();

static double hbn_time_diff(const struct timeval* begin, const struct timeval* end)
{
    double d = end->tv_sec - begin->tv_sec;
    d += 1.0 * (end->tv_usec - begin->tv_usec) / 1e6;
    return d;
}

HbnProgramInfo::~HbnProgramInfo() {
    gettimeofday(&M_end, NULL);
    double dur = hbn_time_diff(&M_begin, &M_end);
    size_t peak_ram = getPeakRSS();
    string size = bytes_to_datasize(peak_ram);
    fprintf(stderr, "\n\n");
    fprintf(stderr, "%s Wallclock time: %.2f seconds.\n", M_program.c_str(), dur);
    fprintf(stderr, "%s Peak RAM usage: %s\n", M_program.c_str(), size.c_str());
}