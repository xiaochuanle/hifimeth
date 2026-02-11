#ifndef __PROGRAM_INFO_HPP
#define __PROGRAM_INFO_HPP

#include <corelib/hbn_aux.hpp>

#include <string>

#include <cstdio>

#define HBN_PACKAGE                         1
#define HBN_PACKAGE_NAME                    "hifimeth"
#define HBN_PACKAGE_VERSION_MAJOR           1
#define HBN_PACKAGE_VERSION_MINOR           1
#define HBN_PACKAGE_VERSION_PATCH           0
#define HBN_PACKAGE_CONFIG                  ""

#define HBN_PACKAGE_VERSION_STRINGIFY(x)    #x
#define HBN_PACKAGE_VERSION_COMPOSE_STR(a, b, c)    \
    HBN_PACKAGE_VERSION_STRINGIFY(a) "."            \
    HBN_PACKAGE_VERSION_STRINGIFY(b) "."            \
    HBN_PACKAGE_VERSION_STRINGIFY(c)

#define HBN_PACKAGE_VERSION             \
    HBN_PACKAGE_VERSION_COMPOSE_STR     \
    (                                   \
        HBN_PACKAGE_VERSION_MAJOR,      \
        HBN_PACKAGE_VERSION_MINOR,      \
        HBN_PACKAGE_VERSION_PATCH       \
    )

typedef void (*program_desc_func)(FILE* out);

struct HbnProgramInfo
{
public:
    HbnProgramInfo(const char* program, program_desc_func desc) {
        M_program = program;
        if (desc) (*desc)(stderr);
        gettimeofday(&M_begin, NULL);
    }

    ~HbnProgramInfo();
    
private:
    std::string     M_program;
    struct timeval  M_begin;
    struct timeval  M_end;
};

#endif // __PROGRAM_INFO_HPP