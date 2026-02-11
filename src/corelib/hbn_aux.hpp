#ifndef __HBN_AUX_HPP
#define __HBN_AUX_HPP

#include <zlib/zlib.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdio.h>
#include <syslog.h>
#include <time.h>
#include <sys/time.h>

#include <random>
#include <string>
#include <vector>

///////////////////////////////////
// type alias

typedef int8_t      i8;
typedef uint8_t     u8;
typedef int16_t     i16;
typedef uint16_t    u16;
typedef int32_t     i32;
typedef uint32_t    u32;
typedef int64_t     i64;
typedef uint64_t    u64;

typedef i8          Int1;
typedef u8          Uint1;
typedef i16         Int2;
typedef u16         Uint2;
typedef i32         Int4;
typedef u32         Uint4;
typedef i64         Int8;
typedef u64         Uint8;

typedef i64         idx;
typedef u64         uidx;
#define PRIdx       PRId64
#define PRUidx      PRIu64

#define I8_MAX      INT8_MAX
#define U8_MAX      UINT8_MAX
#define I16_MIN     INT16_MIN
#define I16_MAX     INT16_MAX
#define U16_MAX     UINT16_MAX
#define I32_MIN     INT32_MIN
#define I32_MAX     INT32_MAX
#define I64_MIN     INT64_MIN
#define I64_MAX     INT64_MAX
#define U64_MAX     UINT64_MAX
#define IDX_MAX     INT64_MAX
#define UIDX_MAX    UINT64_MAX

#define U8_ONE      ((u8)1)
#define U32_ONE		((u32)1)
#define U64_ONE		((u64)1)
#define UIDX_ONE    U64_ONE

#define FWD         (0)
#define REV         (1)
#define F_R         (2)
#define REVERSED_STRAND(__s)   (1-(__s))

#define GAP_CHAR            ('-')
#define GAP_CODE            (4)
#define DECODE_RESIDUE(__r) ("ACGT-"[(u64)(__r)])

extern const Uint1 IUPACNA_TO_BLASTNA[];

#define HBN_MAX_PATH_LEN    2000

/////////////////////////////
// log

#define HBN_LOG_ARGS_DEFAULT    __FILE__, __FUNCTION__, __LINE__
#define HBN_LOG_ARGS_GENERIC    file, func, line
#define HBN_LOG_PARAMS_GENERIC  const char* file, const char* func, const int line

void
what_is_time_now(char now[]);

void
hbn_dump_message(HBN_LOG_PARAMS_GENERIC, const int level, const char* fmt, ...);

void
hbn_dump_message_short(const int level, const char* fmt, ...);

#if 0
#define HBN_LOG(fmt, args...) \
    hbn_dump_message(HBN_LOG_ARGS_DEFAULT, LOG_INFO, fmt, ##args)
#else
#define HBN_LOG(fmt, args...) \
    hbn_dump_message_short(LOG_INFO, fmt, ##args)
#endif

#define HBN_WARN(fmt, args...) \
    hbn_dump_message(HBN_LOG_ARGS_DEFAULT, LOG_WARNING, fmt, ##args)

#define HBN_ERR(fmt, args...) \
    do { hbn_dump_message(HBN_LOG_ARGS_DEFAULT, LOG_ERR, fmt, ##args); abort(); } while(0)

#define HBN_ERR_GENERIC(fmt, args...) \
    do { hbn_dump_message(HBN_LOG_ARGS_GENERIC, LOG_ERR, fmt, ##args); abort(); } while(0)

/// gzFile

int safe_gzread(HBN_LOG_PARAMS_GENERIC, gzFile stream, void* buf, unsigned int len);
int err_gzread(gzFile stream, void* buf, unsigned int len);
gzFile safe_gzopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode);
int safe_gzclose(HBN_LOG_PARAMS_GENERIC, gzFile stream);

#define hbn_gzopen(stream, path, mode)      stream = safe_gzopen(HBN_LOG_ARGS_DEFAULT, path, mode)
#define hbn_dgzopen(stream, path, mode)     gzFile stream; hbn_gzopen(stream, path, mode)
#define hbn_gzclose(stream)                 safe_gzclose(HBN_LOG_ARGS_DEFAULT, stream)
#define hbn_gzread(stream, buf, len)        safe_gzread(HBN_LOG_ARGS_DEFAULT, stream, buf, len)

/// file

#define HBN_SCANF(input_func, stream, nread, fmt, ...) \
do { \
    int __nread = input_func(stream, fmt, __VA_ARGS__); \
    if (nread != __nread) { \
        HBN_ERR("scanf error: should read %d items, but only %d have been read", nread, __nread); \
    } \
} while(0)

FILE* safe_fopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode);
size_t safe_fwrite(HBN_LOG_PARAMS_GENERIC, const void* buf, size_t size, size_t nmemb, FILE* stream);
size_t safe_fread(HBN_LOG_PARAMS_GENERIC, void* buf, size_t size, size_t nmemb, FILE* stream);
int safe_fclose(HBN_LOG_PARAMS_GENERIC, FILE* stream);
off_t hbn_get_file_size(HBN_LOG_PARAMS_GENERIC, const char* path);

#define hbn_fopen(stream, path, mode)           stream = safe_fopen(HBN_LOG_ARGS_DEFAULT, path, mode)
#define hbn_dfopen(stream, path, mode)          FILE* stream; hbn_fopen(stream, path, mode)
#define hbn_fwrite(buf, size, nmemb, stream)    safe_fwrite(HBN_LOG_ARGS_DEFAULT, buf, size, nmemb, stream)
#define hbn_fread(buf, size, nmemb, stream)     safe_fread(HBN_LOG_ARGS_DEFAULT, buf, size, nmemb, stream)
#define hbn_fclose(stream)                      safe_fclose(HBN_LOG_ARGS_DEFAULT, stream)
#define hbn_file_size(path)                     hbn_get_file_size(HBN_LOG_ARGS_DEFAULT, path)

/// assert

void
hbn_exception(const char* expr, HBN_LOG_PARAMS_GENERIC, const char* fmt, ...);

#define __hbn_assert(expr, ...) \
    do { \
        if (!(expr)) { \
            hbn_exception(#expr, __VA_ARGS__, NULL); \
            abort(); \
        } \
    } while(0)

#define hbn_assert(expr, args...) __hbn_assert(expr, HBN_LOG_ARGS_DEFAULT, ##args)

//////////////////////////////////////

size_t get_cpu_core_count();
uint64_t get_system_memory_bytes();
void create_directory(const char* path);

///////////////////////

std::string format_i64_with_commas(int64_t num);
std::string format_u64_with_fixed_string(uint64_t num, size_t length);
std::vector<std::pair<const char*, size_t>> 
split_cstr(const char* s, size_t sl, char delimiter, bool skip_empty = true);

/// Which end to truncate a string.
enum ETruncateEnd {
    eTrunc_Begin,  ///< Truncate leading whitespace only
    eTrunc_End,    ///< Truncate trailing whitespace only
    eTrunc_Both    ///< Truncate whitespace at both begin and end of string
};
std::pair<const char*, size_t> truncate_spaces(const char* str, const size_t length, const ETruncateEnd where);

std::string bytes_to_datasize(uintmax_t bytes);
size_t datasize_to_bytes(const std::string& size_str);

/////////////////////

class RandomFloatNumberGenerator
{
public:
    RandomFloatNumberGenerator(double min_v = 0.0, double max_v = 1.0)
    {
        M_rd = new std::random_device;
        M_gen = new std::mt19937((*M_rd)());
        M_dist = new std::uniform_real_distribution<double>(min_v, max_v);
    }

    ~RandomFloatNumberGenerator()
    {
        delete M_rd;
        delete M_gen; 
        delete M_dist;
    }

    double operator() () 
    {
        return (*M_dist)(*M_gen);
    }

    std::mt19937& gen() {
        return *M_gen;
    }

private:
    std::random_device* M_rd;
    std::mt19937*       M_gen;
    std::uniform_real_distribution<double>*  M_dist;
};

#endif // __HBN_AUX_HPP
