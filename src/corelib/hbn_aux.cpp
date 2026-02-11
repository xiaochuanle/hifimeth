#include "hbn_aux.hpp"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <locale.h>

#include <iostream>
#include <string>
#include <cctype>
#include <cmath>

#include <algorithm>
#include <iomanip>
#include <map>
#include <sstream>

//////////////////////////////////////////

/*
the IUPACna codes are single letters for nucleic acids and the value is the same 
    as the ASCII value of the recommended IUPAC letter. 
Value	Symbol	Name
65	A	Adenine
66	B	G or T or C
67	C	Cytosine
68	D	G or A or T
71	G	Guanine
72	H	A or C or T
75	K	G or T
77	M	A or C
78	N	A or G or C or T
82	R	G or A
83	S	G or C
84	T	Thymine
86	V	G or C or A
87	W	A or T
89	Y	T or C
*/

const Uint1 IUPACNA_TO_BLASTNA[128]={
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};

/////////////////////////////////////////

static const char* hbn_log_level_name[] = {
    [LOG_EMERG]     = "EMERG",
    [LOG_ALERT]     = "ALERT",
    [LOG_CRIT]      = "CRIT",
    [LOG_ERR]       = "ERROR",
    [LOG_WARNING]   = "WARNING",
    [LOG_NOTICE]    = "NOTICE",
    [LOG_INFO]      = "INFO",
    [LOG_DEBUG]     = "DEBUG",
};

void
what_is_time_now(char now[])
{
    time_t ltime;
    time(&ltime);
    ctime_r(&ltime, now);
    size_t n = strlen(now);
    now[n-1] = '\0';
}

static void 
_hbn_dump_message(HBN_LOG_PARAMS_GENERIC, const int level, const char* fmt, va_list args)
{
    char buf[HBN_MAX_PATH_LEN + 1];
    char time_str[256];
    what_is_time_now(time_str);
    vsnprintf(buf, sizeof(buf), fmt, args);
    fprintf(stderr, "[%s %s:%s:%d] <%s> %s\n", time_str, HBN_LOG_ARGS_GENERIC, hbn_log_level_name[level], buf);
}

void
hbn_dump_message(HBN_LOG_PARAMS_GENERIC, const int level, const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    _hbn_dump_message(HBN_LOG_ARGS_GENERIC, level, fmt, args);
    va_end(args);
}

static void
_hbn_dump_message_short(const int level, const char* fmt, va_list args)
{
    char buf[HBN_MAX_PATH_LEN + 1];
    char time_str[256];
    what_is_time_now(time_str);
    vsnprintf(buf, sizeof(buf), fmt, args);
    fprintf(stderr, "[%s] <%s> %s\n", time_str, hbn_log_level_name[level], buf);
}

void
hbn_dump_message_short(const int level, const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    _hbn_dump_message_short(level, fmt, args);
    va_end(args);    
}

////////////////////////////////////////////////////////////////////////////////////////////

int safe_gzread(HBN_LOG_PARAMS_GENERIC, gzFile stream, void* buf, unsigned int len)
{
    int ret = gzread(stream, buf, len);
    if (ret < 0) {
        int errnum = 0;
        const char* msg = gzerror(stream, &errnum);
        const char* why = (Z_ERRNO == errnum) ? strerror(errno) : msg;
        HBN_ERR_GENERIC("%s", why);
    }
    return ret;
}

int err_gzread(gzFile stream, void* buf, unsigned int len)
{
    return safe_gzread(HBN_LOG_ARGS_DEFAULT, stream, buf, len);
}

gzFile safe_gzopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode)
{
    gzFile stream;
    if (strcmp(path, "-") == 0) {
        stream = gzdopen(fileno(strstr(mode, "r") ? stdin : stdout), mode);
        /* according to zlib.h, this is the only reason gzdopen can fail */
        if (!stream) HBN_ERR_GENERIC("fail to open file '%s' with mode '%s': out of memory", path, mode);
        return stream;
    }
    if ((stream = gzopen(path, mode)) == 0) {
        const char* why = errno ? strerror(errno) : "out of memory";
        HBN_ERR_GENERIC("fail to open file '%s' with mode '%s': %s", path, mode, why);
    }
    return stream;
}

int safe_gzclose(HBN_LOG_PARAMS_GENERIC, gzFile stream)
{
    int ret = gzclose(stream);
    if (Z_OK != ret) {
        const char* why = (Z_ERRNO == ret) ? strerror(errno) : zError(ret);
        HBN_ERR_GENERIC("%s", why);
    }
    return ret;
}

FILE* safe_fopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode)
{
    FILE* stream = 0;
    if (strcmp(path, "-") == 0) return strstr(mode, "r") ? stdin : stdout;

    if ((stream = fopen(path, mode)) == 0) {
        const char* y = strerror(errno);
        HBN_ERR_GENERIC("fail to open file '%s' with mode '%s': %s", path, mode, y);
    }
    return stream;
}

size_t safe_fwrite(HBN_LOG_PARAMS_GENERIC, const void* buf, size_t size, size_t nmemb, FILE* stream)
{
    size_t ret = fwrite(buf, size, nmemb, stream);
    if (ret != nmemb) {
        HBN_ERR_GENERIC("%s", strerror(errno));
    }
    return ret;
}

size_t safe_fread(HBN_LOG_PARAMS_GENERIC, void* buf, size_t size, size_t nmemb, FILE* stream)
{
    size_t ret = fread(buf, size, nmemb, stream);
    if (ret != nmemb) {
        const char* y = ferror(stream) ? strerror(errno) : "Unexpected end of file";
        HBN_ERR_GENERIC("%s, expected: %zu, but only read %zu", y, nmemb, ret);
    }
    return ret;
}

int safe_fclose(HBN_LOG_PARAMS_GENERIC, FILE* stream)
{
    int ret = fclose(stream);
    if (ret != 0) {
        HBN_ERR_GENERIC("%s", strerror(errno));
    }
    return ret;
}

void hbn_exception(const char* expr, HBN_LOG_PARAMS_GENERIC, const char* fmt, ...)
{
    fprintf(stderr, "Assertion Failed At '%s:%s:%d'\n", file, func, line);
    fprintf(stderr, "\tExpression: '%s'\n", expr);
    if (!fmt) return;
    fprintf(stderr, "Context Information:\n");
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    abort();
}

off_t hbn_get_file_size(HBN_LOG_PARAMS_GENERIC, const char* path)
{
    struct stat sbuf;
    if (stat(path, &sbuf) == -1) {
        const char* y = strerror(errno);
        HBN_ERR_GENERIC("fail to stat file '%s': %s", path, y);
    }
    return sbuf.st_size;
}

/////////////////////////////////////

#include <cstddef>
#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__) || defined(__APPLE__)
#include <unistd.h>
#endif

size_t get_cpu_core_count() 
{
#if defined(_WIN32)
    // Windows 实现
    SYSTEM_INFO sys_info;
    GetNativeSystemInfo(&sys_info);  // 获取物理硬件信息[5](@ref)
    return sys_info.dwNumberOfProcessors;
#elif defined(__linux__) || defined(__APPLE__)
    // Linux/macOS 实现
    return sysconf(_SC_NPROCESSORS_ONLN);  // 获取当前可用核心数[1,3](@ref)
#else
    return 1;  // 不支持的系统返回默认值
#endif
}

#include <cstdint>

#if defined(_WIN32)
    #include <windows.h>
#elif defined(__linux__)
    #include <sys/sysinfo.h>
    #include <fstream>
    #include <string>
#elif defined(__APPLE__)
    #include <sys/sysctl.h>
    #include <mach/mach_host.h>
#endif

uint64_t get_system_memory_bytes() 
{
    uint64_t total_memory = 0;

    // Windows 实现
    #if defined(_WIN32)
        MEMORYSTATUSEX mem_info;
        mem_info.dwLength = sizeof(mem_info);
        if (GlobalMemoryStatusEx(&mem_info)) {
            total_memory = mem_info.ullTotalPhys; // 直接返回物理内存总量（字节）
        }

    // Linux 实现（优先使用 sysinfo，失败则解析 /proc/meminfo）
    #elif defined(__linux__)
        struct sysinfo info;
        if (sysinfo(&info) == 0) {
            total_memory = info.totalram * info.mem_unit; // 计算总内存字节数
        } else {
            // 备选方案：解析 /proc/meminfo
            std::ifstream meminfo("/proc/meminfo");
            std::string line;
            while (std::getline(meminfo, line)) {
                if (line.find("MemTotal:") != std::string::npos) {
                    uint64_t kb_size = 0;
                    sscanf(line.c_str(), "MemTotal: %lu kB", &kb_size);
                    total_memory = kb_size * 1024; // KB → 字节
                    break;
                }
            }
        }

    // macOS 实现
    #elif defined(__APPLE__)
        int mib[] = { CTL_HW, HW_MEMSIZE }; // 硬件内存大小参数
        size_t len = sizeof(total_memory);
        if (sysctl(mib, 2, &total_memory, &len, nullptr, 0) != 0) {
            total_memory = 0; // 失败时返回 0
        }
    #endif

    return total_memory;
}

void create_directory(const char* path)
{
    if ((access(path, F_OK) != 0)
        &&
        (mkdir(path, S_IRWXU) != 0)) {
        HBN_ERR("Failed to create directory %s: %s", path, strerror(errno));
    }
}

/////////////////////////////

#include <string>
#include <algorithm>
#include <cstdint>

std::string format_i64_with_commas(int64_t num) 
{
    // 处理零值
    if (num == 0) return "0";

    std::string result;
    bool isNegative = (num < 0);
    uint64_t absNum = isNegative ? static_cast<uint64_t>(-num) : static_cast<uint64_t>(num);

    // 反转处理：从低位向高位每3位插入逗号
    int count = 0;
    while (absNum > 0) {
        if (count > 0 && count % 3 == 0) {
            result.push_back(',');
        }
        result.push_back('0' + (absNum % 10));
        absNum /= 10;
        count++;
    }
    
    // 添加负号并反转字符串恢复正确顺序
    if (isNegative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

std::string format_u64_with_fixed_string(uint64_t num, size_t length) 
{
    std::string str = std::to_string(num);
    if (str.length() >= length) return str; // 无需补零
    
    // 补足前导零
    std::string zeros(length - str.length(), '0');
    return zeros + str;
}

std::pair<const char*, size_t> truncate_spaces(const char* str, const size_t length, const ETruncateEnd where)
{
    const std::pair<const char*, size_t> empty_str(nullptr, 0);
    if (length == 0) {
        return empty_str;
    }
    size_t beg = 0;
    if (where == eTrunc_Begin  ||  where == eTrunc_Both) {
        hbn_assert(beg < length);
        while ( isspace((unsigned char) str[beg]) ) {
            if (++beg == length) {
                return empty_str;
            }
        }
    }
    size_t end = length;
    if ( where == eTrunc_End  ||  where == eTrunc_Both ) {
        hbn_assert(beg < end);
        while (isspace((unsigned char) str[--end])) {
            if (beg == end) {
                return empty_str;
            }
        }
        hbn_assert(beg <= end  &&  !isspace((unsigned char) str[end]));
        ++end;
    }
    hbn_assert(beg < end  &&  end <= length);
    if ( beg | (end - length) ) { // if either beg != 0 or end != length
        return std::pair<const char*, size_t>(str + beg, end - beg);
    }
    else {
        return std::pair<const char*, size_t>(str, length);
    }
}

size_t datasize_to_bytes(const std::string& size_str) 
{
    if (size_str.empty()) return 0;

    // 单位到乘数的映射（不区分大小写）
    static const std::map<std::string, uint64_t> units = {
        {"B", 1ULL},
        {"KB", 1ULL << 10}, {"K", 1ULL << 10}, {"KIB", 1ULL << 10},
        {"MB", 1ULL << 20}, {"M", 1ULL << 20}, {"MIB", 1ULL << 20},
        {"GB", 1ULL << 30}, {"G", 1ULL << 30}, {"GIB", 1ULL << 30},
        {"TB", 1ULL << 40}, {"T", 1ULL << 40}, {"TIB", 1ULL << 40},
        {"PB", 1ULL << 50}, {"P", 1ULL << 50}, {"PIB", 1ULL << 50}
    };

    // 查找第一个非数字字符的位置
    size_t num_end = 0;
    while (num_end < size_str.size() && 
          (std::isdigit(size_str[num_end]) || size_str[num_end] == '.')) {
        num_end++;
    }

    // 提取数字部分
    double value;
    try {
        value = std::stod(size_str.substr(0, num_end));
    } catch (...) {
        return 0; // 数字转换失败
    }

    // 提取并处理单位部分
    std::string unit = size_str.substr(num_end);
    
    // 移除单位中的空格
    unit.erase(std::remove_if(unit.begin(), unit.end(), 
              [](unsigned char c) { return std::isspace(c); }), 
              unit.end());
    
    // 转换为大写
    std::transform(unit.begin(), unit.end(), unit.begin(),
                  [](unsigned char c) { return std::toupper(c); });

    // 查找匹配的单位
    auto it = units.find(unit);
    if (it != units.end()) {
        return static_cast<size_t>(value * it->second);
    }

    // 特殊处理无单位的情况
    if (unit.empty()) {
        return static_cast<size_t>(value);
    }

    return 0; // 未知单位
}

std::string bytes_to_datasize(uintmax_t bytes) 
{
    // 定义单位和对应的字节大小（基于1024）
    const char* const units[] = {"B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
    const uintmax_t max_unit = sizeof(units) / sizeof(units[0]) - 1;
    
    // 处理0字节情况
    if (bytes == 0) {
        return "0B";
    }

    // 计算合适的单位
    int unit_index = 0;
    double size = static_cast<double>(bytes);
    while (size >= 1024.0 && unit_index < max_unit) {
        size /= 1024.0;
        unit_index++;
    }

    // 格式化输出
    std::ostringstream oss;
    
    // 根据数值大小确定小数精度
    if (unit_index == 0 || size >= 100) {
        oss << static_cast<uintmax_t>(std::round(size));  // 整数部分
    } else if (size >= 10) {
        oss << std::fixed << std::setprecision(1) << size;  // 保留1位小数
    } else {
        oss << std::fixed << std::setprecision(2) << size;  // 保留2位小数
    }
    
    // 移除多余的小数0（如10.00 -> 10，2.50 -> 2.5）
    std::string result = oss.str();
    size_t pos = result.find('.');
    if (pos != std::string::npos) {
        // 移除末尾多余的0
        result.erase(result.find_last_not_of('0') + 1, std::string::npos);
        // 如果小数部分完全移除，移除小数点
        if (result.back() == '.') {
            result.pop_back();
        }
    }
    
    return result + units[unit_index];
}

std::vector<std::pair<const char*, size_t>> 
split_cstr(const char* s, size_t sl, char delimiter, bool skip_empty) 
{
    std::vector<std::pair<const char*, size_t>> result;
    if (s == nullptr || sl == 0) return result; // 处理空输入
    
    const char* start = s;       // 子串起始位置
    size_t count = 0;            // 当前子串长度计数器

    for (size_t i = 0; i < sl; ++i) {
        if (s[i] == delimiter) {
            // 遇到分隔符时保存当前子串
            if (!skip_empty || count > 0) {
                result.emplace_back(start, count);
            }
            start = &s[i + 1];  // 新子串起始位置（分隔符后）
            count = 0;
        } else {
            ++count;
        }
    }
    
    // 处理最后一个子串
    if (count > 0 || (!skip_empty && start < s + sl)) {
        result.emplace_back(start, count);
    }
    
    return result;
}
