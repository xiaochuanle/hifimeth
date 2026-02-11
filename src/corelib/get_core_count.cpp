#include "get_core_count.hpp"

#include <iostream>
#include <thread>

// 根据操作系统引入必要的头文件
#if defined(_WIN32)
    #include <windows.h>
    #include <malloc.h>
#elif defined(__APPLE__)
    #include <sys/sysctl.h>
    #include <sys/types.h>
#elif defined(__linux__)
    #include <fstream>
    #include <string>
    #include <set>
    #include <sstream>
#endif

// 获取物理核心数的函数
unsigned int get_physical_core_count() {
#if defined(_WIN32)
    // === Windows 实现 ===
    DWORD length = 0;
    GetLogicalProcessorInformation(nullptr, &length);
    
    if (GetLastError() != ERROR_INSUFFICIENT_BUFFER) {
        return std::thread::hardware_concurrency(); // 出错回退
    }

    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION buffer = 
        (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION)malloc(length);
    
    if (GetLogicalProcessorInformation(buffer, &length) == FALSE) {
        free(buffer);
        return std::thread::hardware_concurrency(); // 出错回退
    }

    unsigned int physical_cores = 0;
    DWORD ptr = 0;
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION current = nullptr;
    DWORD offset = 0;

    // 遍历处理器信息结构体
    while (offset + sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION) <= length) {
        current = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION)((char*)buffer + offset);
        offset += sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);

        // RelationProcessorCore 代表物理核心
        if (current->Relationship == RelationProcessorCore) {
            physical_cores++;
        }
    }
    
    free(buffer);
    return (physical_cores > 0) ? physical_cores : std::thread::hardware_concurrency();

#elif defined(__APPLE__)
    // === macOS 实现 ===
    int count;
    size_t size = sizeof(count);
    // "hw.physicalcpu" 获取物理核心数
    if (sysctlbyname("hw.physicalcpu", &count, &size, nullptr, 0) == 0) {
        return static_cast<unsigned int>(count);
    }
    return std::thread::hardware_concurrency(); // 出错回退

#elif defined(__linux__)
    // === Linux 实现 ===
    // Linux 下最可靠的方法是解析 /proc/cpuinfo
    // 我们需要统计不重复的 (physical id, core id) 对
    std::ifstream cpuinfo("/proc/cpuinfo");
    if (!cpuinfo.is_open()) {
        return std::thread::hardware_concurrency();
    }

    std::set<std::pair<int, int>> unique_cores;
    std::string line;
    int current_physical_id = -1;
    int current_core_id = -1;

    while (std::getline(cpuinfo, line)) {
        if (line.find("physical id") == 0) {
            // 解析 physical id
            size_t colon = line.find(':');
            if (colon != std::string::npos) {
                current_physical_id = std::stoi(line.substr(colon + 1));
            }
        } else if (line.find("core id") == 0) {
            // 解析 core id
            size_t colon = line.find(':');
            if (colon != std::string::npos) {
                current_core_id = std::stoi(line.substr(colon + 1));
            }
        }
        
        // 每次遇到一个新的 block 结束（通常是空行或新 processor），
        // 或者是解析完一对 id 后，尝试插入 set
        // 但 /proc/cpuinfo 格式是每个逻辑核一段。
        // 只要我们拿到一对有效值，存入 set 即可去重。
        if (current_physical_id != -1 && current_core_id != -1) {
            unique_cores.insert({current_physical_id, current_core_id});
            // 这里的逻辑简化处理：不用每次清空，因为 cpuinfo 顺序通常固定
            // 为了严谨，应该在遇到 "processor" 字段时重置，
            // 但 set 自动去重机制允许我们多次插入同样的值。
        }
    }
    
    // 如果解析失败（比如容器环境屏蔽了 cpuinfo 细节），回退
    if (unique_cores.empty()) {
        return std::thread::hardware_concurrency();
    }
    
    return static_cast<unsigned int>(unique_cores.size());

#else
    // === 其他平台 ===
    // 无法区分，直接返回逻辑线程数
    return std::thread::hardware_concurrency();
#endif
}

#if 0
int main() {
    unsigned int logical = std::thread::hardware_concurrency();
    unsigned int physical = get_physical_core_count();

    std::cout << "逻辑核心数 (Logical): " << logical << std::endl;
    std::cout << "物理核心数 (Physical): " << physical << std::endl;

    if (logical > physical) {
        std::cout << "检测到超线程 (Hyper-Threading) 或 SMT 已开启。" << std::endl;
    } else {
        std::cout << "未检测到超线程，或系统无法区分。" << std::endl;
    }

    return 0;
}
#endif

/////////////////////////////////////////


#include <iostream>
#include <string>
#include <vector>
#include <algorithm> // for find_last_of

// 根据不同平台引入头文件
#if defined(_WIN32)
    #include <windows.h>
#elif defined(__linux__)
    #include <unistd.h>
    #include <limits.h>
#elif defined(__APPLE__)
    #include <mach-o/dyld.h>
    #include <limits.h>
#endif

// 获取可执行文件所在的目录
std::string get_executable_dir() {
    std::string full_path;
    
    // === Windows 实现 ===
#if defined(_WIN32)
    char buffer[MAX_PATH];
    // GetModuleFileName 获取当前进程的完整路径
    // 这里的 NULL 表示获取当前进程模块
    DWORD length = GetModuleFileNameA(NULL, buffer, MAX_PATH);
    if (length == 0) return ""; // 获取失败
    full_path = std::string(buffer, length);

    // Windows 路径分隔符是 '\'
    size_t last_slash_idx = full_path.find_last_of("\\");
    if (std::string::npos != last_slash_idx) {
        // 截取到最后一个分隔符之前
        return full_path.substr(0, last_slash_idx);
    }

    // === Linux 实现 ===
#elif defined(__linux__)
    char buffer[PATH_MAX];
    // /proc/self/exe 是一个指向当前可执行文件的符号链接
    ssize_t length = readlink("/proc/self/exe", buffer, PATH_MAX);
    if (length <= 0) return ""; // 获取失败
    full_path = std::string(buffer, length);

    // Linux 路径分隔符是 '/'
    size_t last_slash_idx = full_path.find_last_of("/");
    if (std::string::npos != last_slash_idx) {
        return full_path.substr(0, last_slash_idx);
    }

    // === macOS 实现 ===
#elif defined(__APPLE__)
    char buffer[PATH_MAX];
    uint32_t size = sizeof(buffer);
    // _NSGetExecutablePath 获取完整路径
    if (_NSGetExecutablePath(buffer, &size) != 0) {
        return ""; // 缓冲区不够或出错
    }
    full_path = std::string(buffer);

    // macOS 路径分隔符是 '/'
    size_t last_slash_idx = full_path.find_last_of("/");
    if (std::string::npos != last_slash_idx) {
        return full_path.substr(0, last_slash_idx);
    }

#else
    // 不支持的平台
    return "";
#endif

    return ""; // 如果没找到分隔符（极其罕见），返回空或原始路径
}

#if 0
int main() {
    std::string dir = get_executable_dir();
    std::cout << "Executable Directory: " << dir << std::endl;
    
    // 测试：构建一个相对路径的文件
    std::string config_path = dir + "/config.ini"; 
    // 注意：Windows下如果需要严谨，应该判断 dir 结尾是否已有斜杠，
    // 或者在上面代码中保留斜杠。这里简单演示拼接。
    
    return 0;
}
#endif

/////////////////////////////////

#include <iostream>
#include <string>

// 根据平台引入头文件
#ifdef _WIN32
    #include <io.h>      // _access
    #include <windows.h> // GetFileAttributes
#else
    #include <unistd.h>  // access
    #include <sys/stat.h> // stat
#endif

bool path_exists(const std::string& name) {
#ifdef _WIN32
    // Windows 方式 1: 使用 _access (类似于 Linux 的 access)
    // 0 表示检查是否存在
    // return _access(name.c_str(), 0) != -1;

    // Windows 方式 2: 使用 Windows API (更推荐，能够区分文件和目录属性)
    DWORD attr = GetFileAttributesA(name.c_str());
    return (attr != INVALID_FILE_ATTRIBUTES); 
#else
    // Linux/macOS 方式: 使用 stat
    struct stat buffer;   
    return (stat(name.c_str(), &buffer) == 0); 
#endif
}

#if 0
int main() {
    std::string path = "test_folder"; 
    
    if (exists(path)) {
        std::cout << "存在!" << std::endl;
    } else {
        std::cout << "不存在!" << std::endl;
    }
    return 0;
}
#endif