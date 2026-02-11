#ifndef __GET_CORE_COUNT_HPP
#define __GET_CORE_COUNT_HPP

#include <string>

unsigned int get_physical_core_count();

std::string get_executable_dir();

bool path_exists(const std::string& name);

#endif // __GET_CORE_COUNT_HPP