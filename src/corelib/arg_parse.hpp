#ifndef __ARG_PARSE_HPP
#define __ARG_PARSE_HPP

#include <cstring>

#include "hbn_aux.hpp"

static inline bool parse_bool_arg_value_x(int argc, char* argv[], int& i, const char* arg_name, bool& s, bool v)
{
    if (strcmp(argv[i], arg_name)) return false;

    s = v;
    i += 1;

    return true;
}

static inline bool parse_bool_arg_value(int argc, char* argv[], int& i, const char* arg_name, bool& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    s = true;
    i += 1;

    return true;
}

static inline bool parse_int_arg_value(int argc, char* argv[], int& i, const char* arg_name, int& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    if (i == argc) {
        fprintf(stderr, "ERROR: value to argument '%s' is missing\n", argv[i]);
        exit(1);
    }
    s = atoll(argv[i+1]);
    i += 2;

    return true;
}

static inline bool parse_real_arg_value(int argc, char* argv[], int& i, const char* arg_name, double& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    if (i == argc) {
        fprintf(stderr, "ERROR: value to argument '%s' is missing\n", argv[i]);
        exit(1);
    }
    s = atof(argv[i+1]);
    i += 2;

    return true;
}

static inline bool parse_data_size_arg_value(int argc, char* argv[], int& i, const char* arg_name, size_t& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    if (i == argc) {
        fprintf(stderr, "ERROR: value to argument '%s' is missing\n", argv[i]);
        exit(1);
    }
    s = datasize_to_bytes(argv[i+1]);
    i += 2;

    return true;
}

static inline bool parse_string_arg_value(int argc, char* argv[], int& i, const char* arg_name, const char*& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    if (i == argc) {
        fprintf(stderr, "ERROR: value to argument '%s' is missing\n", argv[i]);
        exit(1);
    }
    s = argv[i+1];
    i += 2;

    return true;
}

#endif // __ARG_PARSE_HPP