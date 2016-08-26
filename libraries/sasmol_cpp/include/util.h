#pragma once

#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <iostream> 
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <netdb.h>
#include <common.h>

#include <vector> 
#include <regex>

namespace util
{
    void Error(const std::string &);

    std::string time_and_user(); ///< (Brief description)

    void print_time() ; ///< (Brief description)

    void compile_time(const char[11], const char[11], const char[11], const char[11]) ; ///< (Brief description)

    void print(std::vector<std::string> &s); ///< (Brief description)

    void print(std::vector<int> &v); ///< (Brief description)

    void print(std::vector<float> &v); ///< (Brief description)

    std::vector<std::string> strip_white_space(std::vector<std::string> &s) ; ///< (Brief description)

    void print_run_details() ; ///< (Brief description)

    void pp(const std::string &s) ; ///< (Brief description)

    bool has_only_spaces(const std::string &str) ; ///< (Brief description)

    std::vector<int> unique_vector_int(const std::vector<int> & in); ///< (Brief description)

    std::vector<std::string> unique_vector_string(const std::vector<std::string> & in); ///< (Brief description)
        
    unsigned long int random_seed(); ///< generate a random seed

    float get_random_float(const float a, const float b); ///< generate a random number between [a,b)

    long long diff_ms(struct timeval t1, struct timeval t2); ///< (Brief description)

    void profiling(struct timeval & t1, struct timeval & t2, const std::string & moduleName); ///< (Brief description)
}
