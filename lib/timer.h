#ifndef TIMER_H
#define TIMER_H

#include <cstdint>
#include <iostream>
#include <string>
#include <sys/time.h>

#define MICRO_SECS_PER_SEC 1'000'000ll

inline uint64_t get_micro_time()
{
    timeval t;
    gettimeofday(&t, 0);
    return t.tv_sec * MICRO_SECS_PER_SEC + t.tv_usec;
}

void measure_time(std::string name, auto F)
{
    std::cout << "Starting " << name << std::endl;
    long long cur = get_micro_time();
    F();
    cur = get_micro_time() - cur;
    std::cout << "Ending " << name << " after " << 1.0 * cur / MICRO_SECS_PER_SEC << " seconds" << std::endl;
}

#endif
