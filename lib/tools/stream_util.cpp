/******************************************************************************
 * stream_util.cpp
 *****************************************************************************/

#include "tools/stream_util.h"

#include <sys/resource.h>

#include <iostream>

long getMaxRSS() {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // The maximum resident set size is in kilobytes.
        return usage.ru_maxrss;
    }
    std::cerr << "Error getting resource usage information." << '\n';
    return -1;
}

std::string extractBaseFilename(const std::string& fullPath) {
    size_t lastSlash = fullPath.find_last_of('/');
    size_t lastDot = fullPath.find_last_of('.');

    if (lastSlash != std::string::npos) {
        return fullPath.substr(lastSlash + 1, lastDot - lastSlash - 1);
    }
    return fullPath.substr(0, lastDot);
}
