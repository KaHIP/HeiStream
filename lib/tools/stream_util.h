/******************************************************************************
 * stream_util.h
 *****************************************************************************/

#ifndef STREAM_UTIL_H_
#define STREAM_UTIL_H_

#include <string>

// Return process peak RSS in kilobytes.
long getMaxRSS();
// Extract filename stem from full path.
std::string extractBaseFilename(const std::string& fullPath);

#endif /* STREAM_UTIL_H_ */
