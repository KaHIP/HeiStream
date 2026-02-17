/****************************************************************************
 * stream_algorithm.h
 *****************************************************************************/

#ifndef CORE_STREAM_ALGORITHM_H_
#define CORE_STREAM_ALGORITHM_H_

#include <string>

#include "core/context/stream_context.h"
#include "partition/partition_config.h"

// Common interface for all streaming graph algorithms.
class IStreamAlgorithm {
   public:
    virtual ~IStreamAlgorithm() {}

    virtual void prepare_run(Config& config, const std::string& graph_filename,
                             StreamContext& ctx) {
        (void)config;
        (void)graph_filename;
        (void)ctx;
    }

    virtual void run(Config& config, const std::string& graph_filename,
                     StreamContext& ctx) = 0;

    virtual void finalize_run(Config& config, const std::string& graph_filename,
                              StreamContext& ctx) {
        (void)config;
        (void)graph_filename;
        (void)ctx;
    }
};

#endif /* CORE_STREAM_ALGORITHM_H_ */
