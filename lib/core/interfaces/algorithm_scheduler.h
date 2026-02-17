/****************************************************************************
 * algorithm_scheduler.h
 *****************************************************************************/

#ifndef CORE_ALGORITHM_SCHEDULER_H_
#define CORE_ALGORITHM_SCHEDULER_H_

#include <string>

#include "core/context/stream_context.h"
#include "partition/partition_config.h"

namespace core::interfaces {

/**
 * Executes the algorithm-specific stream pipeline.
 */
class IAlgorithmScheduler {
   public:
    virtual ~IAlgorithmScheduler() {}

    /**
     * Run the algorithm over the configured input stream.
     *
     * @param config Mutable run configuration.
     * @param graph_filename Input graph path.
     * @param ctx Shared runtime context (timing, pass state, memory stats).
     */
    virtual void execute(Config& config, const std::string& graph_filename,
                         StreamContext& ctx) const = 0;
};

} // namespace core::interfaces


#endif /* CORE_ALGORITHM_SCHEDULER_H_ */
