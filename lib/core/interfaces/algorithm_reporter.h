/****************************************************************************
 * algorithm_reporter.h
 *****************************************************************************/

#ifndef CORE_ALGORITHM_REPORTER_H_
#define CORE_ALGORITHM_REPORTER_H_

#include <string>

#include "core/context/stream_context.h"
#include "partition/partition_config.h"

class FlatBufferWriter;


namespace core::interfaces {

/**
 * Emits run summaries and final logging artifacts.
 */
class IAlgorithmReporter {
   public:
    virtual ~IAlgorithmReporter() {}

    /**
     * Print human-readable run configuration prior to execution.
     *
     * @param config Immutable configuration snapshot.
     * @param graph_filename Input graph path.
     */
    virtual void print_run_configuration(const Config& config,
                                         const std::string& graph_filename) const {
        (void)config;
        (void)graph_filename;
    }

    /**
     * Emit end-of-run artifacts and summary output.
     *
     * @param config Mutable run configuration.
     * @param graph_filename Input graph path.
     * @param ctx Shared runtime context.
     * @param fb_writer Flatbuffer writer with run/evaluation fields.
     */
    virtual void report(Config& config, const std::string& graph_filename, StreamContext& ctx,
                        FlatBufferWriter& fb_writer) const = 0;
};

} // namespace core::interfaces


#endif /* CORE_ALGORITHM_REPORTER_H_ */
