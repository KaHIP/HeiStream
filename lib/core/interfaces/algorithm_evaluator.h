/****************************************************************************
 * algorithm_evaluator.h
 *****************************************************************************/

#ifndef CORE_ALGORITHM_EVALUATOR_H_
#define CORE_ALGORITHM_EVALUATOR_H_

#include <string>

#include "core/context/stream_context.h"
#include "partition/partition_config.h"

class FlatBufferWriter;


namespace core::interfaces {

/**
 * Evaluates algorithm outputs and updates reporting payloads if enabled.
 */
class IAlgorithmEvaluator {
   public:
    virtual ~IAlgorithmEvaluator() {}

    /**
     * Evaluate algorithm output and update quality fields in the flatbuffer writer.
     *
     * @param config Mutable run configuration.
     * @param graph_filename Input graph path (used for quality evaluation).
     * @param ctx Shared runtime context containing algorithm state.
     * @param fb_writer Flatbuffer payload writer to update with metrics.
     */
    virtual void evaluate(Config& config, const std::string& graph_filename, StreamContext& ctx,
                          FlatBufferWriter& fb_writer) const = 0;
};

} // namespace core::interfaces


#endif /* CORE_ALGORITHM_EVALUATOR_H_ */
