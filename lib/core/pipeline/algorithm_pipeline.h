/****************************************************************************
 * algorithm_pipeline.h
 *****************************************************************************/

#ifndef CORE_ALGORITHM_PIPELINE_H_
#define CORE_ALGORITHM_PIPELINE_H_

#include <string>

#include "core/context/stream_context.h"
#include "partition/partition_config.h"

namespace core {
namespace interfaces {
class IAlgorithmScheduler;
class IAlgorithmEvaluator;
class IAlgorithmReporter;
}  // namespace interfaces
namespace pipeline {

// Shared lifecycle contract for all registered streaming algorithms.
void execute_algorithm_pipeline(Config& config, const std::string& graph_filename,
                                StreamContext& ctx, const interfaces::IAlgorithmScheduler& scheduler,
                                const interfaces::IAlgorithmEvaluator& evaluator,
                                const interfaces::IAlgorithmReporter& reporter);

}  // namespace pipeline
}  // namespace core

#endif /* CORE_ALGORITHM_PIPELINE_H_ */
