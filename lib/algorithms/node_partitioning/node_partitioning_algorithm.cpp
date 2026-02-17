/****************************************************************************
 * node_partitioning_algorithm.cpp
 *****************************************************************************/

#include "algorithms/node_partitioning/node_partitioning_algorithm.h"

#include "core/pipeline/algorithm_pipeline.h"

namespace node_partitioning {

NodePartitioningAlgorithm::NodePartitioningAlgorithm(
    std::unique_ptr<core::interfaces::IAlgorithmSetup> setup,
    std::unique_ptr<core::interfaces::IAlgorithmScheduler> scheduler,
    std::unique_ptr<core::interfaces::IAlgorithmEvaluator> evaluator,
    std::unique_ptr<core::interfaces::IAlgorithmReporter> reporter)
    : setup_(std::move(setup)),
      scheduler_(std::move(scheduler)),
      evaluator_(std::move(evaluator)),
      reporter_(std::move(reporter)) {}

void NodePartitioningAlgorithm::prepare_run(Config& config, const std::string& graph_filename,
                                            StreamContext& ctx) {
    (void)ctx;
    reporter_->print_run_configuration(config, graph_filename);
    setup_->prepare_run(config, graph_filename);
}

void NodePartitioningAlgorithm::run(Config& config, const std::string& graph_filename,
                                    StreamContext& ctx) {
    core::pipeline::execute_algorithm_pipeline(config, graph_filename, ctx, *scheduler_, *evaluator_,
                                               *reporter_);
}

void NodePartitioningAlgorithm::finalize_run(Config& config, const std::string& graph_filename,
                                             StreamContext& ctx) {
    (void)config;
    (void)graph_filename;
    ctx.node_pass_state.clear_runtime();
}

}  // namespace node_partitioning
