/****************************************************************************
 * node_partitioning_algorithm_factory.cpp
 *****************************************************************************/

#include "algorithms/node_partitioning/factory/node_partitioning_algorithm_factory.h"

#include <memory>

#include "algorithms/node_partitioning/evaluation/node_partitioning_evaluator.h"
#include "algorithms/node_partitioning/node_partitioning_algorithm.h"
#include "algorithms/node_partitioning/reporting/node_partitioning_reporter.h"
#include "algorithms/node_partitioning/scheduler/node_partitioning_scheduler.h"
#include "algorithms/node_partitioning/setup/node_partitioning_setup.h"


namespace node_partitioning::factory {

std::unique_ptr<IStreamAlgorithm> create_algorithm() {
    return std::unique_ptr<IStreamAlgorithm>(new NodePartitioningAlgorithm(
        std::unique_ptr<core::interfaces::IAlgorithmSetup>(new setup::NodePartitioningSetup()),
        std::unique_ptr<core::interfaces::IAlgorithmScheduler>(
            new scheduler::NodePartitioningScheduler()),
        std::unique_ptr<core::interfaces::IAlgorithmEvaluator>(
            new evaluation::NodePartitioningEvaluator()),
        std::unique_ptr<core::interfaces::IAlgorithmReporter>(
            new reporting::NodePartitioningReporter())));
}

} // namespace node_partitioning::factory

