/****************************************************************************
 * edge_partitioning_algorithm_factory.cpp
 *****************************************************************************/

#include "algorithms/edge_partitioning/factory/edge_partitioning_algorithm_factory.h"

#include <memory>

#include "algorithms/edge_partitioning/edge_partitioning_algorithm.h"
#include "algorithms/edge_partitioning/evaluation/edge_partitioning_evaluator.h"
#include "algorithms/edge_partitioning/reporting/edge_partitioning_reporter.h"
#include "algorithms/edge_partitioning/scheduler/edge_partitioning_scheduler.h"
#include "algorithms/edge_partitioning/setup/edge_partitioning_setup.h"


namespace edge_partitioning::factory {

std::unique_ptr<IStreamAlgorithm> create_algorithm() {
    return std::unique_ptr<IStreamAlgorithm>(new EdgePartitioningAlgorithm(
        std::unique_ptr<core::interfaces::IAlgorithmSetup>(new setup::EdgePartitioningSetup()),
        std::unique_ptr<core::interfaces::IAlgorithmScheduler>(
            new scheduler::EdgePartitioningScheduler()),
        std::unique_ptr<core::interfaces::IAlgorithmEvaluator>(
            new evaluation::EdgePartitioningEvaluator()),
        std::unique_ptr<core::interfaces::IAlgorithmReporter>(
            new reporting::EdgePartitioningReporter())));
}

} // namespace edge_partitioning::factory

