/****************************************************************************
 * edge_partitioning_algorithm.h
 *****************************************************************************/

#ifndef EDGE_PARTITIONING_ALGORITHM_H_
#define EDGE_PARTITIONING_ALGORITHM_H_

#include <memory>
#include <string>

#include "core/interfaces/algorithm_evaluator.h"
#include "core/interfaces/algorithm_reporter.h"
#include "core/interfaces/algorithm_scheduler.h"
#include "core/interfaces/algorithm_setup.h"
#include "core/interfaces/stream_algorithm.h"

namespace edge_partitioning {

// Edge-partitioning algorithm adapter for runtime dispatch.
class EdgePartitioningAlgorithm : public IStreamAlgorithm {
   public:
    EdgePartitioningAlgorithm(std::unique_ptr<core::interfaces::IAlgorithmSetup> setup,
                              std::unique_ptr<core::interfaces::IAlgorithmScheduler> scheduler,
                              std::unique_ptr<core::interfaces::IAlgorithmEvaluator> evaluator,
                              std::unique_ptr<core::interfaces::IAlgorithmReporter> reporter);

    void prepare_run(Config& config, const std::string& graph_filename,
                     StreamContext& ctx) override;
    void run(Config& config, const std::string& graph_filename, StreamContext& ctx) override;
    void finalize_run(Config& config, const std::string& graph_filename,
                      StreamContext& ctx) override;

   private:
    std::unique_ptr<core::interfaces::IAlgorithmSetup> setup_;
    std::unique_ptr<core::interfaces::IAlgorithmScheduler> scheduler_;
    std::unique_ptr<core::interfaces::IAlgorithmEvaluator> evaluator_;
    std::unique_ptr<core::interfaces::IAlgorithmReporter> reporter_;
};

}  // namespace edge_partitioning

#endif /* EDGE_PARTITIONING_ALGORITHM_H_ */
