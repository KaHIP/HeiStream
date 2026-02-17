/****************************************************************************
 * edge_partitioning_evaluator.h
 *****************************************************************************/

#ifndef EDGE_PARTITIONING_EVALUATOR_H_
#define EDGE_PARTITIONING_EVALUATOR_H_

#include <string>

#include "core/interfaces/algorithm_evaluator.h"
#include "partition/partition_config.h"

namespace edge_partitioning {
namespace evaluation {

class EdgePartitioningEvaluator : public core::interfaces::IAlgorithmEvaluator {
   public:
    void evaluate(Config& config, const std::string& graph_filename, StreamContext& ctx,
                  FlatBufferWriter& fb_writer) const override;
};

}  // namespace evaluation
}  // namespace edge_partitioning

#endif /* EDGE_PARTITIONING_EVALUATOR_H_ */
