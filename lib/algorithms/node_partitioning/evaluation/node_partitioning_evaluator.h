/****************************************************************************
 * node_partitioning_evaluator.h
 *****************************************************************************/

#ifndef NODE_PARTITIONING_EVALUATOR_H_
#define NODE_PARTITIONING_EVALUATOR_H_

#include <string>

#include "core/interfaces/algorithm_evaluator.h"
#include "partition/partition_config.h"

namespace node_partitioning {
namespace evaluation {

class NodePartitioningEvaluator : public core::interfaces::IAlgorithmEvaluator {
   public:
    void evaluate(Config& config, const std::string& graph_filename, StreamContext& ctx,
                  FlatBufferWriter& fb_writer) const override;
};

}  // namespace evaluation
}  // namespace node_partitioning

#endif /* NODE_PARTITIONING_EVALUATOR_H_ */
