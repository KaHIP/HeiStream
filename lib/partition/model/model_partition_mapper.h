/****************************************************************************
 * model_partition_mapper.h
 *****************************************************************************/

#ifndef PARTITION_MODEL_MODEL_PARTITION_MAPPER_H_
#define PARTITION_MODEL_MODEL_PARTITION_MAPPER_H_

#include <vector>

#include "definitions.h"
#include "partition/partition_config.h"

namespace partition_model {

/**
 * Resolve the partition assignment to attach to a model node.
 *
 * @param config Immutable run configuration.
 * @param stream_nodes_assign Optional global node-to-partition assignments.
 *        Required only for restream recovery paths.
 * @param lower_global_node Global id corresponding to local node 0.
 * @param node Local model node id.
 * @param node_counter Total node count in the materialized model.
 * @return Partition id to assign to the model node.
 */
PartitionID recover_partition_for_model_node(const Config& config,
                                             const std::vector<PartitionID>* stream_nodes_assign,
                                             LongNodeID lower_global_node, NodeID node,
                                             NodeID node_counter, int restream_number);

}  // namespace partition_model

#endif /* PARTITION_MODEL_MODEL_PARTITION_MAPPER_H_ */
