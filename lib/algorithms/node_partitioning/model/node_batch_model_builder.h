/****************************************************************************
 * node_batch_model_builder.h
 *****************************************************************************/

#ifndef ALGORITHMS_NODE_PARTITIONING_MODEL_NODE_BATCH_MODEL_BUILDER_H_
#define ALGORITHMS_NODE_PARTITIONING_MODEL_NODE_BATCH_MODEL_BUILDER_H_

#include <cstddef>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "partition/partition_config.h"
#include "partition/state/node_partitioner_pass_state.h"

namespace node_partitioning {
namespace model {

/**
 * Build a node-partitioning model graph from stream-order batch lines.
 *
 * @param config Mutable run configuration.
 * @param pass_state Mutable node pass state.
 * @param graph Output model graph for the current batch.
 * @param input_lines Parsed stream lines for current batch nodes.
 * @return Number of model nodes materialized in `graph`.
 */
NodeID build_batch_model_from_stream_input(Config& config,
                                           partition::state::NodePartitionerPassState& pass_state,
                                           graph_access& graph,
                                           std::vector<std::vector<LongNodeID>>& input_lines);

/**
 * Build a node-partitioning model graph from a reordered batch.
 *
 * @param config Mutable run configuration.
 * @param pass_state Mutable node pass state.
 * @param graph Output model graph for the current batch.
 * @param batch_nodes Reordered `(global_node_id, adjacency_line)` records.
 * @param batch_id Logical batch identifier.
 * @return Number of model nodes materialized in `graph`.
 */
NodeID build_batch_model_from_reordered_nodes(
    Config& config, partition::state::NodePartitionerPassState& pass_state, graph_access& graph,
    std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>& batch_nodes, size_t batch_id);

}  // namespace model
}  // namespace node_partitioning

#endif /* ALGORITHMS_NODE_PARTITIONING_MODEL_NODE_BATCH_MODEL_BUILDER_H_ */
