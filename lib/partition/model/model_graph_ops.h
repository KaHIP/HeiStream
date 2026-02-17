/****************************************************************************
 * model_graph_ops.h
 *****************************************************************************/

#ifndef PARTITION_MODEL_MODEL_GRAPH_OPS_H_
#define PARTITION_MODEL_MODEL_GRAPH_OPS_H_

#include <functional>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "partition/partition_config.h"


namespace partition_model::graph_ops {

/**
 * Insert an undirected edge into batch adjacency storage.
 *
 * @param all_edges Per-node adjacency storage.
 * @param node Source node.
 * @param target Target node.
 * @param edge_weight Edge weight.
 * @return Number of directed edges inserted (always 2).
 */
EdgeID include_edge(std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges, NodeID node,
                    NodeID target, EdgeWeight edge_weight);

/**
 * Insert one regular in-batch edge based on global target id.
 *
 * @param config Mutable run configuration.
 * @param all_edges Per-node adjacency storage.
 * @param node Local source node id.
 * @param global_target Global stream target node id.
 * @param lower_global_node Global id of local node 0 in this batch.
 * @param edge_weight Edge weight.
 * @return Number of directed edges inserted.
 */
EdgeID insert_regular_edge(Config& config,
                           std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges,
                           NodeID node, LongNodeID global_target, LongNodeID lower_global_node,
                           EdgeWeight edge_weight);

/**
 * Append quotient nodes representing current block state.
 *
 * @param config Mutable run configuration.
 * @param stream_blocks_weight Current block weights.
 * @param all_nodes Output node-weight vector to append to.
 * @param base_nodes Number of non-quotient nodes currently in `all_nodes`.
 * @param node_counter Running node counter (advanced in-place).
 */
void insert_quotient_nodes(Config& config, const std::vector<NodeWeight>& stream_blocks_weight,
                           std::vector<NodeWeight>& all_nodes, NodeID base_nodes,
                           NodeID& node_counter);

/**
 * Insert edges from regular nodes to quotient nodes.
 *
 * @param config Mutable run configuration.
 * @param edge_block_nodes Per-block incident node/weight lists.
 * @param all_edges Per-node adjacency storage to append to.
 * @param base_nodes Number of non-quotient nodes.
 * @return Number of directed edges inserted.
 */
EdgeID insert_quotient_edges(
    Config& config, const std::vector<std::vector<std::pair<NodeID, NodeWeight>>>& edge_block_nodes,
    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges, NodeID base_nodes);

/**
 * Materialize a model graph from adjacency/weight vectors.
 *
 * @param graph Destination graph object.
 * @param node_count Number of nodes to materialize.
 * @param edge_count Number of directed edges expected.
 * @param all_edges Per-node adjacency lists.
 * @param all_nodes Node weights.
 * @param implicit_ghost_nodes Per-node implicit ghost counts.
 * @param partition_for_node Callback returning initial partition for each node.
 */
void materialize_graph(graph_access& graph, NodeID node_count, EdgeID edge_count,
                       const std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges,
                       const std::vector<NodeWeight>& all_nodes,
                       const std::vector<NodeWeight>& implicit_ghost_nodes,
                       const std::function<PartitionID(NodeID)>& partition_for_node);

} // namespace partition_model::graph_ops


#endif /* PARTITION_MODEL_MODEL_GRAPH_OPS_H_ */
