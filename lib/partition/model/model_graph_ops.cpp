/****************************************************************************
 * model_graph_ops.cpp
 *****************************************************************************/

#include "partition/model/model_graph_ops.h"

#include <iostream>


namespace partition_model::graph_ops {

EdgeID include_edge(std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges, NodeID node,
                    NodeID target, EdgeWeight edge_weight) {
    // Stream models are stored as undirected graphs.
    all_edges[node].emplace_back(target, edge_weight);
    all_edges[target].emplace_back(node, edge_weight);
    return 2;
}

EdgeID insert_regular_edge(Config& config,
                           std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges,
                           NodeID node, LongNodeID global_target, LongNodeID lower_global_node,
                           EdgeWeight edge_weight) {
    // Convert from global stream ID to current batch-local ID.
    NodeID target = (NodeID)(global_target - lower_global_node);
    if (target == node) {
        std::cerr << "The graph file contains self-loops, which are not supported. "
                     "Please remove them from the file."
                  << '\n';
        exit(0);
    }
    if (target > node) {
        // Store each undirected edge once; symmetric edge is added in include_edge.
        return 0;
    }
    return include_edge(all_edges, node, target, (1 + config.double_non_ghost_edges) * edge_weight);
}

void insert_quotient_nodes(Config& config, const std::vector<NodeWeight>& stream_blocks_weight,
                           std::vector<NodeWeight>& all_nodes, NodeID base_nodes,
                           NodeID& node_counter) {
    // Quotient node weight reflects currently assigned weight of that block.
    while (node_counter < base_nodes + config.quotient_nodes) {
        NodeID node = node_counter++;
        NodeID targetPar = node - base_nodes;
        NodeWeight weight = stream_blocks_weight[targetPar];
        all_nodes[node] = weight;
    }
}

EdgeID insert_quotient_edges(
    Config& config, const std::vector<std::vector<std::pair<NodeID, NodeWeight>>>& edge_block_nodes,
    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges, NodeID base_nodes) {
    EdgeID inserted_edges = 0;
    for (PartitionID block = 0; block < config.k; block++) {
        NodeID target = base_nodes + block;
        if (edge_block_nodes[block].size() < 1) {
            continue;
        }
        // Connect each batch node to incident block summary.
        for (auto& element : edge_block_nodes[block]) {
            inserted_edges += include_edge(all_edges, element.first, target,
                                           (1 + config.double_non_ghost_edges) * element.second);
        }
    }
    return inserted_edges;
}

void materialize_graph(graph_access& graph, NodeID node_count, EdgeID edge_count,
                       const std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges,
                       const std::vector<NodeWeight>& all_nodes,
                       const std::vector<NodeWeight>& implicit_ghost_nodes,
                       const std::function<PartitionID(NodeID)>& partition_for_node) {
    graph.stream_repeat_construction(node_count, edge_count);
    graph.resizeImplicitGhostNodesList(node_count);
    std::vector<EdgeWeight> neighbor_edgeweight(node_count, 0);
    std::vector<EdgeID> neighbor_edgeid(node_count, 0);
    std::vector<NodeID> neighbors;

    for (NodeID node = 0; node < node_count; ++node) {
        NodeID gnode = graph.new_node();
        graph.setNodeWeight(gnode, all_nodes[node]);
        graph.setImplicitGhostNodes(gnode, implicit_ghost_nodes[node]);
        graph.setPartitionIndex(gnode, partition_for_node(node));

        for (auto& [target, edge_weight] : all_edges[node]) {
            EdgeID e;
            if (neighbor_edgeweight[target] == 0) {
                e = graph.new_edge(gnode, target);
                neighbor_edgeid[target] = e;
                neighbors.push_back(target);
            } else {
                e = neighbor_edgeid[target];
            }
            neighbor_edgeweight[target] += edge_weight;
            graph.setEdgeWeight(e, neighbor_edgeweight[target]);
        }
        for (auto target : neighbors) {
            neighbor_edgeweight[target] = 0;
        }
        neighbors.clear();
    }

    graph.stream_finish_construction();
}

} // namespace partition_model::graph_ops

