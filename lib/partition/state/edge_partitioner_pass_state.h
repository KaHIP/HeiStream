/****************************************************************************
 * edge_partitioner_pass_state.h
 *****************************************************************************/

#ifndef EDGE_PARTITIONER_PASS_STATE_H_
#define EDGE_PARTITIONER_PASS_STATE_H_

#include <fstream>
#include <memory>
#include <sparsehash/dense_hash_set>
#include <vector>

#include "definitions.h"


namespace partition::state {

// Mutable state that is scoped to one edge-stream pass.
struct EdgePartitionerPassState {
    int restream_number = 0;
    int curr_batch = 0;
    LongNodeID total_stream_nodecounter = 0;
    LongNodeID total_stream_nodeweight = 0;
    LongNodeID lower_global_node = 1;
    LongNodeID upper_global_node = 0;

    bool stream_input_binary = false;
    LongEdgeID remaining_stream_partition_nodes = 0;
    LongEdgeID remaining_stream_nodes_og = 0;
    LongEdgeID remaining_stream_graph_nodes = 0;
    LongNodeID lower_global_node_conv = 0;
    LongNodeID upper_global_node_conv = 0;
    NodeID incremental_edge_id = 0;
    NodeID prev_batch_edge_id = 0;
    NodeID last_edge_count = 0;

    // Owned pass-scoped runtime storage (Config stores non-owning pointers).
    std::unique_ptr<std::ofstream> stream_out;
    std::unique_ptr<std::vector<PartitionID>> stream_nodes_assign;
    std::unique_ptr<std::vector<NodeWeight>> stream_blocks_weight;
    std::unique_ptr<std::vector<NodeWeight>> add_blocks_weight;
    std::unique_ptr<std::vector<std::vector<std::pair<NodeID, NodeWeight>>>> edge_block_nodes;
    std::unique_ptr<std::vector<std::vector<NodeID>>> nodes_on_edge_conv;
    std::unique_ptr<std::vector<google::dense_hash_set<PartitionID>>> blocks_on_node;
    std::unique_ptr<std::vector<NodeID>> blocks_on_node_minimal;

    void reset_for_pass(LongNodeID original_nodes, LongEdgeID original_edges, bool binary_input) {
        total_stream_nodecounter = 0;
        total_stream_nodeweight = 0;
        stream_input_binary = binary_input;
        remaining_stream_partition_nodes = original_edges;
        remaining_stream_nodes_og = original_nodes;
        remaining_stream_graph_nodes = original_nodes;
        lower_global_node_conv = 0;
        upper_global_node_conv = 0;
        incremental_edge_id = 0;
        prev_batch_edge_id = 0;
        last_edge_count = 0;
        curr_batch = 0;
        lower_global_node = 1;
        upper_global_node = 0;
    }

    void clear_runtime() {
        stream_out.reset();
        stream_nodes_assign.reset();
        stream_blocks_weight.reset();
        add_blocks_weight.reset();
        edge_block_nodes.reset();
        nodes_on_edge_conv.reset();
        blocks_on_node.reset();
        blocks_on_node_minimal.reset();
    }
};

} // namespace partition::state


#endif /* EDGE_PARTITIONER_PASS_STATE_H_ */
