/****************************************************************************
 * node_partitioner_pass_state.h
 *****************************************************************************/

#ifndef NODE_PARTITIONER_PASS_STATE_H_
#define NODE_PARTITIONER_PASS_STATE_H_

#include <memory>
#include <vector>

#include "data_structure/buffered_map.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "definitions.h"


namespace partition::state {

// Mutable state scoped to one node-stream pass.
struct NodePartitionerPassState {
    int restream_number = 0;
    int curr_batch = 0;
    LongNodeID total_stream_nodecounter = 0;
    LongNodeID total_stream_nodeweight = 0;
    LongNodeID lower_global_node = 1;
    LongNodeID upper_global_node = 0;

    LongNodeID remaining_stream_nodes = 0;
    LongEdgeID remaining_stream_edges = 0;
    int remaining_stream_ew = 0;
    bool stream_input_binary = false;
    LongNodeID total_nodes_loaded = 0;

    // Owned pass-scoped runtime storage (Config stores non-owning pointers).
    std::unique_ptr<std::vector<PartitionID>> stream_nodes_assign;
    std::unique_ptr<std::vector<PartitionID>> stream_nodes_batch_marker;
    std::unique_ptr<std::vector<NodeWeight>> stream_blocks_weight;
    std::unique_ptr<std::vector<NodeWeight>> add_blocks_weight;
    std::unique_ptr<std::vector<std::vector<std::pair<NodeID, NodeWeight>>>> edge_block_nodes;
    std::unique_ptr<buffered_map> ghostglobal_to_ghostkey;
    std::unique_ptr<std::vector<NodeID>> ghostkey_to_node;
    std::unique_ptr<std::vector<std::vector<std::pair<NodeID, ShortEdgeWeight>>>> ghostkey_to_edges;
    std::unique_ptr<maxNodeHeap> stream_block_min_queue;

    void reset_from_header(LongNodeID nodes, LongEdgeID edges, int edge_weight_flag,
                           bool binary_input) {
        total_stream_nodecounter = 0;
        total_stream_nodeweight = 0;
        remaining_stream_nodes = nodes;
        remaining_stream_edges = edges;
        remaining_stream_ew = edge_weight_flag;
        stream_input_binary = binary_input;
        total_nodes_loaded = 0;
        curr_batch = 0;
        lower_global_node = 1;
        upper_global_node = 0;
    }

    void clear_runtime() {
        stream_nodes_assign.reset();
        stream_nodes_batch_marker.reset();
        stream_blocks_weight.reset();
        add_blocks_weight.reset();
        edge_block_nodes.reset();
        ghostglobal_to_ghostkey.reset();
        ghostkey_to_node.reset();
        ghostkey_to_edges.reset();
        stream_block_min_queue.reset();
    }
};

} // namespace partition::state


#endif /* NODE_PARTITIONER_PASS_STATE_H_ */
