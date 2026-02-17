/******************************************************************************
 * stream_edge_assign.cpp
 *****************************************************************************/

#include "partition/assignment/edge_stream_assign.h"

namespace stream_edge_assign {

void generalize_stream_partition(Config& config,
                                 partition::state::EdgePartitionerPassState& pass_state,
                                 graph_access& G_local) {
    // Merge local edge-batch assignments into global stream structures.
    for (NodeID node = 0, end = config.nmbNodes; node < end; node++) {
        PartitionID block = G_local.getPartitionIndex(node);
        LongNodeID global_node = (LongNodeID)node + pass_state.lower_global_node - 1;

        if (pass_state.stream_nodes_assign) {
            (*pass_state.stream_nodes_assign)[global_node] = block;
        }
        (*pass_state.stream_blocks_weight)[block] +=
            G_local.getNodeWeight(node) - G_local.getImplicitGhostNodes(node);

        std::vector<NodeID> incident_node_list;
        incident_node_list = (*pass_state.nodes_on_edge_conv)[node];

        if (config.past_subset_size != 0) {
            // Track which edge partitions touch each original node; this is
            // used later to compute replication and vertex-cut statistics.
            for (auto& incident_node : incident_node_list) {
                if (config.minimal_mode == false) {
                    if ((*pass_state.blocks_on_node)[incident_node].empty() == true) {
                        (*pass_state.blocks_on_node)[incident_node].set_empty_key(-1);
                        (*pass_state.blocks_on_node)[incident_node].insert(block);
                        continue;
                    }
                    if ((*pass_state.blocks_on_node)[incident_node].empty() == false) {
                        if ((*pass_state.blocks_on_node)[incident_node].find(block) ==
                            (*pass_state.blocks_on_node)[incident_node].end()) {
                            (*pass_state.blocks_on_node)[incident_node].insert(block);
                        }
                    }
                } else {
                    (*pass_state.blocks_on_node_minimal)[incident_node] = block;
                }
            }
        }
    }
}

}  // namespace stream_edge_assign
