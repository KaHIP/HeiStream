/******************************************************************************
 * stream_node_assign.cpp
 *****************************************************************************/

#include "partition/assignment/node_stream_assign.h"
#include "partition/assignment/node_fennel_assign.h"

namespace stream_node_assign {

void generalize_stream_partition(Config& config,
                                 partition::state::NodePartitionerPassState& pass_state,
                                 graph_access& G_local) {
    // Translate each local batch node to its global node ID and persist
    // assignment/weights into stream-global arrays.
    for (NodeID node = 0, end = config.nmbNodes; node < end; node++) {
        PartitionID block = G_local.getPartitionIndex(node);
        LongNodeID global_node = (LongNodeID)node + pass_state.lower_global_node - 1;
        if (config.local_to_global_map) {
            global_node = (*config.local_to_global_map)[node] - 1;
        }
        (*pass_state.stream_nodes_assign)[global_node] = block;
        if (config.sep_batch_marker && pass_state.stream_nodes_batch_marker) {
            (*pass_state.stream_nodes_batch_marker)[global_node] = INVALID_PARTITION;
        }
        node_fennel_assign::add_block_weight(
            config, pass_state, block,
            static_cast<NodeWeight>(G_local.getNodeWeight(node) - G_local.getImplicitGhostNodes(node)));

        if (config.ghost_neighbors_enabled && config.batch_unpartitioned_neighbors) {
            // Ghost neighbors that were deferred in the local model receive
            // soft assignments in the companion [k,2k) marker range.
            for (auto& ghost_target : (*config.batch_unpartitioned_neighbors)[node]) {
                PartitionID cur_ghost_par = (*pass_state.stream_nodes_assign)[ghost_target - 1];
                if ((cur_ghost_par >= config.k && cur_ghost_par < 2 * config.k) ||
                    cur_ghost_par == INVALID_PARTITION) {
                    (*pass_state.stream_nodes_assign)[ghost_target - 1] = block + config.k;
                }
            }
        }
    }
}

}  // namespace stream_node_assign
