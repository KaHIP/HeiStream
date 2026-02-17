/******************************************************************************
 * stream_driver_helpers.cpp
 *****************************************************************************/

#include "partition/stages/batch_partition_execution.h"

#include "core/timing/scoped_stage_timer.h"
#include "data_structure/graph_access.h"
#include "partition/assignment/edge_stream_assign.h"
#include "partition/assignment/node_fennel_assign.h"
#include "partition/assignment/node_stream_assign.h"
#include "partition/balance_configuration.h"
#include "partition/heuristics/fennel_scoring.h"
#include "partition/stages/stream_partition_stages.h"

void perform_node_batch_partition(Config& partition_config,
                                  partition::state::NodePartitionerPassState& pass_state,
                                  graph_access& G, balance_configuration& bc) {
    partition_stages::prepare_batch_partition(partition_config, pass_state.restream_number,
                                              pass_state.curr_batch,
                                              *pass_state.stream_blocks_weight, G, bc);
    // Keep stream queue acceleration orthogonal to KaHIP internal partition search.
    const bool saved_use_queue = partition_config.use_queue;
    partition_config.use_queue = false;
    partition_stages::execute_partition(partition_config, G, pass_state.restream_number);
    partition_config.use_queue = saved_use_queue;
}

void postprocess_node_batch_partition(Config& partition_config,
                                      partition::state::NodePartitionerPassState& pass_state,
                                      graph_access& G) {
    // Push local batch assignment back to global stream assignment vectors.
    stream_node_assign::generalize_stream_partition(partition_config, pass_state, G);
}

void perform_edge_batch_partition(Config& partition_config,
                                  partition::state::EdgePartitionerPassState& pass_state,
                                  graph_access& G, balance_configuration& bc) {
    partition_stages::prepare_batch_partition(partition_config, pass_state.restream_number,
                                              pass_state.curr_batch,
                                              *pass_state.stream_blocks_weight, G, bc);
    partition_stages::execute_partition(partition_config, G, pass_state.restream_number);
}

void postprocess_edge_batch_partition(Config& partition_config,
                                      partition::state::EdgePartitionerPassState& pass_state,
                                      graph_access& G) {
    stream_edge_assign::generalize_stream_partition(partition_config, pass_state, G);
}

void assign_isolated_node(Config& partition_config,
                          partition::state::NodePartitionerPassState& pass_state,
                          LongNodeID global_node_id) {
    // Isolated nodes are placed greedily to maintain balance.
    const PartitionID best_partition = partition::heuristics::fennel::find_min_feasible_block(
        *pass_state.stream_blocks_weight, partition_config.max_block_weight);
    (*pass_state.stream_nodes_assign)[global_node_id - 1] = best_partition;
    node_fennel_assign::add_block_weight(partition_config, pass_state, best_partition, 1);
}
