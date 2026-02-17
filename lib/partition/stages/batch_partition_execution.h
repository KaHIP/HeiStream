/******************************************************************************
 * stream_driver_helpers.h
 *****************************************************************************/

#ifndef STREAM_DRIVER_HELPERS_H_
#define STREAM_DRIVER_HELPERS_H_

#include <string>

#include "partition/partition_config.h"
#include "partition/state/edge_partitioner_pass_state.h"
#include "partition/state/node_partitioner_pass_state.h"

class graph_access;
class balance_configuration;

// Run node-batch partitioning stage (including balance setup).
void perform_node_batch_partition(Config& partition_config,
                                  partition::state::NodePartitionerPassState& pass_state,
                                  graph_access& G, balance_configuration& bc);

// Run node-batch postprocess stage that writes local IDs back to global stream IDs.
void postprocess_node_batch_partition(Config& partition_config,
                                      partition::state::NodePartitionerPassState& pass_state,
                                      graph_access& G);

// Run common post-model edge batch steps and accumulate partition/postprocess timers.
void perform_edge_batch_partition(Config& partition_config,
                                  partition::state::EdgePartitionerPassState& pass_state,
                                  graph_access& G, balance_configuration& bc);

// Run edge-batch postprocess stage that writes local IDs back to global stream IDs.
void postprocess_edge_batch_partition(Config& partition_config,
                                      partition::state::EdgePartitionerPassState& pass_state,
                                      graph_access& G);

// Assign an isolated node to the lightest block.
void assign_isolated_node(Config& partition_config,
                          partition::state::NodePartitionerPassState& pass_state,
                          LongNodeID global_node_id);

#endif /* STREAM_DRIVER_HELPERS_H_ */
