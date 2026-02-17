/******************************************************************************
 * stream_node_assign.h
 *****************************************************************************/

#ifndef STREAM_NODE_ASSIGN_H_
#define STREAM_NODE_ASSIGN_H_

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "partition/state/node_partitioner_pass_state.h"

namespace stream_node_assign {

// Applies a locally computed batch partition back to global stream state.
void generalize_stream_partition(Config& config,
                                 partition::state::NodePartitionerPassState& pass_state,
                                 graph_access& G_local);

}  // namespace stream_node_assign

#endif /* STREAM_NODE_ASSIGN_H_ */
