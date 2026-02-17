/******************************************************************************
 * stream_edge_assign.h
 *****************************************************************************/

#ifndef STREAM_EDGE_ASSIGN_H_
#define STREAM_EDGE_ASSIGN_H_

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "partition/state/edge_partitioner_pass_state.h"

namespace stream_edge_assign {

// Applies a locally computed edge-batch partition back to global stream state.
void generalize_stream_partition(Config& config,
                                 partition::state::EdgePartitionerPassState& pass_state,
                                 graph_access& G_local);

}  // namespace stream_edge_assign

#endif /* STREAM_EDGE_ASSIGN_H_ */
