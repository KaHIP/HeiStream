/******************************************************************************
 * stream_edge_driver.h
 *****************************************************************************/

#ifndef STREAM_EDGE_DRIVER_H_
#define STREAM_EDGE_DRIVER_H_

#include <string>

#include "core/context/stream_context.h"
#include "partition/partition_config.h"
#include "partition/state/edge_partitioner_pass_state.h"

// Executes the edge-stream partitioning pipeline end-to-end.
void run_edge_stream(Config& config, const std::string& graph_filename, StreamContext& ctx,
                     partition::state::EdgePartitionerPassState& pass_state);

#endif /* STREAM_EDGE_DRIVER_H_ */
