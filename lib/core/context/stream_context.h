/****************************************************************************
 * stream_context.h
 *****************************************************************************/

#ifndef CORE_STREAM_CONTEXT_H_
#define CORE_STREAM_CONTEXT_H_

#include "core/timing/stream_timing_recorder.h"
#include "definitions.h"
#include "partition/partition_config.h"
#include "partition/state/edge_partitioner_pass_state.h"
#include "partition/state/node_partitioner_pass_state.h"

struct StreamCounts {
    LongNodeID total_nodes = 0;
    LongEdgeID total_edges = 0;
    LongNodeID remaining_nodes = 0;
    LongEdgeID remaining_edges = 0;

    // Snapshot mutable counters from partition config.
    void sync_from_config(const Config& config) {
        total_nodes = config.total_nodes;
        total_edges = config.total_edges;
        remaining_nodes = config.total_nodes - config.stream_assigned_nodes;
        remaining_edges = 0;
    }
};

struct StreamContext {
    Config& config;
    StreamMode mode;
    StreamCounts counts;
    core::timing::StreamTimingRecorder timing;
    long max_rss_kb = 0;
    partition::state::NodePartitionerPassState node_pass_state;
    partition::state::EdgePartitionerPassState edge_pass_state;

    StreamContext(Config& config_in, StreamMode mode_in) : config(config_in), mode(mode_in) {}
};

#endif /* CORE_STREAM_CONTEXT_H_ */
