/****************************************************************************
 * stream_node_direct_fennel.cpp
 *****************************************************************************/

#include <memory>

#include "algorithms/node_partitioning/execution/node_execution_internal.h"
#include "core/timing/scoped_stage_timer.h"
#include "io/node_stream_reader.h"
#include "partition/assignment/node_fennel_assign.h"

namespace stream_node_modes {

namespace {

void run_direct_fennel_node(Config& config, partition::state::NodePartitionerPassState& pass_state,
                            const std::vector<LongNodeID>& line_numbers, StreamContext& ctx) {
    // Direct no-buffer fast path: partition one stream node with Fennel.
    LongNodeID global_node_id = pass_state.total_stream_nodecounter + 1;
    pass_state.lower_global_node = global_node_id;
    pass_state.upper_global_node = global_node_id;
    pass_state.curr_batch++;

    NodeWeight node_weight = 1;
    LongEdgeID used_edges = 0;
    std::vector<LongNodeID> adjacents;
    adjacents.reserve(line_numbers.size());
    extract_fennel_neighbors(pass_state.remaining_stream_ew, line_numbers, global_node_id,
                             node_weight, adjacents, used_edges);

    if (pass_state.restream_number) {
        PartitionID prev_partition = (*pass_state.stream_nodes_assign)[global_node_id - 1];
        if (prev_partition < config.k) {
            node_fennel_assign::sub_block_weight(config, pass_state, prev_partition, node_weight);
        }
    }

    {
        TIMED_SCOPE(ctx.timing, core::timing::StreamTimingStage::Partition);
        assign_single_node_direct(config, pass_state, global_node_id, adjacents);
    }

    if (config.sep_batch_marker && pass_state.stream_nodes_batch_marker) {
        (*pass_state.stream_nodes_batch_marker)[global_node_id - 1] =
            static_cast<PartitionID>(pass_state.curr_batch);
    }

    pass_state.total_stream_nodecounter += 1;
    pass_state.total_stream_nodeweight += node_weight;
    pass_state.remaining_stream_nodes -= 1;
    pass_state.remaining_stream_edges -= used_edges;

    write_stream_progress_if_enabled(config, pass_state);
}

}  // namespace

void run_node_stream_direct_fennel(Config& config, const std::string& graph_filename,
                                   StreamContext& ctx,
                                   partition::state::NodePartitionerPassState& pass_state) {
    auto input_cfg = config.stream_input_config();
    NodeStreamReader reader(config);

    int& passes = input_cfg.num_streams_passes;
    for (pass_state.restream_number = 0; pass_state.restream_number < passes;
         ++pass_state.restream_number) {
        ctx.timing.add(core::timing::StreamTimingStage::ReadInput,
                       begin_node_stream_pass(reader, pass_state, graph_filename));

        std::vector<LongNodeID> line_numbers;
        line_numbers.reserve(64);
        while (pass_state.remaining_stream_nodes) {
            config.nmbNodes = 1;
            const double io_before = reader.disk_read_time();
            if (!reader.next_node_line(line_numbers)) {
                break;
            }
            ctx.timing.add_delta(core::timing::StreamTimingStage::ReadInput, io_before,
                                 reader.disk_read_time());
            run_direct_fennel_node(config, pass_state, line_numbers, ctx);
        }

        reader.end_pass();
    }
}

}  // namespace stream_node_modes
