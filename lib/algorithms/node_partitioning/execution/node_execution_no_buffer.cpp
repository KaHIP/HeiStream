/****************************************************************************
 * stream_node_mode_no_buffer.cpp
 *****************************************************************************/

#include <memory>

#include "algorithms/node_partitioning/execution/node_execution_internal.h"
#include "algorithms/node_partitioning/model/node_batch_model_builder.h"
#include "core/timing/scoped_stage_timer.h"
#include "io/node_stream_reader.h"
#include "partition/balance_configuration.h"
#include "partition/stages/batch_partition_execution.h"

namespace stream_node_modes {

void run_node_stream_no_priority_buffer(Config& config, const std::string& graph_filename,
                                        StreamContext& ctx,
                                        partition::state::NodePartitionerPassState& pass_state) {
    auto input_cfg = config.stream_input_config();
    balance_configuration bc;
    NodeStreamReader reader(config);
    std::unique_ptr<graph_access> G(new graph_access());

    int& passes = input_cfg.num_streams_passes;
    for (pass_state.restream_number = 0; pass_state.restream_number < passes;
         ++pass_state.restream_number) {
        // Baseline stream path without priority buffering.
        ctx.timing.add(core::timing::StreamTimingStage::ReadInput,
                       begin_node_stream_pass(reader, pass_state, graph_filename));

        std::vector<std::vector<LongNodeID>> batch_lines;
        while (pass_state.remaining_stream_nodes != 0U) {
            LongNodeID requested_batch =
                std::min(input_cfg.batch_size, pass_state.remaining_stream_nodes);
            config.nmbNodes = requested_batch;
            const double io_before = reader.disk_read_time();
            if (!reader.next_batch(requested_batch, batch_lines)) {
                break;
            }
            ctx.timing.add_delta(core::timing::StreamTimingStage::ReadInput, io_before,
                                 reader.disk_read_time());

            G->set_partition_count(config.k);
            {
                TIMED_SCOPE(ctx.timing, core::timing::StreamTimingStage::ModelBuild);
                node_partitioning::model::build_batch_model_from_stream_input(config, pass_state,
                                                                              *G, batch_lines);
            }

            {
                TIMED_SCOPE(ctx.timing, core::timing::StreamTimingStage::Partition);
                perform_node_batch_partition(config, pass_state, *G, bc);
            }
            {
                TIMED_SCOPE(ctx.timing, core::timing::StreamTimingStage::Postprocess);
                postprocess_node_batch_partition(config, pass_state, *G);
            }
            write_stream_progress_if_enabled(config, pass_state);
        }

        reader.end_pass();
    }
}

}  // namespace stream_node_modes
