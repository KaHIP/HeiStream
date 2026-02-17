/******************************************************************************
 * stream_edge_driver.cpp
 *****************************************************************************/

#include "algorithms/edge_partitioning/execution/edge_execution.h"

#include <algorithm>
#include <memory>

#include "algorithms/edge_partitioning/model/edge_batch_model_builder.h"
#include "algorithms/edge_partitioning/setup/edge_partitioning_setup.h"
#include "core/context/stream_context.h"
#include "core/timing/scoped_stage_timer.h"
#include "data_structure/graph_access.h"
#include "io/edge_stream_reader.h"
#include "io/stream_reader.h"
#include "partition/balance_configuration.h"
#include "partition/stages/batch_partition_execution.h"

namespace {

void append_edge_batch_assignments(const Config& config,
                                   partition::state::EdgePartitionerPassState& pass_state,
                                   graph_access& graph) {
    for (NodeID node = 0, end = config.nmbNodes; node < end; ++node) {
        (*pass_state.stream_out) << graph.getPartitionIndex(node) << "\n";
    }
}

}  // namespace

void run_edge_stream(Config& config, const std::string& graph_filename, StreamContext& ctx,
                     partition::state::EdgePartitionerPassState& pass_state) {
    auto input_cfg = config.stream_input_config();
    const auto out_cfg = config.stream_output_config();
    balance_configuration bc;

    std::unique_ptr<graph_access> G(new graph_access());
    io::stream::StreamReader stream_reader;
    edge_partitioning::setup::EdgePartitioningSetup setup;

    if (input_cfg.num_streams_passes > 1 && !config.minimal_mode) {
        std::cerr << "Edge restreaming is currently supported only with --minimal_mode=true."
                  << '\n';
        exit(1);
    }

    // Every stream pass rebuilds the split-graph batches from input order.
    int& passes = input_cfg.num_streams_passes;
    for (pass_state.restream_number = 0; pass_state.restream_number < passes;
         ++pass_state.restream_number) {
        const double header_io_before = stream_reader.disk_read_time();
        io::stream::StreamSourceConfig source_config;
        source_config.graph_filename = graph_filename;
        source_config.allow_binary = true;
        const io::stream::StreamHeader header = stream_reader.read_header(source_config);
        ctx.timing.add_delta(core::timing::StreamTimingStage::ReadInput, header_io_before,
                             stream_reader.disk_read_time());
        setup.initialize_pass(config, pass_state, header);

        while (pass_state.remaining_stream_partition_nodes) {
            // Batch granularity is defined over remaining original graph nodes.
            config.nmbNodes =
                std::min(input_cfg.batch_size, pass_state.remaining_stream_graph_nodes);

            G->set_partition_count(config.k);
            const double io_before = stream_reader.disk_read_time();
            {
                TIMED_SCOPE(ctx.timing, core::timing::StreamTimingStage::ModelBuild);
                edge_partitioning::model::build_batch_model(config, pass_state, *G, stream_reader);
            }
            const double io_after = stream_reader.disk_read_time();
            ctx.timing.add_delta(core::timing::StreamTimingStage::ReadInput, io_before, io_after);
            // The model stage currently includes stream reads; remove read delta
            // from model time so stage accounting matches node pipeline semantics.
            ctx.timing.add(core::timing::StreamTimingStage::ModelBuild, -(io_after - io_before));

            {
                TIMED_SCOPE(ctx.timing, core::timing::StreamTimingStage::Partition);
                perform_edge_batch_partition(config, pass_state, *G, bc);
            }
            {
                TIMED_SCOPE(ctx.timing, core::timing::StreamTimingStage::Postprocess);
                postprocess_edge_batch_partition(config, pass_state, *G);
            }
            if (out_cfg.stream_output_progress) {
                append_edge_batch_assignments(config, pass_state, *G);
            }
        }

        pass_state.lower_global_node = 0;
        pass_state.upper_global_node = 0;
    }
    // graph object is released automatically
}
