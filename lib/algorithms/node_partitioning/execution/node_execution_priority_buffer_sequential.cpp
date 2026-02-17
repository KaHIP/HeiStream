/****************************************************************************
 * stream_node_mode_priority_buffer_sequential.cpp
 *****************************************************************************/

#include <memory>
#include <utility>
#include <vector>

#include "algorithms/node_partitioning/execution/batch_id_manager.h"
#include "algorithms/node_partitioning/execution/node_execution_internal.h"
#include "core/timing/scoped_stage_timer.h"
#include "data_structure/buffer.h"
#include "io/node_stream_reader.h"

namespace stream_node_modes {

namespace {

using BatchNodes = std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>;

struct SequentialBatchState {
    size_t current_batch_id;
    PartitionID current_batch_marker;
    BatchNodes current_batch;
};

struct BufferMaintenanceTotals {
    double add_node_time = 0.0;
    double extract_time = 0.0;
    double update_adj_time = 0.0;
};

struct FirstPassState {
    Buffer& buffer;
    SequentialBatchState& batch_state;
    bool always_top_extraction;
};

struct RestreamPassState {
    SequentialBatchState& batch_state;
};

void advance_batch(Config& config, SequentialBatchState& batch_state) {
    advance_batch_window(config, batch_state.current_batch_id, batch_state.current_batch_marker);
}

void flush_current_batch_mlp(Config& config, StreamContext& ctx,
                             partition::state::NodePartitionerPassState& pass_state,
                             SequentialBatchState& batch_state) {
    if (batch_state.current_batch.empty()) {
        return;
    }

    config.nmbNodes = batch_state.current_batch.size();
    double local_model_time = 0.0;
    double local_partition_time = 0.0;
    double local_postprocess_time = 0.0;
    perform_mlp_on_batch(config, pass_state, batch_state.current_batch, batch_state.current_batch_id,
                         local_model_time, local_partition_time, local_postprocess_time);

    ctx.timing.add(core::timing::StreamTimingStage::ModelBuild, local_model_time);
    ctx.timing.add(core::timing::StreamTimingStage::Partition, local_partition_time);
    ctx.timing.add(core::timing::StreamTimingStage::Postprocess, local_postprocess_time);

    batch_state.current_batch.clear();
    advance_batch(config, batch_state);
}

void partition_single_node_with_timing(Config& config, StreamContext& ctx,
                                       partition::state::NodePartitionerPassState& pass_state,
                                       LongNodeID node_id, std::vector<LongNodeID>& neighbors) {
    TIMED_SCOPE(ctx.timing, core::timing::StreamTimingStage::Partition);
    assign_single_node_direct(config, pass_state, node_id, neighbors);
}

void add_top_buffer_node_to_batch(Config& config,
                                  partition::state::NodePartitionerPassState& pass_state,
                                  FirstPassState& first_state) {
    BatchNode top_node = extract_top_buffer_node_for_batch(config, pass_state, first_state.buffer,
                                                           first_state.batch_state.current_batch_marker,
                                                           false);
    first_state.batch_state.current_batch.emplace_back(std::move(top_node));
}

bool read_next_stream_node(NodeStreamReader& reader, StreamContext& ctx,
                           partition::state::NodePartitionerPassState& pass_state,
                           std::vector<LongNodeID>& neighbors, LongNodeID& global_node_id) {
    double io_delta = 0.0;
    if (!read_next_node_from_stream(reader, pass_state, neighbors, global_node_id, io_delta)) {
        return false;
    }
    ctx.timing.add(core::timing::StreamTimingStage::ReadInput, io_delta);
    return true;
}

void process_first_pass_with_mlp_node(Config& config, StreamContext& ctx,
                                      partition::state::NodePartitionerPassState& pass_state,
                                      const Config::StreamInputConfig& input_cfg,
                                      const Config::NodePartitioningTuningConfig& tuning_cfg,
                                      FirstPassState& first_state, LongNodeID global_node_id,
                                      std::vector<LongNodeID>& neighbors) {
    const unsigned degree = static_cast<unsigned>(neighbors.size());

    if (degree >= tuning_cfg.d_max || degree == 0) {
        partition_single_node_with_timing(config, ctx, pass_state, global_node_id, neighbors);
        if (degree != 0) {
            first_state.buffer.update_neighbours_priority(neighbors);
        }
        return;
    }

    const bool added_to_buffer = first_state.buffer.addNode(global_node_id, neighbors);
    if (!added_to_buffer) {
        partition_single_node_with_timing(config, ctx, pass_state, global_node_id, neighbors);
        first_state.buffer.update_neighbours_priority(neighbors);
        return;
    }

    if (first_state.buffer.size() < tuning_cfg.max_buffer_size) {
        return;
    }

    if (first_state.always_top_extraction) {
        add_top_buffer_node_to_batch(config, pass_state, first_state);
        if (first_state.batch_state.current_batch.size() >= static_cast<size_t>(input_cfg.batch_size)) {
            flush_current_batch_mlp(config, ctx, pass_state, first_state.batch_state);
        }
        return;
    }

    first_state.buffer.loadTopNodesToBatch(first_state.batch_state.current_batch, input_cfg.batch_size,
                                           first_state.batch_state.current_batch_id);
    flush_current_batch_mlp(config, ctx, pass_state, first_state.batch_state);
}

void process_first_pass_without_mlp_node(Config& config, StreamContext& ctx,
                                         partition::state::NodePartitionerPassState& pass_state,
                                         const Config::NodePartitioningTuningConfig& tuning_cfg,
                                         Buffer& buffer, LongNodeID global_node_id,
                                         std::vector<LongNodeID>& neighbors) {
    const unsigned degree = static_cast<unsigned>(neighbors.size());

    if (degree >= tuning_cfg.d_max || degree == 0) {
        partition_single_node_with_timing(config, ctx, pass_state, global_node_id, neighbors);
        if (degree != 0) {
            buffer.update_neighbours_priority(neighbors);
        }
        return;
    }

    const bool added_to_buffer = buffer.addNode(global_node_id, neighbors);
    if (!added_to_buffer) {
        partition_single_node_with_timing(config, ctx, pass_state, global_node_id, neighbors);
        buffer.update_neighbours_priority(neighbors);
        return;
    }

    if (buffer.size() >= tuning_cfg.max_buffer_size) {
        buffer.partitionTopNode();
    }
}

void process_restream_node(Config& config, StreamContext& ctx,
                           partition::state::NodePartitionerPassState& pass_state,
                           const Config::StreamInputConfig& input_cfg,
                           RestreamPassState& restream_state, LongNodeID global_node_id,
                           std::vector<LongNodeID>& neighbors) {
    restream_state.batch_state.current_batch.emplace_back(global_node_id, neighbors);
    if (restream_state.batch_state.current_batch.size() >= static_cast<size_t>(input_cfg.batch_size) ||
        pass_state.remaining_stream_nodes == 0) {
        flush_current_batch_mlp(config, ctx, pass_state, restream_state.batch_state);
    }
}

void flush_first_pass_tail_with_mlp(Config& config, StreamContext& ctx,
                                    partition::state::NodePartitionerPassState& pass_state,
                                    const Config::StreamInputConfig& input_cfg,
                                    FirstPassState& first_state) {
    // Match legacy BuffCut semantics: tail phase always extracts with
    // loadTopNodesToBatch, independent of extraction strategy.
    while (!first_state.buffer.isEmpty()) {
        first_state.buffer.loadTopNodesToBatch(first_state.batch_state.current_batch,
                                               input_cfg.batch_size,
                                               first_state.batch_state.current_batch_id);
        if (!first_state.batch_state.current_batch.empty()) {
            config.nmbNodes = first_state.batch_state.current_batch.size();
            double local_model_time = 0.0;
            double local_partition_time = 0.0;
            double local_postprocess_time = 0.0;
            perform_mlp_on_batch(config, pass_state, first_state.batch_state.current_batch,
                                 first_state.batch_state.current_batch_id, local_model_time,
                                 local_partition_time, local_postprocess_time);
            ctx.timing.add(core::timing::StreamTimingStage::ModelBuild, local_model_time);
            ctx.timing.add(core::timing::StreamTimingStage::Partition, local_partition_time);
            ctx.timing.add(core::timing::StreamTimingStage::Postprocess, local_postprocess_time);
            first_state.batch_state.current_batch.clear();
        }
    }
}

void flush_first_pass_tail_without_mlp(Buffer& buffer) {
    while (!buffer.isEmpty()) {
        buffer.partitionTopNode();
    }
}

}  // namespace

void run_node_stream_priority_buffer_sequential(
    Config& config, const std::string& graph_filename, StreamContext& ctx,
    partition::state::NodePartitionerPassState& pass_state) {
    const auto input_cfg = config.stream_input_config();
    const auto tuning_cfg = config.node_tuning_config();

    if (tuning_cfg.max_buffer_size <= 1) {
        run_node_stream_no_priority_buffer(config, graph_filename, ctx, pass_state);
        return;
    }

    auto batch_manager = std::make_unique<BatchIDManager>(config.max_active_batches);
    config.batch_manager = batch_manager.get();

    SequentialBatchState batch_state{0, INVALID_PARTITION, {}};
    acquire_batch_window(config, batch_state.current_batch_id, batch_state.current_batch_marker);

    BufferMaintenanceTotals maintenance_totals;
    auto buffer = std::make_unique<Buffer>(config, pass_state, tuning_cfg.max_buffer_size);
    NodeStreamReader reader(config);

    const int passes = input_cfg.num_streams_passes;
    for (pass_state.restream_number = 0; pass_state.restream_number < passes;
         ++pass_state.restream_number) {
        const bool first_pass = (pass_state.restream_number == 0);
        const bool use_mlp = use_mlp_batching(config);
        const bool always_top =
            (config.batch_extraction_strategy == BATCH_EXTRACTION_STRATEGY_ALWAYS_TOP_NODE);

        ctx.timing.add(core::timing::StreamTimingStage::ReadInput,
                       begin_node_stream_pass(reader, pass_state, graph_filename));
        ensure_priority_buffer_supports_input(pass_state);

        std::vector<LongNodeID> cur_line;
        cur_line.reserve(1000);
        pass_state.total_nodes_loaded = 0;

        if (first_pass) {
            FirstPassState first_state{*buffer, batch_state, always_top};
            if (use_mlp) {
                while (pass_state.remaining_stream_nodes) {
                    LongNodeID global_node_id = 0;
                    if (!read_next_stream_node(reader, ctx, pass_state, cur_line, global_node_id)) {
                        break;
                    }
                    process_first_pass_with_mlp_node(config, ctx, pass_state, input_cfg, tuning_cfg,
                                                     first_state, global_node_id, cur_line);
                }
                flush_first_pass_tail_with_mlp(config, ctx, pass_state, input_cfg, first_state);
            } else {
                while (pass_state.remaining_stream_nodes) {
                    LongNodeID global_node_id = 0;
                    if (!read_next_stream_node(reader, ctx, pass_state, cur_line, global_node_id)) {
                        break;
                    }
                    process_first_pass_without_mlp_node(config, ctx, pass_state, tuning_cfg,
                                                        first_state.buffer, global_node_id,
                                                        cur_line);
                }
                flush_first_pass_tail_without_mlp(first_state.buffer);
            }

            collect_buffer_maintenance_totals(*buffer, maintenance_totals.add_node_time,
                                              maintenance_totals.extract_time,
                                              maintenance_totals.update_adj_time);
            buffer.reset();
        } else {
            RestreamPassState restream_state{batch_state};
            while (pass_state.remaining_stream_nodes) {
                LongNodeID global_node_id = 0;
                if (!read_next_stream_node(reader, ctx, pass_state, cur_line, global_node_id)) {
                    break;
                }
                if (should_skip_restream_node(tuning_cfg.restream_include_high_degree_nodes,
                                              tuning_cfg.d_max, cur_line)) {
                    cur_line.clear();
                    continue;
                }
                process_restream_node(config, ctx, pass_state, input_cfg, restream_state,
                                      global_node_id, cur_line);
            }
        }

        reader.end_pass();
    }

    config.batch_manager = nullptr;

    ctx.timing.add(core::timing::StreamTimingStage::BufferMaintenance,
                   maintenance_totals.add_node_time + maintenance_totals.extract_time +
                       maintenance_totals.update_adj_time);
}

}  // namespace stream_node_modes
