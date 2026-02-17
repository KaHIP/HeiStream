/****************************************************************************
 * stream_node_mode_priority_buffer_parallel.cpp
 *****************************************************************************/

#include <atomic>
#include <memory>
#include <thread>

#include "algorithms/node_partitioning/execution/batch_id_manager.h"
#include "algorithms/node_partitioning/execution/node_execution_internal.h"
#include "algorithms/node_partitioning/reporting/node_partitioning_reporter.h"
#include "core/timing/scoped_stage_timer.h"
#include "data_structure/buffer.h"
#include "io/node_stream_reader.h"
#include "macros_assertions.h"
#include "readerwriterqueue.h"

namespace stream_node_modes {

namespace {

struct ParallelPipelineStats {
    double io_time = 0.0;
    double buffer_add_node_time = 0.0;
    double buffer_extract_time = 0.0;
    double first_phase_time = 0.0;
    double second_phase_time = 0.0;
    double part_single_node_time = 0.0;
    double mlp_time = 0.0;
    double batch_model_time = 0.0;
    double batch_partition_time = 0.0;
    double batch_postprocess_time = 0.0;

    double updating_adj_time = 0.0;
    double io_thread_time = 0.0;
    double pq_thread_time = 0.0;
    double partition_thread_time = 0.0;
    double io_wait_time = 0.0;
    double pq_wait_time = 0.0;
    double partition_wait_time = 0.0;
    double thread_io_time = 0.0;
    double thread_part_single_node_time = 0.0;
    double thread_mlp_time = 0.0;
    double thread_batch_model_time = 0.0;
    double thread_batch_partition_time = 0.0;
    double thread_batch_postprocess_time = 0.0;

    std::atomic<size_t> nodes_processed_io{0};
    std::atomic<size_t> nodes_processed_pq{0};
    std::atomic<size_t> tasks_processed_partition{0};
    std::atomic<size_t> max_input_queue_size{0};
    std::atomic<size_t> max_partition_queue_size{0};
};

template <typename EnqueueFn>
void spin_enqueue_until_success(EnqueueFn&& try_enqueue) {
    while (!try_enqueue()) {
    }
}

node_partitioning::reporting::ParallelRuntimeTableData build_parallel_runtime_report(
    const ParallelPipelineStats& stats, double total_time) {
    node_partitioning::reporting::ParallelRuntimeTableData report_data;
    report_data.total_time = total_time;
    report_data.first_phase_time = stats.first_phase_time;
    report_data.second_phase_time = stats.second_phase_time;
    report_data.io_time = stats.io_time;
    report_data.buffer_add_node_time = stats.buffer_add_node_time;
    report_data.buffer_extract_time = stats.buffer_extract_time;
    report_data.updating_adj_time = stats.updating_adj_time;
    report_data.part_single_node_time = stats.part_single_node_time;
    report_data.mlp_time = stats.mlp_time;
    report_data.io_thread_time = stats.io_thread_time;
    report_data.pq_thread_time = stats.pq_thread_time;
    report_data.partition_thread_time = stats.partition_thread_time;
    report_data.io_wait_time = stats.io_wait_time;
    report_data.pq_wait_time = stats.pq_wait_time;
    report_data.partition_wait_time = stats.partition_wait_time;
    report_data.nodes_processed_io = TIMING_LOAD(stats.nodes_processed_io);
    report_data.nodes_processed_pq = TIMING_LOAD(stats.nodes_processed_pq);
    report_data.tasks_processed_partition = TIMING_LOAD(stats.tasks_processed_partition);
    report_data.max_input_queue_size = TIMING_LOAD(stats.max_input_queue_size);
    report_data.max_partition_queue_size = TIMING_LOAD(stats.max_partition_queue_size);
    return report_data;
}

void run_io_stage(Config& config, partition::state::NodePartitionerPassState& pass_state,
                  const Config::NodePartitioningTuningConfig& node_tuning_cfg,
                  NodeStreamReader& reader,
                  moodycamel::ReaderWriterQueue<ParsedLine>& input_queue,
                  std::atomic<bool>& input_closed, ParallelPipelineStats& stats) {
    TIMED_SCOPE(stats.io_thread_time, "parallel_io_thread_total");

    std::vector<LongNodeID> adjacents;
    adjacents.reserve(1000);
    pass_state.total_nodes_loaded = 0;

    while (pass_state.remaining_stream_nodes) {
        LongNodeID global_node_id = 0;
        double io_delta = 0.0;
        if (!read_next_node_from_stream(reader, pass_state, adjacents, global_node_id, io_delta)) {
            break;
        }
        assert(global_node_id <= config.number_of_nodes);
        stats.thread_io_time += io_delta;

        if (pass_state.restream_number &&
            should_skip_restream_node(node_tuning_cfg.restream_include_high_degree_nodes,
                                      node_tuning_cfg.d_max, adjacents)) {
            adjacents.clear();
        }

        ParsedLine parsed_line{global_node_id, adjacents};
        spin_enqueue_until_success([&]() { return input_queue.try_enqueue(std::move(parsed_line)); });
        TIMING_INCREMENT(stats.nodes_processed_io);
    }

    input_closed = true;
}

void run_priority_queue_stage(Config& config, partition::state::NodePartitionerPassState& pass_state,
                              const Config::StreamInputConfig& input_cfg,
                              const Config::NodePartitioningTuningConfig& node_tuning_cfg,
                              Buffer*& buffer, bool use_mlp, bool use_priority_buffer,
                              moodycamel::ReaderWriterQueue<ParsedLine>& input_queue,
                              std::atomic<bool>& input_closed,
                              moodycamel::ReaderWriterQueue<PartitionTask>& task_queue,
                              std::atomic<bool>& task_closed,
                              ParallelPipelineStats& stats) {
    TIMED_SCOPE(stats.pq_thread_time, "parallel_pq_thread_total");

    size_t cur_batch_id = 0;
    PartitionID cur_batch_marker = INVALID_PARTITION;
    acquire_batch_window(config, cur_batch_id, cur_batch_marker);
    std::vector<BatchNode> current_batch;

    auto push_single_node_task = [&](LongNodeID node_id, std::vector<LongNodeID>&& adjacents) {
        if (pass_state.restream_number == 0 && use_priority_buffer) {
            const bool active_ghost_neighbor =
                is_active_ghost_neighbor(config, pass_state, node_id);
            buffer->update_neighbours_priority(adjacents, true, true, active_ghost_neighbor);
            (*pass_state.stream_nodes_assign)[node_id - 1] = TO_BE_PARTITIONED;
        }
        PartitionTask task(-1, {{node_id, std::move(adjacents)}});
        spin_enqueue_until_success([&]() { return task_queue.try_enqueue(std::move(task)); });
    };

    auto add_node_to_batch = [&](LongNodeID node_id, std::vector<LongNodeID>&& adjacents) {
        if (pass_state.restream_number == 0) {
            mark_node_batch_assignment(config, pass_state, node_id, cur_batch_marker);
        }
        current_batch.emplace_back(node_id, std::move(adjacents));
        if (current_batch.size() >= static_cast<size_t>(input_cfg.batch_size)) {
            PartitionTask task;
            if (build_batch_task_and_advance(config, current_batch, cur_batch_id, cur_batch_marker,
                                             task)) {
                spin_enqueue_until_success([&]() { return task_queue.try_enqueue(std::move(task)); });
            }
        }
    };

    auto emit_batch_if_any = [&]() {
        PartitionTask task;
        if (build_batch_task_and_advance(config, current_batch, cur_batch_id, cur_batch_marker,
                                         task)) {
            spin_enqueue_until_success([&]() { return task_queue.try_enqueue(std::move(task)); });
        }
    };

    auto extract_top_node_from_buffer = [&]() {
        return extract_top_buffer_node_for_batch(config, pass_state, *buffer, cur_batch_marker,
                                                 true);
    };

    ParsedLine parsed_line;
    while (true) {
        if (!input_queue.try_dequeue(parsed_line)) {
            if (input_closed) {
                break;
            }
            continue;
        }

        LongNodeID global_node_id = parsed_line.node_id;
        std::vector<LongNodeID>& adjacents = parsed_line.neighbors;
        const unsigned degree = static_cast<unsigned>(adjacents.size());

        if (pass_state.restream_number) {
            if (degree == 0) {
                push_single_node_task(global_node_id, std::vector<LongNodeID>());
                TIMING_INCREMENT(stats.nodes_processed_pq);
                continue;
            }
            add_node_to_batch(global_node_id, std::move(adjacents));
        } else {
            if (degree >= node_tuning_cfg.d_max || degree == 0) {
                push_single_node_task(global_node_id, std::move(adjacents));
                TIMING_INCREMENT(stats.nodes_processed_pq);
                continue;
            }

            if (!use_priority_buffer) {
                if (use_mlp) {
                    add_node_to_batch(global_node_id, std::move(adjacents));
                } else {
                    push_single_node_task(global_node_id, std::move(adjacents));
                }
            } else {
                const bool added_to_buffer = buffer->addNode(global_node_id, adjacents);
                if (!added_to_buffer) {
                    push_single_node_task(global_node_id, std::move(adjacents));
                }

                if (buffer->size() > node_tuning_cfg.max_buffer_size) {
                    if (use_mlp) {
                        if (config.batch_extraction_strategy !=
                            BATCH_EXTRACTION_STRATEGY_ALWAYS_TOP_NODE) {
                            buffer->loadTopNodesAndNeighborsToBatch(current_batch, input_cfg.batch_size,
                                                                    cur_batch_id);
                            emit_batch_if_any();
                        } else {
                            auto [node_id, node_adjacents] = extract_top_node_from_buffer();
                            add_node_to_batch(node_id, std::move(node_adjacents));
                        }
                    } else {
                        auto [node_id, node_adjacents] = extract_top_node_from_buffer();
                        push_single_node_task(node_id, std::move(node_adjacents));
                    }
                }
            }
        }

        TIMING_INCREMENT(stats.nodes_processed_pq);
    }

    if (pass_state.restream_number == 0 && use_priority_buffer) {
        if (use_mlp) {
            while (!buffer->isEmpty()) {
                if (config.batch_extraction_strategy != BATCH_EXTRACTION_STRATEGY_ALWAYS_TOP_NODE) {
                    buffer->loadTopNodesAndNeighborsToBatch(current_batch, input_cfg.batch_size,
                                                            cur_batch_id);
                    emit_batch_if_any();
                } else {
                    auto [node_id, node_adjacents] = extract_top_node_from_buffer();
                    add_node_to_batch(node_id, std::move(node_adjacents));
                }
            }
        } else {
            while (!buffer->isEmpty()) {
                auto [node_id, node_adjacents] = extract_top_node_from_buffer();
                push_single_node_task(node_id, std::move(node_adjacents));
            }
        }
    }

    emit_batch_if_any();
    task_closed = true;
}

void run_partition_stage(Config& config, partition::state::NodePartitionerPassState& pass_state,
                         moodycamel::ReaderWriterQueue<PartitionTask>& task_queue,
                         std::atomic<bool>& task_closed,
                         ParallelPipelineStats& stats) {
    TIMED_SCOPE(stats.partition_thread_time, "parallel_partition_thread_total");

    PartitionTask task;
    while (true) {
        if (!task_queue.try_dequeue(task)) {
            if (task_closed) {
                break;
            }
            continue;
        }

        if (task.batch_id == -1) {
            LongNodeID node_id = task.nodes[0].first;
            auto& adjacents = task.nodes[0].second;

            if (!pass_state.restream_number) {
                TIMED_SCOPE(stats.thread_part_single_node_time, "parallel_single_node_partition");
                assign_single_node_direct(config, pass_state, node_id, adjacents);
            }
        } else {
            std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>> batch_nodes(
                std::move(task.nodes));
            config.nmbNodes = batch_nodes.size();

            double local_model_time = 0.0;
            double local_partition_time = 0.0;
            double local_postprocess_time = 0.0;
            perform_mlp_on_batch(config, pass_state, batch_nodes, task.batch_id, local_model_time,
                                 local_partition_time, local_postprocess_time);
            config.batch_manager->release_id(task.batch_id);

            stats.thread_batch_model_time += local_model_time;
            stats.thread_batch_partition_time += local_partition_time;
            stats.thread_batch_postprocess_time += local_postprocess_time;
            stats.thread_mlp_time +=
                local_model_time + local_partition_time + local_postprocess_time;
            task.nodes.clear();
        }

        TIMING_INCREMENT(stats.tasks_processed_partition);
    }
}

void run_parallel_priority_buffer_pass(Config& config, const std::string& graph_filename,
                                       partition::state::NodePartitionerPassState& pass_state,
                                       Buffer*& buffer, ParallelPipelineStats& stats) {
    auto input_cfg = config.stream_input_config();
    const auto node_tuning_cfg = config.node_tuning_config();
    NodeStreamReader reader(config);

    stats.io_time += begin_node_stream_pass(reader, pass_state, graph_filename);
    ensure_priority_buffer_supports_input(pass_state);

    const bool use_mlp = use_mlp_batching(config);
    const bool use_priority_buffer = (node_tuning_cfg.max_buffer_size > 1);
    TIMED_SCOPE(stats.first_phase_time, "parallel_pass_first_phase");

    // Unified queue backends keep stage code independent from queue implementation.
    moodycamel::ReaderWriterQueue<ParsedLine> input_queue(node_tuning_cfg.max_input_q_size);
    moodycamel::ReaderWriterQueue<PartitionTask> task_queue(100);
    std::atomic<bool> input_closed{false};
    std::atomic<bool> task_closed{false};

    auto push_to_partition_queue = [&](PartitionTask&& task) {
        spin_enqueue_until_success([&]() { return task_queue.try_enqueue(std::move(task)); });
    };
    if (buffer != nullptr) {
        buffer->set_push_task_callback(push_to_partition_queue);
    }

    std::thread io_reader([&]() {
        run_io_stage(config, pass_state, node_tuning_cfg, reader, input_queue, input_closed, stats);
    });
    std::thread pq_handler([&]() {
        run_priority_queue_stage(config, pass_state, input_cfg, node_tuning_cfg, buffer, use_mlp,
                                 use_priority_buffer, input_queue, input_closed, task_queue,
                                 task_closed, stats);
    });
    std::thread partition_worker([&]() {
        run_partition_stage(config, pass_state, task_queue, task_closed, stats);
    });

    io_reader.join();
    pq_handler.join();
    partition_worker.join();

    if (pass_state.restream_number == 0 && buffer != nullptr) {
        collect_buffer_maintenance_totals(*buffer, stats.buffer_add_node_time,
                                          stats.buffer_extract_time, stats.updating_adj_time);
        delete buffer;
        buffer = nullptr;
    }

    pass_state.ghostkey_to_edges.reset();
    pass_state.ghostkey_to_node.reset();
    pass_state.ghostglobal_to_ghostkey.reset();
    reader.end_pass();
}

}  // namespace

void run_node_stream_priority_buffer_parallel(
    Config& config, const std::string& graph_filename, StreamContext& ctx,
    partition::state::NodePartitionerPassState& pass_state) {
    double total_time = 0.0;
    auto input_cfg = config.stream_input_config();
    const auto node_tuning_cfg = config.node_tuning_config();
    ParallelPipelineStats stats;

    auto batch_manager = std::make_unique<BatchIDManager>(config.max_active_batches);
    config.batch_manager = batch_manager.get();

    {
        TIMED_SCOPE(total_time, "parallel_priority_buffer_total");
        Buffer* buffer = nullptr;
        if (node_tuning_cfg.max_buffer_size > 1) {
            buffer = new Buffer(config, pass_state, node_tuning_cfg.max_buffer_size);
        }
        int& passes = input_cfg.num_streams_passes;
        for (pass_state.restream_number = 0; pass_state.restream_number < passes;
             ++pass_state.restream_number) {
            run_parallel_priority_buffer_pass(config, graph_filename, pass_state, buffer, stats);
        }
    }

    config.batch_manager = nullptr;

    stats.io_time += stats.thread_io_time;
    stats.part_single_node_time += stats.thread_part_single_node_time;
    stats.mlp_time += stats.thread_mlp_time;
    stats.batch_model_time += stats.thread_batch_model_time;
    stats.batch_partition_time += stats.thread_batch_partition_time;
    stats.batch_postprocess_time += stats.thread_batch_postprocess_time;

    if (config.print_times) {
        const auto report_data = build_parallel_runtime_report(stats, total_time);
        node_partitioning::reporting::print_parallel_runtime_table(report_data);
    }

    for (LongNodeID i = 0; i < config.number_of_nodes; i++) {
        ASSERT_TRUE((*pass_state.stream_nodes_assign)[i] < TO_BE_PARTITIONED - 10000);
    }

    ctx.timing.set(core::timing::StreamTimingStage::ReadInput, stats.io_time);
    ctx.timing.set(core::timing::StreamTimingStage::ModelBuild, stats.batch_model_time);
    ctx.timing.set(core::timing::StreamTimingStage::Postprocess, stats.batch_postprocess_time);
    ctx.timing.set(core::timing::StreamTimingStage::Partition,
                   stats.part_single_node_time + stats.batch_partition_time);
    ctx.timing.set(core::timing::StreamTimingStage::BufferMaintenance,
                   stats.buffer_add_node_time + stats.buffer_extract_time +
                       stats.updating_adj_time);
}

}  // namespace stream_node_modes
