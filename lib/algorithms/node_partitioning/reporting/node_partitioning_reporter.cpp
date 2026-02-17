/****************************************************************************
 * node_partitioning_reporter.cpp
 *****************************************************************************/

#include "algorithms/node_partitioning/reporting/node_partitioning_reporter.h"

#include <iomanip>
#include <iostream>

#include "io/output/flatbuffer_log_writer.h"
#include "io/output/partition_output_writer.h"
#include "tools/flat_buffer_writer.h"


namespace node_partitioning::reporting {

namespace {

void print_node_partitioning_run_configuration(const Config& config,
                                               const std::string& graph_filename) {
    const auto exec_cfg = config.stream_execution_config();
    const auto input_cfg = config.stream_input_config();
    const auto out_cfg = config.stream_output_config();
    const auto node_tuning_cfg = config.node_tuning_config();
    if (!exec_cfg.evaluate) {
        return;
    }

    const int key_width = 26;
    const int val_width = 54;
    auto sep = [&]() {
        std::cout << "+" << std::string(key_width + 2, '-') << "+"
                  << std::string(val_width + 2, '-') << "+\n";
    };
    auto row = [&](const std::string& key, const std::string& value) {
        std::cout << "| " << std::left << std::setw(key_width) << key << " | " << std::right
                  << std::setw(val_width) << value << " |\n";
    };

    std::cout << "\n================ Stream Run Configuration ================\n";
    sep();
    row("Graph file", graph_filename);
    row("Pipeline", "node partitioning");
    row("k", std::to_string(config.k));
    row("Seed", std::to_string(config.seed));
    row("Batch size", std::to_string(input_cfg.batch_size));
    row("Buffer size", std::to_string(node_tuning_cfg.max_buffer_size));
    row("Stream passes", std::to_string(input_cfg.num_streams_passes));
    row("Run parallel", exec_cfg.run_parallel ? "true" : "false");
    row("Output progress", out_cfg.stream_output_progress ? "true" : "false");
    row("Write flatbuffer log", out_cfg.write_log ? "true" : "false");
    row("Evaluate", exec_cfg.evaluate ? "true" : "false");
    sep();
}

}  // namespace

void print_parallel_runtime_table(const ParallelRuntimeTableData& data) {
#ifdef ENABLE_TIME_MEASUREMENTS
    const double total_time = data.total_time;
    const double sum_detailed =
        data.io_time + data.buffer_add_node_time + data.buffer_extract_time +
        data.updating_adj_time + data.mlp_time + data.part_single_node_time;
    const double total_thread_time =
        data.io_thread_time + data.pq_thread_time + data.partition_thread_time;

    std::cout << "┌─────────────────────────┬───────────────┬───────────────┐" << '\n';
    std::cout << "│ Metric                  │ Time (s)      │ Percentage    │" << '\n';
    std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << '\n';
    std::cout << "│ Total time              │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << total_time << " │ " << std::setw(13) << "100%" << " │"
              << '\n';
    std::cout << "│ First phase time        │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.first_phase_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0) << (data.first_phase_time / total_time * 100)
              << "%" << " │" << '\n';
    std::cout << "│ Second phase time       │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.second_phase_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0) << (data.second_phase_time / total_time * 100)
              << "%" << " │" << '\n';
    std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << '\n';
    std::cout << "│ IO time                 │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.io_time << " │ " << std::setw(12) << std::fixed
              << std::setprecision(0) << (data.io_time / total_time * 100) << "%" << " │" << '\n';
    std::cout << "│ Buffer add node time    │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.buffer_add_node_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0)
              << (data.buffer_add_node_time / total_time * 100) << "%" << " │" << '\n';
    std::cout << "│ Buffer extract time     │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.buffer_extract_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0)
              << (data.buffer_extract_time / total_time * 100) << "%"
              << " │" << '\n';
    std::cout << "│ Updating adj time       │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.updating_adj_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0) << (data.updating_adj_time / total_time * 100)
              << "%" << " │" << '\n';
    std::cout << "│ Part single node time   │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.part_single_node_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0)
              << (data.part_single_node_time / total_time * 100) << "%" << " │" << '\n';
    std::cout << "│ MLP time                │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.mlp_time << " │ " << std::setw(12) << std::fixed
              << std::setprecision(0) << (data.mlp_time / total_time * 100) << "%" << " │" << '\n';
    std::cout << "│ Sum of detailed times   │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << sum_detailed << " │ " << std::setw(12) << std::fixed
              << std::setprecision(0) << (sum_detailed / total_time * 100) << "%" << " │" << '\n';
    std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << '\n';
    std::cout << "│ THREAD RUNTIME TRACKING │               │               │" << '\n';
    std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << '\n';
    std::cout << "│ IOReader thread time    │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.io_thread_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0) << (data.io_thread_time / total_time * 100)
              << "%" << " │" << '\n';
    std::cout << "│ PQHandler thread time   │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.pq_thread_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0) << (data.pq_thread_time / total_time * 100)
              << "%" << " │" << '\n';
    std::cout << "│ PartitionWorker time    │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.partition_thread_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0)
              << (data.partition_thread_time / total_time * 100) << "%" << " │" << '\n';
    std::cout << "│ Total thread time       │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << total_thread_time << " │ " << std::setw(12) << std::fixed
              << std::setprecision(0) << (total_thread_time / total_time * 100) << "%" << " │"
              << '\n';
    std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << '\n';
    std::cout << "│ BOTTLENECK ANALYSIS     │               │               │" << '\n';
    std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << '\n';
    std::cout << "│ IO wait time            │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.io_wait_time << " │ " << std::setw(12) << std::fixed
              << std::setprecision(0) << (data.io_wait_time / total_time * 100) << "%" << " │"
              << '\n';
    std::cout << "│ PQ wait time            │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.pq_wait_time << " │ " << std::setw(12) << std::fixed
              << std::setprecision(0) << (data.pq_wait_time / total_time * 100) << "%" << " │"
              << '\n';
    std::cout << "│ Partition wait time     │ " << std::setw(13) << std::fixed
              << std::setprecision(3) << data.partition_wait_time << " │ " << std::setw(12)
              << std::fixed << std::setprecision(0)
              << (data.partition_wait_time / total_time * 100) << "%" << " │" << '\n';
    std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << '\n';
    std::cout << "│ Nodes/sec IO            │ " << std::setw(13) << std::fixed
              << std::setprecision(0) << (data.nodes_processed_io / total_time)
              << " │               │" << '\n';
    std::cout << "│ Nodes/sec PQ            │ " << std::setw(13) << std::fixed
              << std::setprecision(0) << (data.nodes_processed_pq / total_time)
              << " │               │" << '\n';
    std::cout << "│ Tasks/sec Partition     │ " << std::setw(13) << std::fixed
              << std::setprecision(0) << (data.tasks_processed_partition / total_time)
              << " │               │" << '\n';
    std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << '\n';
    std::cout << "│ Max input queue size    │ " << std::setw(13) << std::fixed
              << std::setprecision(0) << data.max_input_queue_size << " │               │" << '\n';
    std::cout << "│ Max partition queue size│ " << std::setw(13) << std::fixed
              << std::setprecision(0) << data.max_partition_queue_size
              << " │               │" << '\n';
    std::cout << "└─────────────────────────┴───────────────┴───────────────┘" << '\n';
#else
    std::cout << "Timing disabled - compile with -DENABLE_TIME_MEASUREMENTS to see detailed timing "
                 "information"
              << '\n';
    (void)data;
#endif
}

void NodePartitioningReporter::print_run_configuration(const Config& config,
                                                       const std::string& graph_filename) const {
    print_node_partitioning_run_configuration(config, graph_filename);
}

void NodePartitioningReporter::report(Config& config, const std::string& graph_filename,
                                      StreamContext& ctx, FlatBufferWriter& fb_writer) const {
    if (!config.suppress_output) {
        io::output::write_partition_output(config, *ctx.node_pass_state.stream_nodes_assign);
    } else {
        std::cout << "No partition will be written as output." << '\n';
    }

    const auto times = ctx.timing.snapshot();
    io::output::write_flatbuffer_log(config, fb_writer, graph_filename, times.read_input,
                                     times.model_build, times.postprocess,
                                     times.buffer_maintenance, times.partition, times.total,
                                     ctx.max_rss_kb);
}

} // namespace node_partitioning::reporting
