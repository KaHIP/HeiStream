/****************************************************************************
 * edge_partitioning_reporter.cpp
 *****************************************************************************/

#include "algorithms/edge_partitioning/reporting/edge_partitioning_reporter.h"

#include <iomanip>
#include <iostream>

#include "io/output/flatbuffer_log_writer.h"
#include "io/output/partition_output_writer.h"
#include "tools/flat_buffer_writer.h"


namespace edge_partitioning::reporting {

namespace {

void print_edge_partitioning_run_configuration(const Config& config,
                                               const std::string& graph_filename) {
    const auto exec_cfg = config.stream_execution_config();
    const auto input_cfg = config.stream_input_config();
    const auto out_cfg = config.stream_output_config();
    const auto edge_tuning_cfg = config.edge_tuning_config();
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
    row("Pipeline", "edge partitioning");
    row("k", std::to_string(config.k));
    row("Seed", std::to_string(config.seed));
    row("Batch size", std::to_string(input_cfg.batch_size));
    row("Stream passes", std::to_string(input_cfg.num_streams_passes));
    row("Minimal mode", edge_tuning_cfg.minimal_mode ? "true" : "false");
    row("Dynamic alpha", edge_tuning_cfg.dynamic_alpha ? "true" : "false");
    row("Output progress", out_cfg.stream_output_progress ? "true" : "false");
    row("Write flatbuffer log", out_cfg.write_log ? "true" : "false");
    row("Evaluate", exec_cfg.evaluate ? "true" : "false");
    sep();
}

}  // namespace

void EdgePartitioningReporter::print_run_configuration(const Config& config,
                                                       const std::string& graph_filename) const {
    print_edge_partitioning_run_configuration(config, graph_filename);
}

void EdgePartitioningReporter::report(Config& config, const std::string& graph_filename,
                                      StreamContext& ctx, FlatBufferWriter& fb_writer) const {
    if (config.stream_output_progress && ctx.edge_pass_state.stream_out) {
        ctx.edge_pass_state.stream_out->close();
    }

    if (!config.benchmark && !config.stream_output_progress) {
        io::output::write_partition_output(config, *ctx.edge_pass_state.stream_nodes_assign);
    }

    const auto times = ctx.timing.snapshot();
    io::output::write_flatbuffer_log(config, fb_writer, graph_filename, times.read_input,
                                     times.model_build, times.postprocess,
                                     times.buffer_maintenance, times.partition, times.total,
                                     ctx.max_rss_kb);
}

} // namespace edge_partitioning::reporting
