/****************************************************************************
 * stream_node_mode_common.cpp
 *****************************************************************************/

#include <cmath>
#include <iomanip>
#include <memory>

#include "algorithms/node_partitioning/execution/node_execution_internal.h"
#include "algorithms/node_partitioning/model/node_batch_model_builder.h"
#include "algorithms/node_partitioning/setup/node_partitioning_setup.h"
#include "core/timing/scoped_stage_timer.h"
#include "data_structure/buffer.h"
#include "io/node_stream_reader.h"
#include "io/output/partition_output_writer.h"
#include "partition/balance_configuration.h"
#include "partition/stages/batch_partition_execution.h"
#include "tools/stream_util.h"

namespace stream_node_modes {

bool use_mlp_batching(const Config& config) {
    // Batch model partitioning is enabled whenever batch size exceeds 1.
    const auto input_cfg = config.stream_input_config();
    return input_cfg.batch_size != 1;
}

double begin_node_stream_pass(NodeStreamReader& reader,
                              partition::state::NodePartitionerPassState& pass_state,
                              const std::string& graph_filename) {
    const double io_before = reader.disk_read_time();
    reader.begin_pass(graph_filename);
    node_partitioning::setup::NodePartitioningSetup setup;
    setup.initialize_pass(reader.config_ref(), pass_state, reader.header());
    return reader.disk_read_time() - io_before;
}

bool read_next_node_from_stream(NodeStreamReader& reader,
                                partition::state::NodePartitionerPassState& pass_state,
                                std::vector<LongNodeID>& neighbors, LongNodeID& global_node_id,
                                double& io_delta_seconds) {
    const double io_before = reader.disk_read_time();
    if (!reader.next_node_line(neighbors)) {
        io_delta_seconds = reader.disk_read_time() - io_before;
        return false;
    }
    io_delta_seconds = reader.disk_read_time() - io_before;
    --pass_state.remaining_stream_nodes;
    global_node_id = ++pass_state.total_nodes_loaded;
    return true;
}

void ensure_priority_buffer_supports_input(
    const partition::state::NodePartitionerPassState& pass_state) {
    if (pass_state.remaining_stream_ew != 0) {
        std::cerr << "Priority buffer mode does not yet support weighted node streams. "
                     "Please run with --buffer_size=0."
                  << '\n';
        exit(1);
    }
}

void write_stream_progress_if_enabled(
    Config& config, const partition::state::NodePartitionerPassState& pass_state) {
    if (config.stream_output_progress) {
        io::output::write_partition_output(config, *pass_state.stream_nodes_assign);
    }
}

bool is_active_ghost_neighbor(const Config& config,
                              const partition::state::NodePartitionerPassState& pass_state,
                              LongNodeID node_id) {
    if (!config.ghost_neighbors_enabled) {
        return false;
    }
    const PartitionID part = (*pass_state.stream_nodes_assign)[node_id - 1];
    return config.k <= part && part < 2 * config.k;
}

bool should_skip_restream_node(bool include_high_degree_nodes, LongNodeID d_max,
                               const std::vector<LongNodeID>& adjacents) {
    return !include_high_degree_nodes &&
           (adjacents.empty() || adjacents.size() >= static_cast<size_t>(d_max));
}

void mark_node_batch_assignment(const Config& config,
                                partition::state::NodePartitionerPassState& pass_state,
                                LongNodeID node_id, PartitionID batch_marker) {
    if (config.sep_batch_marker) {
        (*pass_state.stream_nodes_batch_marker)[node_id - 1] = batch_marker;
    } else {
        (*pass_state.stream_nodes_assign)[node_id - 1] = batch_marker;
    }
}

BatchNode extract_top_buffer_node_for_batch(
    Config& config, partition::state::NodePartitionerPassState& pass_state, Buffer& buffer,
    PartitionID batch_marker, bool use_parallel_update) {
    LongNodeID node_id = buffer.deleteMax();
    std::vector<LongNodeID> node_adjacents = std::move(buffer.get_adjacents(node_id));
    const bool active_ghost_neighbor = is_active_ghost_neighbor(config, pass_state, node_id);

    if (use_parallel_update) {
        buffer.update_neighbours_priority(node_adjacents, false, true, active_ghost_neighbor);
    } else {
        buffer.update_neighbours_priority(node_adjacents, false);
    }

    mark_node_batch_assignment(config, pass_state, node_id, batch_marker);
    buffer.completely_remove_node(node_id);
    return std::make_pair(node_id, std::move(node_adjacents));
}

void advance_batch_window(Config& config, size_t& current_batch_id,
                          PartitionID& current_batch_marker) {
    config.batch_manager->release_id(current_batch_id);
    current_batch_id = config.batch_manager->acquire_id();
    current_batch_marker = config.batch_manager->get_batch_marker(current_batch_id);
}

void acquire_batch_window(Config& config, size_t& current_batch_id,
                          PartitionID& current_batch_marker) {
    current_batch_id = config.batch_manager->acquire_id();
    current_batch_marker = config.batch_manager->get_batch_marker(current_batch_id);
}

bool build_batch_task_and_advance(Config& config, std::vector<BatchNode>& current_batch,
                                  size_t& current_batch_id, PartitionID& current_batch_marker,
                                  PartitionTask& task) {
    if (current_batch.empty()) {
        return false;
    }
    task = PartitionTask(static_cast<int>(current_batch_id), std::move(current_batch));
    current_batch.clear();
    advance_batch_window(config, current_batch_id, current_batch_marker);
    return true;
}

void collect_buffer_maintenance_totals(Buffer& buffer, double& add_node_time, double& extract_time,
                                       double& update_adj_time) {
    add_node_time += buffer.get_insert_time();
    extract_time += buffer.get_extract_time();
    update_adj_time += buffer.get_update_adj_time();
}

void assign_single_node_direct(Config& config,
                               partition::state::NodePartitionerPassState& pass_state,
                               LongNodeID global_node_id, std::vector<LongNodeID>& adjacents) {
    // Baseline single-node path: isolated nodes use dedicated handling.
    if (adjacents.empty()) {
        assign_isolated_node(config, pass_state, global_node_id);
    } else {
        partition_single_node(config, pass_state, global_node_id, adjacents);
    }
}

void extract_fennel_neighbors(int remaining_stream_ew, const std::vector<LongNodeID>& line_numbers,
                              LongNodeID current_global_node, NodeWeight& node_weight,
                              std::vector<LongNodeID>& adjacents, LongEdgeID& used_edges) {
    // Parse one serialized line and keep only already-seen neighbors for
    // one-pass Fennel assignment (edges to future nodes are ignored here).
    bool read_ew = false;
    bool read_nw = false;
    switch (remaining_stream_ew) {
        case 1:
            read_ew = true;
            break;
        case 10:
            read_nw = true;
            break;
        case 11:
            read_ew = true;
            read_nw = true;
            break;
    }

    LongNodeID col_counter = 0;
    node_weight = 1;
    if (read_nw && col_counter < static_cast<LongNodeID>(line_numbers.size())) {
        node_weight = static_cast<NodeWeight>(line_numbers[col_counter++]);
    }

    while (col_counter < static_cast<LongNodeID>(line_numbers.size())) {
        LongNodeID target = line_numbers[col_counter++];
        if (read_ew && col_counter < static_cast<LongNodeID>(line_numbers.size())) {
            ++col_counter;  // Skip edge weight.
        }
        if (target < current_global_node) {
            adjacents.push_back(target);
            ++used_edges;
        }
    }
}

void perform_mlp_on_batch(Config& partition_config,
                          partition::state::NodePartitionerPassState& pass_state,
                          std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>& batch_nodes,
                          size_t batch_id, double& model_time, double& partition_time,
                          double& postprocess_time) {
    // Build a temporary local model graph for this batch and immediately
    // project partition assignments back to global stream state.
    graph_access G;
    balance_configuration bc;
    G.set_partition_count(partition_config.k);
    std::vector<NodeID> local_to_global(partition_config.nmbNodes);
    partition_config.local_to_global_map = &local_to_global;
    std::vector<std::vector<LongNodeID>> batch_unpartitioned;
    if (partition_config.ghost_neighbors_enabled) {
        batch_unpartitioned.resize(static_cast<size_t>(partition_config.nmbNodes));
        partition_config.batch_unpartitioned_neighbors = &batch_unpartitioned;
    }

    {
        TIMED_SCOPE(model_time, "node_batch_model_build");
        node_partitioning::model::build_batch_model_from_reordered_nodes(
            partition_config, pass_state, G, batch_nodes, batch_id);
    }

    {
        TIMED_SCOPE(partition_time, "node_batch_partition");
        perform_node_batch_partition(partition_config, pass_state, G, bc);
    }
    {
        TIMED_SCOPE(postprocess_time, "node_batch_postprocess");
        postprocess_node_batch_partition(partition_config, pass_state, G);
    }

    partition_config.local_to_global_map = nullptr;
    partition_config.batch_unpartitioned_neighbors = nullptr;
}

}  // namespace stream_node_modes
