/****************************************************************************
 * edge_partitioning_setup.cpp
 *****************************************************************************/

#include "algorithms/edge_partitioning/setup/edge_partitioning_setup.h"

#include <cmath>
#include <fstream>
#include <memory>
#include <limits>
#include <sparsehash/dense_hash_set>

#include "partition/setup/partition_stream_common_setup.h"


namespace edge_partitioning::setup {

void EdgePartitioningSetup::prepare_run(Config& config, const std::string& graph_filename) const {
    (void)graph_filename;
    const auto edge_tuning_cfg = config.edge_tuning_config();
    (void)edge_tuning_cfg;
    // Edge pipeline run defaults are centralized in setup.
    config.edge_partition = true;
    config.global_cycle_iterations = config.reps;
    config.use_queue = true;
    config.initial_partitioning_repetitions = 4;
}

void EdgePartitioningSetup::initialize_pass(Config& partition_config,
                                            partition::state::EdgePartitionerPassState& pass_state,
                                            const io::stream::StreamHeader& header) const {
    const auto part_cfg = partition_config.partition_algorithm_config();
    const auto out_cfg = partition_config.stream_output_config();
    const auto exec_cfg = partition_config.stream_execution_config();
    const auto edge_tuning_cfg = partition_config.edge_tuning_config();
    partition::setup::stream_pass_setup::ensure_output_filename(partition_config);
    if (out_cfg.stream_output_progress) {
        pass_state.stream_out = std::make_unique<std::ofstream>(out_cfg.filename_output.c_str());
    } else {
        pass_state.stream_out.reset();
    }

    partition_config.total_nodes = header.nodes;
    partition_config.total_edges = header.edges;

    pass_state.reset_for_pass(header.nodes, header.edges, header.binary_input);

    // Edge streaming partitions edge-nodes (split graph).
    partition_config.total_stream_edges = header.edges;

    if (!pass_state.stream_nodes_assign &&
        ((!exec_cfg.benchmark && !out_cfg.stream_output_progress) ||
         partition_config.num_streams_passes > 1)) {
        pass_state.stream_nodes_assign =
            std::make_unique<std::vector<PartitionID>>(header.edges, INVALID_PARTITION);
    }

    if (!pass_state.blocks_on_node && !edge_tuning_cfg.minimal_mode) {
        auto blocks = std::make_unique<std::vector<google::dense_hash_set<PartitionID>>>(
            pass_state.remaining_stream_nodes_og);
        for (auto& set : *blocks) {
            set.set_empty_key(INVALID_PARTITION);
        }
        pass_state.blocks_on_node = std::move(blocks);
    }
    if (edge_tuning_cfg.minimal_mode) {
        pass_state.blocks_on_node.reset();
    }

    if (!pass_state.blocks_on_node_minimal && edge_tuning_cfg.minimal_mode) {
        constexpr NodeID kUnsetBlock = std::numeric_limits<NodeID>::max();
        pass_state.blocks_on_node_minimal =
            std::make_unique<std::vector<NodeID>>(pass_state.remaining_stream_nodes_og,
                                                  kUnsetBlock);
    }
    if (!edge_tuning_cfg.minimal_mode) {
        pass_state.blocks_on_node_minimal.reset();
    }

    partition::setup::stream_pass_setup::init_common_stream_state(
        partition_config, pass_state.restream_number, pass_state.remaining_stream_partition_nodes,
        pass_state.remaining_stream_partition_nodes, pass_state.stream_blocks_weight,
        pass_state.add_blocks_weight);

    // Keep historical alpha semantics used by original edge pipeline.
    if (edge_tuning_cfg.dynamic_alpha) {
        partition_config.fennel_edges = 6 * header.edges;
    } else {
        partition_config.fennel_edges = 5 * header.edges;
    }
    if (edge_tuning_cfg.num_split_edges == std::numeric_limits<NodeID>::max()) {
        partition_config.fennel_alpha =
            partition_config.fennel_edges *
            std::pow(part_cfg.k, partition_config.fennel_gamma - 1) /
            (std::pow(pass_state.remaining_stream_partition_nodes, partition_config.fennel_gamma));
    } else {
        partition_config.fennel_alpha =
            header.edges * std::pow(part_cfg.k, partition_config.fennel_gamma - 1) /
            (std::pow(pass_state.remaining_stream_nodes_og, partition_config.fennel_gamma));
    }

    partition_config.fennel_alpha_gamma =
        partition_config.fennel_alpha * partition_config.fennel_gamma;
}

} // namespace edge_partitioning::setup
