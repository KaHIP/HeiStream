/****************************************************************************
 * node_partitioning_setup.cpp
 *****************************************************************************/

#include "algorithms/node_partitioning/setup/node_partitioning_setup.h"

#include <cmath>
#include <memory>

#include "partition/assignment/node_fennel_assign.h"
#include "partition/setup/partition_stream_common_setup.h"


namespace node_partitioning::setup {

void NodePartitioningSetup::prepare_run(Config& config, const std::string& graph_filename) const {
    (void)graph_filename;
    const auto exec_cfg = config.stream_execution_config();
    (void)exec_cfg;
    config.edge_partition = false;
}

void NodePartitioningSetup::initialize_pass(Config& partition_config,
                                            partition::state::NodePartitionerPassState& pass_state,
                                            const io::stream::StreamHeader& header) const {
    const auto part_cfg = partition_config.partition_algorithm_config();
    pass_state.reset_from_header(header.nodes, header.edges, header.edge_weight_flag,
                                 header.binary_input);

    // Node stream keeps assignments by original global node IDs.
    partition_config.total_nodes = header.nodes;
    partition_config.total_edges = header.edges;
    partition_config.number_of_nodes = header.nodes;

    if (!pass_state.stream_nodes_assign) {
        pass_state.stream_nodes_assign =
            std::make_unique<std::vector<PartitionID>>(header.nodes, INVALID_PARTITION);
    }

    if (partition_config.sep_batch_marker && !pass_state.stream_nodes_batch_marker) {
        pass_state.stream_nodes_batch_marker =
            std::make_unique<std::vector<PartitionID>>(header.nodes, INVALID_PARTITION);
    }
    if (!partition_config.sep_batch_marker) {
        pass_state.stream_nodes_batch_marker.reset();
    }
    partition_config.total_stream_edges = header.nodes;

    // In node streaming, balancing may use node-only or node+edge mass.
    auto total_weight =
        (partition_config.balance_edges) ? (header.nodes + (2 * header.edges)) : header.nodes;
    partition::setup::stream_pass_setup::init_common_stream_state(
        partition_config, pass_state.restream_number, pass_state.remaining_stream_nodes,
        total_weight, pass_state.stream_blocks_weight, pass_state.add_blocks_weight);
    partition_config.max_block_weight =
        partition_config.stream_total_upperbound;

    partition_config.fennel_alpha = header.edges *
                                    std::pow(part_cfg.k, partition_config.fennel_gamma - 1) /
                                    (std::pow(header.nodes, partition_config.fennel_gamma));

    partition_config.fennel_alpha_gamma =
        partition_config.fennel_alpha * partition_config.fennel_gamma;

    // Keep node-side queue mode in sync with current pass block weights.
    node_fennel_assign::rebuild_block_min_queue(partition_config, pass_state);

    // Ghost-neighbor approximation is only used during the first pass.
    partition_config.ghost_neighbors_enabled =
        partition_config.ghost_neighbors_enabled && pass_state.restream_number == 0;
    if (partition_config.ghost_neighbors_enabled) {
        partition_config.ghost_weight = 1.0f / partition_config.default_weight_non_ghost;
        partition_config.inv_ghost_weight = 1.0f - partition_config.ghost_weight;
    }
}

} // namespace node_partitioning::setup
