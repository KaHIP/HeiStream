/****************************************************************************
 * stream_partition_stages.cpp
 *****************************************************************************/

#include "partition/stages/stream_partition_stages.h"

#include <algorithm>

#include "data_structure/graph_access.h"
#include "partition/balance_configuration.h"
#include "partition/graph_partitioner.h"

namespace {

void count_assigned_nodes(Config& config, const std::vector<NodeWeight>& stream_blocks_weight) {
    config.stream_assigned_nodes = 0;
    for (int i = 0; i < static_cast<int>(stream_blocks_weight.size()); i++) {
        config.stream_assigned_nodes += stream_blocks_weight[i];
    }
}

void prescribe_buffer_imbalance(Config& config, int restream_number) {
    double& global_epsilon = config.stream_global_epsilon;
    int& passes = config.num_streams_passes;
    if (restream_number && restream_number >= config.num_streams_passes - 1) {
        config.imbalance = 100 * global_epsilon;
    } else {
        double current_nodes = static_cast<double>(config.stream_assigned_nodes) + config.nmbNodes +
                               config.ghost_nodes;
        config.imbalance = 100 * config.stream_n_nodes * (1 + global_epsilon) / current_nodes - 100;
        if (passes > 1) {
            config.imbalance = std::min(config.imbalance, config.batch_inbalance);
        } else {
            config.imbalance = std::min(std::max(100 * global_epsilon, 0.75 * config.imbalance),
                                        config.batch_inbalance);
        }
    }
}

}  // namespace

namespace partition_stages {

void config_multibfs_initial_partitioning(Config& config, int curr_batch) {
    // Multibfs is only activated after at least one batch was processed.
    if (config.initial_part_multi_bfs && curr_batch >= 2) {
        config.initial_partitioning_type = INITIAL_PARTITIONING_MULTIBFS;
    }
}

void prepare_batch_partition(Config& config, int restream_number, int curr_batch,
                             const std::vector<NodeWeight>& stream_blocks_weight, graph_access& G,
                             balance_configuration& bc) {
    count_assigned_nodes(config, stream_blocks_weight);
    prescribe_buffer_imbalance(config, restream_number);

    bool already_fully_partitioned = (config.restream_vcycle && restream_number);
    bc.configurate_balance(config, G,
                           already_fully_partitioned || !config.stream_initial_bisections);
    partition_stages::config_multibfs_initial_partitioning(config, curr_batch);
}

void execute_partition(Config& config, graph_access& G, int restream_number) {
    graph_partitioner partitioner;
    partitioner.perform_partitioning(config, G, restream_number);
}

}  // namespace partition_stages
