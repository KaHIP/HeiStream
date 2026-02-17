/****************************************************************************
 * partition_stream_common_setup.h
 *****************************************************************************/

#ifndef PARTITION_STREAM_COMMON_SETUP_H_
#define PARTITION_STREAM_COMMON_SETUP_H_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "partition/partition_config.h"



namespace partition::setup::stream_pass_setup {

// Ensure output filename is always initialized.
inline void ensure_output_filename(Config& partition_config) {
    if (!partition_config.filename_output.compare("")) {
        partition_config.filename_output = "tmp_output.txt";
    }
}

// Lazily initialize block-weight vectors used by stream assignment and balancing.
inline void ensure_block_weights(
    Config& partition_config, std::unique_ptr<std::vector<NodeWeight>>& stream_blocks_weight_owner,
    std::unique_ptr<std::vector<NodeWeight>>& add_blocks_weight_owner) {
    (void)partition_config;
    if (!stream_blocks_weight_owner) {
        stream_blocks_weight_owner =
            std::make_unique<std::vector<NodeWeight>>(partition_config.k, 0);
    }
    if (!add_blocks_weight_owner) {
        add_blocks_weight_owner = std::make_unique<std::vector<NodeWeight>>(partition_config.k, 0);
    }
}

// Reset per-pass stream counters before processing batches.
inline void reset_stream_counters(Config& partition_config, LongNodeID stream_nodes) {
    partition_config.stream_assigned_nodes = 0;
    partition_config.stream_n_nodes = stream_nodes;
}

// Compute per-block upper bound used by stream assignment heuristics.
inline void compute_stream_total_upperbound(Config& partition_config, int restream_number,
                                            double total_weight) {
    if (partition_config.num_streams_passes > 1 + restream_number) {
        partition_config.stream_total_upperbound =
            std::ceil(((100 + 1.5 * partition_config.imbalance) / 100.) *
                      (total_weight / static_cast<double>(partition_config.k)));
    } else {
        partition_config.stream_total_upperbound =
            std::ceil(((100 + partition_config.imbalance) / 100.) *
                      (total_weight / static_cast<double>(partition_config.k)));
    }
}

// Derive batch size and number of batches for current stream pass.
inline void init_stream_batching(Config& partition_config, LongNodeID remaining_nodes,
                                 LongNodeID& batch_len) {
    if (batch_len == 0) {
        batch_len = static_cast<LongNodeID>(
            std::ceil(remaining_nodes / static_cast<double>(partition_config.k)));
    }
    partition_config.nmbNodes = std::min(batch_len, remaining_nodes);
    partition_config.n_batches =
        std::ceil(remaining_nodes / static_cast<double>(partition_config.nmbNodes));
}

inline void init_stream_batching(Config& partition_config, LongNodeID remaining_nodes) {
    init_stream_batching(partition_config, remaining_nodes, partition_config.batch_size);
}

// Shared post-header initialization for partition stream readers.
inline void init_common_stream_state(
    Config& partition_config, int restream_number, LongNodeID remaining_nodes, double total_weight,
    std::unique_ptr<std::vector<NodeWeight>>& stream_blocks_weight_owner,
    std::unique_ptr<std::vector<NodeWeight>>& add_blocks_weight_owner) {
    ensure_output_filename(partition_config);
    ensure_block_weights(partition_config, stream_blocks_weight_owner, add_blocks_weight_owner);
    reset_stream_counters(partition_config, remaining_nodes);
    compute_stream_total_upperbound(partition_config, restream_number, total_weight);
    partition_config.quotient_nodes = partition_config.k;
    init_stream_batching(partition_config, remaining_nodes);
}

} // namespace partition::setup::stream_pass_setup



#endif /* PARTITION_STREAM_COMMON_SETUP_H_ */
