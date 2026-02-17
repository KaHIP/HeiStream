/****************************************************************************
 * node_partitioning_mode_prepare.cpp
 *****************************************************************************/

#include "algorithms/node_partitioning/mode/node_partitioning_mode_prepare.h"

namespace node_partitioning::mode {

void prepare_execution_config_for_mode(Config& config, NodePartitioningMode mode) {
    const auto node_tuning_cfg = config.node_tuning_config();
    if (config.bb_ratio == UNDEFINED_BB_RATIO) {
        return;
    }

    const bool priority_buffer_active = node_tuning_cfg.max_buffer_size > 1;
    const bool apply_for_mode = (mode == NodePartitioningMode::BufferedPrioritySequential) ||
                                (mode == NodePartitioningMode::BufferedPriorityParallel &&
                                 priority_buffer_active);
    if (!apply_for_mode) {
        return;
    }

    auto input_cfg = config.stream_input_config();
    input_cfg.batch_size = node_tuning_cfg.max_buffer_size / config.bb_ratio;
}

}  // namespace node_partitioning::mode
