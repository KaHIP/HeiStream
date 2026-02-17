/****************************************************************************
 * node_partitioning_mode_resolver.cpp
 *****************************************************************************/

#include "algorithms/node_partitioning/mode/node_partitioning_mode_resolver.h"


namespace node_partitioning::mode {

NodePartitioningMode NodePartitioningModeResolver::resolve(const Config& config) const {
    const auto exec_cfg = config.stream_execution_config();
    const auto input_cfg = config.stream_input_config();
    const auto node_tuning_cfg = config.node_tuning_config();
    const bool use_parallel_mode = exec_cfg.run_parallel;
    const bool use_priority_buffer = (node_tuning_cfg.max_buffer_size > 1);
    const bool direct_fennel =
        (!use_parallel_mode && !use_priority_buffer && input_cfg.batch_size == 1);

    if (use_parallel_mode) {
        return NodePartitioningMode::BufferedPriorityParallel;
    }
    if (use_priority_buffer) {
        return NodePartitioningMode::BufferedPrioritySequential;
    }
    if (direct_fennel) {
        return NodePartitioningMode::DirectOnePass;
    }
    return NodePartitioningMode::BufferedNoPriority;
}

} // namespace node_partitioning::mode

