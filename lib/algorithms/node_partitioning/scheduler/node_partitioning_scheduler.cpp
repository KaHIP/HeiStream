/****************************************************************************
 * node_partitioning_scheduler.cpp
 *****************************************************************************/

#include "algorithms/node_partitioning/scheduler/node_partitioning_scheduler.h"

#include "algorithms/node_partitioning/execution/node_execution_internal.h"
#include "algorithms/node_partitioning/mode/node_partitioning_mode_prepare.h"
#include "algorithms/node_partitioning/mode/node_partitioning_mode_resolver.h"

namespace node_partitioning::scheduler {

void NodePartitioningScheduler::execute(Config& config, const std::string& graph_filename,
                                        StreamContext& ctx) const {
    node_partitioning::mode::NodePartitioningModeResolver mode_resolver;
    const node_partitioning::mode::NodePartitioningMode mode = mode_resolver.resolve(config);
    node_partitioning::mode::prepare_execution_config_for_mode(config, mode);

    switch (mode) {
        case node_partitioning::mode::NodePartitioningMode::BufferedPriorityParallel:
            stream_node_modes::run_node_stream_priority_buffer_parallel(config, graph_filename, ctx,
                                                                        ctx.node_pass_state);
            break;
        case node_partitioning::mode::NodePartitioningMode::BufferedPrioritySequential:
            stream_node_modes::run_node_stream_priority_buffer_sequential(
                config, graph_filename, ctx, ctx.node_pass_state);
            break;
        case node_partitioning::mode::NodePartitioningMode::DirectOnePass:
            stream_node_modes::run_node_stream_direct_fennel(config, graph_filename, ctx,
                                                             ctx.node_pass_state);
            break;
        case node_partitioning::mode::NodePartitioningMode::BufferedNoPriority:
            stream_node_modes::run_node_stream_no_priority_buffer(config, graph_filename, ctx,
                                                                  ctx.node_pass_state);
            break;
    }
}

} // namespace node_partitioning::scheduler
