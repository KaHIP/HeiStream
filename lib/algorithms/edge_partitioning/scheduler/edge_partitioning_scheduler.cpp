/****************************************************************************
 * edge_partitioning_scheduler.cpp
 *****************************************************************************/

#include "algorithms/edge_partitioning/scheduler/edge_partitioning_scheduler.h"

#include "algorithms/edge_partitioning/execution/edge_execution.h"
#include "partition/state/edge_partitioner_pass_state.h"


namespace edge_partitioning::scheduler {

void EdgePartitioningScheduler::execute(Config& config, const std::string& graph_filename,
                                        StreamContext& ctx) const {
    mode::EdgePartitioningModeResolver mode_resolver;
    const mode::EdgePartitioningMode mode = mode_resolver.resolve(config);
    switch (mode) {
        case mode::EdgePartitioningMode::Buffered:
            run_edge_stream(config, graph_filename, ctx, ctx.edge_pass_state);
            break;
    }
}

} // namespace edge_partitioning::scheduler

