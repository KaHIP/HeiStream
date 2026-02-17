/****************************************************************************
 * edge_partitioning_scheduler.h
 *****************************************************************************/

#ifndef EDGE_PARTITIONING_SCHEDULER_H_
#define EDGE_PARTITIONING_SCHEDULER_H_

#include <string>

#include "algorithms/edge_partitioning/mode/edge_partitioning_mode_resolver.h"
#include "core/context/stream_context.h"
#include "core/interfaces/algorithm_scheduler.h"


namespace edge_partitioning::scheduler {

class EdgePartitioningScheduler : public core::interfaces::IAlgorithmScheduler {
   public:
    void execute(Config& config, const std::string& graph_filename,
                 StreamContext& ctx) const override;
};

} // namespace edge_partitioning::scheduler


#endif /* EDGE_PARTITIONING_SCHEDULER_H_ */
