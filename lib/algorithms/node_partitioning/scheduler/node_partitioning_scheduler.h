/****************************************************************************
 * node_partitioning_scheduler.h
 *****************************************************************************/

#ifndef NODE_PARTITIONING_SCHEDULER_H_
#define NODE_PARTITIONING_SCHEDULER_H_

#include <string>

#include "core/context/stream_context.h"
#include "core/interfaces/algorithm_scheduler.h"

namespace node_partitioning::scheduler {

class NodePartitioningScheduler : public core::interfaces::IAlgorithmScheduler {
   public:
    void execute(Config& config, const std::string& graph_filename,
                 StreamContext& ctx) const override;
};

} // namespace node_partitioning::scheduler


#endif /* NODE_PARTITIONING_SCHEDULER_H_ */
