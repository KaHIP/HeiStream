/****************************************************************************
 * node_partitioning_mode_resolver.h
 *****************************************************************************/

#ifndef NODE_PARTITIONING_MODE_RESOLVER_H_
#define NODE_PARTITIONING_MODE_RESOLVER_H_

#include "partition/partition_config.h"


namespace node_partitioning::mode {

enum class NodePartitioningMode {
    DirectOnePass,
    BufferedNoPriority,
    BufferedPrioritySequential,
    BufferedPriorityParallel,
};

class NodePartitioningModeResolver {
   public:
    // Resolve execution mode from config without mutating any configuration.
    NodePartitioningMode resolve(const Config& config) const;
};

} // namespace node_partitioning::mode


#endif /* NODE_PARTITIONING_MODE_RESOLVER_H_ */
