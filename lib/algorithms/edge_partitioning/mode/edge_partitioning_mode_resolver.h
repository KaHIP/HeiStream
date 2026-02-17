/****************************************************************************
 * edge_partitioning_mode_resolver.h
 *****************************************************************************/

#ifndef EDGE_PARTITIONING_MODE_RESOLVER_H_
#define EDGE_PARTITIONING_MODE_RESOLVER_H_

#include "partition/partition_config.h"


namespace edge_partitioning::mode {

enum class EdgePartitioningMode {
    Buffered,
};

class EdgePartitioningModeResolver {
   public:
    // Resolve execution mode from config without mutating any configuration.
    EdgePartitioningMode resolve(const Config& config) const;
};

} // namespace edge_partitioning::mode


#endif /* EDGE_PARTITIONING_MODE_RESOLVER_H_ */
