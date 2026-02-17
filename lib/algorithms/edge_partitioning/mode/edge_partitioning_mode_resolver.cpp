/****************************************************************************
 * edge_partitioning_mode_resolver.cpp
 *****************************************************************************/

#include "algorithms/edge_partitioning/mode/edge_partitioning_mode_resolver.h"


namespace edge_partitioning::mode {

EdgePartitioningMode EdgePartitioningModeResolver::resolve(const Config& config) const {
    (void)config;
    return EdgePartitioningMode::Buffered;
}

} // namespace edge_partitioning::mode

