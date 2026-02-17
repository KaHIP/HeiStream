/****************************************************************************
 * node_partitioning_mode_prepare.h
 *****************************************************************************/

#ifndef NODE_PARTITIONING_MODE_PREPARE_H_
#define NODE_PARTITIONING_MODE_PREPARE_H_

#include "algorithms/node_partitioning/mode/node_partitioning_mode_resolver.h"
#include "partition/partition_config.h"

namespace node_partitioning::mode {

// Apply mode-specific derived config values after mode resolution.
void prepare_execution_config_for_mode(Config& config, NodePartitioningMode mode);

}  // namespace node_partitioning::mode

#endif /* NODE_PARTITIONING_MODE_PREPARE_H_ */
