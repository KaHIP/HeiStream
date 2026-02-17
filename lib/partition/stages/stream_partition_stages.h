/****************************************************************************
 * stream_partition_stages.h
 *****************************************************************************/

#ifndef STREAM_PARTITION_STAGES_H_
#define STREAM_PARTITION_STAGES_H_

#include "partition/partition_config.h"
#include "partition/state/edge_partitioner_pass_state.h"
#include "partition/state/node_partitioner_pass_state.h"

class graph_access;
class balance_configuration;

namespace partition_stages {

// Switch initial partitioning strategy after the first batch when multibfs is enabled.
void config_multibfs_initial_partitioning(Config& config, int curr_batch);

// Shared pre-partition preparation used by streaming decomposition algorithms.
void prepare_batch_partition(Config& config, int restream_number, int curr_batch,
                             const std::vector<NodeWeight>& stream_blocks_weight, graph_access& G,
                             balance_configuration& bc);

// Shared partition-execution stage.
void execute_partition(Config& config, graph_access& G, int restream_number);

}  // namespace partition_stages

#endif /* STREAM_PARTITION_STAGES_H_ */
