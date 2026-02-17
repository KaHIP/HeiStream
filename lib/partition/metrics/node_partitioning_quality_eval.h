/****************************************************************************
 * node_partitioning_quality_eval.h
 *****************************************************************************/

#ifndef PARTITION_METRICS_NODE_PARTITIONING_QUALITY_EVAL_H_
#define PARTITION_METRICS_NODE_PARTITIONING_QUALITY_EVAL_H_

#include <string>
#include <vector>

#include "definitions.h"
#include "partition/partition_config.h"



namespace partition::metrics::node {

// Evaluates edge cut for a node partition result.
void evaluate_partition(const Config& config, const std::vector<PartitionID>& stream_nodes_assign,
                        const std::string& filename, EdgeWeight& edgeCut);

} // namespace partition::metrics::node



#endif /* PARTITION_METRICS_NODE_PARTITIONING_QUALITY_EVAL_H_ */
