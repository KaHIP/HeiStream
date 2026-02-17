/******************************************************************************
 * edge_partitioning_quality_eval.h
 *****************************************************************************/

#ifndef PARTITION_METRICS_EDGE_PARTITIONING_QUALITY_EVAL_H_
#define PARTITION_METRICS_EDGE_PARTITIONING_QUALITY_EVAL_H_

#include <string>
#include <vector>

#include "definitions.h"
#include "partition/partition_config.h"



namespace partition::metrics::edge {

// Evaluates edge partition quality by loading the full graph (batch mode).
void evaluate_partition_batch(Config& config, const std::string& filename,
                              std::vector<PartitionID>* stream_nodes_assign, NodeID& vertex_cut,
                              NodeID& replicas, double& replication_factor, double& balance);

// Evaluates edge partition quality by streaming the input graph.
void evaluate_partition(Config& config, const std::string& filename,
                        std::vector<PartitionID>* stream_nodes_assign, NodeID& vertex_cut,
                        NodeID& replicas, double& replication_factor, double& balance);

} // namespace partition::metrics::edge



#endif /* PARTITION_METRICS_EDGE_PARTITIONING_QUALITY_EVAL_H_ */
