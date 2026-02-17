/****************************************************************************
 * edge_batch_model_builder.h
 *****************************************************************************/

#ifndef ALGORITHMS_EDGE_PARTITIONING_MODEL_EDGE_BATCH_MODEL_BUILDER_H_
#define ALGORITHMS_EDGE_PARTITIONING_MODEL_EDGE_BATCH_MODEL_BUILDER_H_

#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "io/stream_reader.h"
#include "partition/partition_config.h"
#include "partition/state/edge_partitioner_pass_state.h"


namespace edge_partitioning::model {

/**
 * Prepare per-batch edge conversion window and edge-id offsets.
 *
 * @param config Mutable run configuration.
 * @param pass_state Mutable edge pass state.
 */
void prepare_edge_batch_window(Config& config,
                               partition::state::EdgePartitionerPassState& pass_state);

/**
 * Build split-graph batch representation and edge-id mappings.
 *
 * @param config Mutable run configuration.
 * @param pass_state Mutable edge pass state.
 * @param subgraph_edges Temporary adjacency storage for split graph.
 * @param mapping Local split-node to original-global-node mapping.
 * @param rev_edge Output reverse-edge mapping.
 * @param contracted_edge_graph_id Output contracted-edge id mapping.
 * @param number_of_deg1_vertices Output degree-1 split vertex count.
 * @param number_of_deg2_vertices Output degree-2 split vertex count.
 * @param g_temp_no_nodes Number of split nodes to materialize.
 * @param g_temp_no_edges Number of split edges to materialize.
 * @param g_temp Output split graph.
 */
void build_edge_split_batch_graph(Config& config,
                                  partition::state::EdgePartitionerPassState& pass_state,
                                  std::vector<std::vector<Edge>>& subgraph_edges,
                                  std::vector<NodeID>& mapping, std::vector<EdgeID>& rev_edge,
                                  std::vector<EdgeID>& contracted_edge_graph_id,
                                  NodeID& number_of_deg1_vertices, NodeID& number_of_deg2_vertices,
                                  NodeID g_temp_no_nodes, NodeID g_temp_no_edges,
                                  graph_access& g_temp);

/**
 * Build the edge-partitioning model graph for the current batch.
 *
 * @param config Mutable run configuration.
 * @param pass_state Mutable edge pass state.
 * @param graph Output model graph for this batch.
 * @param stream_reader Stream reader positioned at batch start.
 */
void build_batch_model(Config& config, partition::state::EdgePartitionerPassState& pass_state,
                       graph_access& graph, io::stream::StreamReader& stream_reader);

} // namespace edge_partitioning::model


#endif /* ALGORITHMS_EDGE_PARTITIONING_MODEL_EDGE_BATCH_MODEL_BUILDER_H_ */
