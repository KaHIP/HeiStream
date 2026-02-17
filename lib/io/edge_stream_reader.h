/****************************************************************************
 * edge_stream_reader.h
 *****************************************************************************/

#ifndef IO_EDGE_STREAM_READER_H_
#define IO_EDGE_STREAM_READER_H_

#include <string>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "io/stream_reader.h"
#include "partition/partition_config.h"
#include "partition/state/edge_partitioner_pass_state.h"


namespace io::edge_stream {

/**
 * Read one text batch and build the temporary split graph input.
 *
 * @param partition_config Mutable run configuration.
 * @param pass_state Mutable edge pass state.
 * @param stream_reader Stream reader positioned at current batch start.
 * @param rev_edge Output reverse-edge mapping for split graph edges.
 * @param contracted_edge_graph_id Output split-edge to contracted-edge id mapping.
 * @param number_of_deg1_vertices Output count of degree-1 vertices in split graph.
 * @param number_of_deg2_vertices Output count of degree-2 vertices in split graph.
 * @param G_temp Output temporary split graph for this batch.
 */
void read_input_as_graph(Config& partition_config,
                         partition::state::EdgePartitionerPassState& pass_state,
                         io::stream::StreamReader& stream_reader, std::vector<EdgeID>& rev_edge,
                         std::vector<EdgeID>& contracted_edge_graph_id,
                         NodeID& number_of_deg1_vertices, NodeID& number_of_deg2_vertices,
                         graph_access& G_temp);

/**
 * Read one binary batch and build the temporary split graph input.
 *
 * @param partition_config Mutable run configuration.
 * @param pass_state Mutable edge pass state.
 * @param stream_reader Stream reader positioned at current batch start.
 * @param rev_edge Output reverse-edge mapping for split graph edges.
 * @param contracted_edge_graph_id Output split-edge to contracted-edge id mapping.
 * @param number_of_deg1_vertices Output count of degree-1 vertices in split graph.
 * @param number_of_deg2_vertices Output count of degree-2 vertices in split graph.
 * @param G_temp Output temporary split graph for this batch.
 */
void read_input_as_graph_binary(Config& partition_config,
                                partition::state::EdgePartitionerPassState& pass_state,
                                io::stream::StreamReader& stream_reader,
                                std::vector<EdgeID>& rev_edge,
                                std::vector<EdgeID>& contracted_edge_graph_id,
                                NodeID& number_of_deg1_vertices, NodeID& number_of_deg2_vertices,
                                graph_access& G_temp);

} // namespace io::edge_stream


#endif /* IO_EDGE_STREAM_READER_H_ */
