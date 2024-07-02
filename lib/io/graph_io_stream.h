/******************************************************************************
 * graph_io_stream.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPHIOSTREAM_H_
#define GRAPHIOSTREAM_H_

#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <set>
#include <unordered_map>
#include <list>
#include <algorithm>
#include <cstdint>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "definitions.h"
#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "timer.h"
#include "random_functions.h"
#include "data_structure/buffered_map.h"
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include "cpi/run_length_compression.hpp"
#include "data_structure/compression_vectors/CompressionDataStructure.h"
#include "data_structure/compression_vectors/RunLengthCompressionVector.h"
#include "data_structure/compression_vectors/BatchRunLengthCompression.h"

typedef std::vector <std::string> *LINE_BUFFER;

class graph_io_stream {
public:
    graph_io_stream();

    virtual ~graph_io_stream();

    static
    NodeID createModel(PartitionConfig &config, const std::shared_ptr <CompressionDataStructure<PartitionID>> &block_assignments, graph_access &G, std::vector <std::vector<LongNodeID>> *&input);

    static
    void
    processNodeWeight(PartitionConfig &config, const std::shared_ptr <CompressionDataStructure<PartitionID>> &block_assignments, std::vector <NodeWeight> &all_nodes, NodeID node, NodeWeight weight);

    static
    void generalizeStreamPartition(PartitionConfig &config, const std::shared_ptr <CompressionDataStructure<PartitionID>> &block_assignments, graph_access &G_local);

    static
    void countAssignedNodes(PartitionConfig &config);

    static
    void onePassPartition(PartitionConfig &config, std::vector <std::vector<EdgeWeight>> &edges_virtualReal,
                          std::vector <PartitionID> &blockVirtualToReal,
                          std::vector <NodeWeight> &weight_VirtualBlocks);

    static
    int onePassDecide(PartitionConfig &config, NodeID node, std::vector <EdgeWeight> &edges_i_real);

    static
    double getFennelWeight(PartitionConfig &partition_config);

    static
    void writePartitionStream(PartitionConfig &config, const std::shared_ptr <CompressionDataStructure<PartitionID>> &block_assignments);

    static
    void writePartitionStream(PartitionConfig &config);

    static
    void readFirstLineStream(PartitionConfig &partition_config, std::string graph_filename, EdgeWeight &total_edge_cut);

    static
    void loadRemainingLines(PartitionConfig &partition_config, LINE_BUFFER &lines);

    static
    void loadBufferLines(PartitionConfig &partition_config, LINE_BUFFER &lines, LongNodeID num_lines);

    static
    std::vector <std::string> *loadLinesFromStream(PartitionConfig &partition_config, LongNodeID num_lines);

    static
    void sortBatchByDegree(PartitionConfig &config);

    static
    void createGraphForBatch(PartitionConfig &config, graph_access &G, NodeID node_counter, EdgeID edge_counter,
                             std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    >& all_edges,
    std::vector <NodeWeight> &all_nodes, std::vector<NodeWeight>
    & all_assigned_ghost_nodes);

    static
    void recoverBlockAssignedToNode(PartitionConfig &config, graph_access &G, NodeID node, NodeID node_counter);

    static
    void setupForGhostNeighbors(PartitionConfig &config);

    static
    void
    processGhostNeighborInBatch(PartitionConfig &config, NodeID node, LongNodeID ghost_target, EdgeWeight edge_weight);

    static
    void
    processQuotientEdgeInBatch(PartitionConfig &config, const std::shared_ptr <CompressionDataStructure<PartitionID>> &block_assignments, NodeID node, LongNodeID global_target, EdgeWeight edge_weight);

    static
    EdgeID insertRegularEdgeInBatch(PartitionConfig &config, std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    >& all_edges,
    NodeID node, LongNodeID
    global_target,
    EdgeWeight edge_weight
    );

    static
    NodeID mapGhostKeysToNodesInBatch(PartitionConfig &config, std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    >& all_edges,
    std::vector <NodeWeight> &all_nodes, std::vector<NodeWeight>
    & all_assigned_ghost_nodes,
    NodeID &node_counter
    );

    static
    NodeID restreamMapGhostKeysToNodes(PartitionConfig &config);

    static
    NodeID greedyMapGhostKeysToNodes(PartitionConfig &config, std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    >& all_edges,
    std::vector <NodeWeight> &all_nodes, std::vector<NodeWeight>
    & all_assigned_ghost_nodes,
    NodeID &node_counter
    );

    static
    EdgeID insertGhostEdgesInBatch(PartitionConfig &config, std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    >& all_edgesInBatch);

    static
    void insertQuotientNodesInBatch(PartitionConfig &config, std::vector <NodeWeight> &all_nodes,
                                    NodeID uncontracted_ghost_nodes, NodeID &node_counter);

    static
    EdgeID insertQuotientEdgesInBatch(PartitionConfig &config, std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    >& all_edges,
    NodeID uncontracted_ghost_nodes
    );


    static
    EdgeID includeEdgeInBatch(std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    >& all_edges,
    NodeID node, NodeID
    target,
    EdgeWeight edge_weight
    );

    static
    void prescribeBufferInbalance(PartitionConfig &partition_config);

    static
    void streamEvaluatePartition(PartitionConfig &config, const std::shared_ptr <CompressionDataStructure<PartitionID>> &block_assignments, const std::string &filename, EdgeWeight &edgeCut);

    static
    void loadRemainingLinesToBinary(PartitionConfig &partition_config, std::vector <std::vector<LongNodeID>> *&input);

    static
    void loadBufferLinesToBinary(PartitionConfig &partition_config, std::vector <std::vector<LongNodeID>> *&input,
                                 LongNodeID num_lines);

    static
    std::vector <std::vector<LongNodeID>> *
    loadLinesFromStreamToBinary(PartitionConfig &partition_config, LongNodeID num_lines);

    template<typename T>
    static
    T return_and_delete_element(std::vector <T> &vec, LongNodeID pos);

    // for edge partitioning

    static void readFirstLineStreamEdge(PartitionConfig &partition_config,
                                        std::string graph_filename,
                                        EdgeWeight &total_edge_cut);

    static void generalizeStreamPartitionEdge(PartitionConfig &config,
                                              graph_access &G_local);

    static bool hasEnding(std::string const &string, std::string const &ending);

    static void createGraphForBatchEdge(
            PartitionConfig &config, graph_access &G, NodeID node_counter,
            EdgeID edge_counter,
            std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    > &all_edges,
    std::vector <NodeWeight> &all_nodes
    );

    static void processQuotientEdgeInBatchEdge(PartitionConfig &config,
                                               NodeID node,
                                               LongNodeID global_target,
                                               EdgeWeight edge_weight);

    static EdgeID insertRegularEdgeInBatchEdge(
            PartitionConfig &config,
            std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    > &all_edges,
    NodeID node, LongNodeID
    global_target,
    EdgeWeight edge_weight
    );

    static void insertQuotientNodesInBatchEdge(PartitionConfig &config,
                                               std::vector <NodeWeight> &all_nodes,
                                               NodeID &node_counter);

    static EdgeID insertQuotientEdgesInBatchEdge(
            PartitionConfig &config,
            std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    > &all_edges);

    static EdgeID includeEdgeInBatchEdge(
            std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

    > &all_edges,
    NodeID node, NodeID
    target,
    EdgeWeight edge_weight
    );

    static void streamEvaluatePartitionEdgeBatch(PartitionConfig &config,
                                                 const std::string &filename,
                                                 NodeID &vertex_cut,
                                                 NodeID &replicas,
                                                 double &replication_factor,
                                                 double &balance);

    static void
    streamEvaluateEdgePartition(PartitionConfig &config, const std::string &filename, NodeID &vertex_cut,
                                NodeID &replicas, double &replication_factor, double &balance);

    static void find_reverse_edges(graph_access &G,
                                   std::vector <EdgeID> &m_reverse_edge);

    static std::vector <PartitionID>
    reservoir_sampling(std::vector <PartitionID> &v, size_t n);


    static void buildGraphModel(graph_access &G, PartitionConfig &partition_config,
                                std::vector <EdgeID> &contracted_edge_graph_id,
                                std::vector <EdgeID> &rev_edge, EdgeID S_nodes,
                                NodeID S_edges, graph_access *G_temp);


    static void updateBatchVariables(PartitionConfig &partition_config);

    static void constructBatchIOGraph(
            PartitionConfig &partition_config,
            std::vector <std::vector<Edge>> &subgraph_edges,
            std::vector <NodeID> &mapping, std::vector <EdgeID> &rev_edge,
            std::vector <EdgeID> &contracted_edge_graph_id,
            NodeID &number_of_deg1_vertices, NodeID &number_of_deg2_vertices,
            NodeID &G_temp_no_nodes, NodeID &G_temp_no_edges, graph_access &G_temp);


    static void readInputAsGraph(PartitionConfig &partition_config,
                                 std::vector <EdgeID> &rev_edge,
                                 std::vector <EdgeID> &contracted_edge_graph_id,
                                 NodeID &number_of_deg1_vertices,
                                 NodeID &number_of_deg2_vertices,
                                 graph_access &G_temp);

    static void readInputAsGraphBinary(PartitionConfig &partition_config,
                                       std::vector <EdgeID> &rev_edge,
                                       std::vector <EdgeID> &contracted_edge_graph_id,
                                       NodeID &number_of_deg1_vertices,
                                       NodeID &number_of_deg2_vertices,
                                       graph_access &G_temp);


    static void constructBatchModel(PartitionConfig &partition_config,
                                    graph_access &G);

    static void addEdgesToIncidentBlock(PartitionConfig &partition_config,
                                        NodeID node, NodeID target,
                                        EdgeWeight edge_weight,
                                        LongEdgeID &used_edges);

    static void writeStreamOutput(PartitionConfig &config, graph_access &G_local);


//                static
//		int readEdgeStream_writeMetisBuffered(const std::string & graph_filename, std::string filename_output, bool relabel_nodes);

};

inline void graph_io_stream::loadRemainingLinesToBinary(PartitionConfig &partition_config,
                                                        std::vector <std::vector<LongNodeID>> *&input) {
    if (partition_config.ram_stream) {
        input = graph_io_stream::loadLinesFromStreamToBinary(partition_config, partition_config.remaining_stream_nodes);
    }
}

inline void graph_io_stream::loadBufferLinesToBinary(PartitionConfig &partition_config,
                                                     std::vector <std::vector<LongNodeID>> *&input,
                                                     LongNodeID num_lines) {
    if (!partition_config.ram_stream) {
        input = graph_io_stream::loadLinesFromStreamToBinary(partition_config, num_lines);
    }
}

inline std::vector <std::vector<LongNodeID>> *
graph_io_stream::loadLinesFromStreamToBinary(PartitionConfig &partition_config, LongNodeID num_lines) {
    std::vector <std::vector<LongNodeID>> *input;
    input = new std::vector <std::vector<LongNodeID>>(num_lines);
    std::vector <std::string> *lines;
    lines = new std::vector<std::string>(1);
    LongNodeID node_counter = 0;
    buffered_input *ss2 = NULL;
    while (node_counter < num_lines) {
        std::getline(*(partition_config.stream_in), (*lines)[0]);
        if ((*lines)[0][0] == '%') { // a comment in the file
            continue;
        }
        ss2 = new buffered_input(lines);
        ss2->simple_scan_line((*input)[node_counter++]);
        (*lines)[0].clear();
        delete ss2;
    }
    delete lines;
    return input;
}


#endif /*GRAPHIOSTREAM_H_*/
