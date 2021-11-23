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

#include "definitions.h"
#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "timer.h"
#include "random_functions.h"
#include "data_structure/buffered_map.h"

typedef std::vector<std::string>* LINE_BUFFER;

class graph_io_stream {
        public:
                graph_io_stream();
                virtual ~graph_io_stream () ;

                static
		NodeID createModel (PartitionConfig & config, graph_access & G, std::vector<std::vector<LongNodeID>>* &input);

                static
		void processNodeWeight(PartitionConfig & config, std::vector<NodeWeight>& all_nodes, NodeID node, NodeWeight weight);

                static
                void generalizeStreamPartition(PartitionConfig & config, graph_access & G_local);

                static
		void countAssignedNodes(PartitionConfig & config);

                static
                void onePassPartition(PartitionConfig & config, std::vector<std::vector<EdgeWeight>> & edges_virtualReal,
					std::vector<PartitionID> & blockVirtualToReal, std::vector<NodeWeight> & weight_VirtualBlocks);

		static
		int onePassDecide(PartitionConfig & config, NodeID node, std::vector<EdgeWeight> & edges_i_real);

		static
		double getFennelWeight(PartitionConfig & partition_config);

                static
		void writePartitionStream(PartitionConfig & config, const std::string & filename);

                static
		void readFirstLineStream(PartitionConfig & partition_config, std::string graph_filename, EdgeWeight& total_edge_cut);
                
                static
		void loadRemainingLines(PartitionConfig & partition_config, LINE_BUFFER &lines);

                static
		void loadBufferLines(PartitionConfig & partition_config, LINE_BUFFER &lines, LongNodeID num_lines);

                static
		std::vector<std::string>* loadLinesFromStream(PartitionConfig & partition_config, LongNodeID num_lines);

                static
		void sortBatchByDegree(PartitionConfig & config);

                static
		void createGraphForBatch(PartitionConfig & config, graph_access & G, NodeID node_counter, EdgeID edge_counter,
			std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes);

                static
		void recoverBlockAssignedToNode(PartitionConfig & config, graph_access & G, NodeID node, NodeID node_counter);

                static
		void setupForGhostNeighbors(PartitionConfig & config);

                static
		void processGhostNeighborInBatch(PartitionConfig & config, NodeID node, LongNodeID ghost_target, EdgeWeight edge_weight);

                static
		void processQuotientEdgeInBatch(PartitionConfig & config, NodeID node, LongNodeID global_target, EdgeWeight edge_weight);

                static
		EdgeID insertRegularEdgeInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, 
										NodeID node, LongNodeID global_target, EdgeWeight edge_weight);

                static
		NodeID mapGhostKeysToNodesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, 
							std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes, NodeID& node_counter);

                static
		NodeID restreamMapGhostKeysToNodes(PartitionConfig & config);

                static
		NodeID greedyMapGhostKeysToNodes(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, 
						std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes, NodeID& node_counter);

                static
		EdgeID insertGhostEdgesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edgesInBatch);

                static
		void insertQuotientNodesInBatch(PartitionConfig & config, std::vector<NodeWeight>& all_nodes, NodeID uncontracted_ghost_nodes, NodeID& node_counter);

                static
		EdgeID insertQuotientEdgesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, NodeID uncontracted_ghost_nodes);


                static
		EdgeID includeEdgeInBatch(std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, NodeID node, NodeID target, EdgeWeight edge_weight);

                static
		void prescribeBufferInbalance(PartitionConfig & partition_config);

                static
		void streamEvaluatePartition(PartitionConfig & config, const std::string & filename, EdgeWeight& edgeCut);

                static
		void loadRemainingLinesToBinary(PartitionConfig & partition_config, std::vector<std::vector<LongNodeID>>* &input);

                static
		void loadBufferLinesToBinary(PartitionConfig & partition_config, std::vector<std::vector<LongNodeID>>* &input, LongNodeID num_lines);

                static
		std::vector<std::vector<LongNodeID>>* loadLinesFromStreamToBinary(PartitionConfig & partition_config, LongNodeID num_lines);

		template< typename T>
                static
		T return_and_delete_element(std::vector<T> & vec, LongNodeID pos);

//                static
//		int readEdgeStream_writeMetisBuffered(const std::string & graph_filename, std::string filename_output, bool relabel_nodes);

};

inline void graph_io_stream::loadRemainingLinesToBinary(PartitionConfig & partition_config, std::vector<std::vector<LongNodeID>>* &input) {
	if (partition_config.ram_stream) {                                                                    
		input = graph_io_stream::loadLinesFromStreamToBinary(partition_config, partition_config.remaining_stream_nodes);
	}                                                                                                     
}                                                                                                             

inline void graph_io_stream::loadBufferLinesToBinary(PartitionConfig & partition_config, std::vector<std::vector<LongNodeID>>* &input, LongNodeID num_lines) {
	if (!partition_config.ram_stream) {                                                                   
		input = graph_io_stream::loadLinesFromStreamToBinary(partition_config, num_lines);            
	}                                                                                                     
}                                                                                                             

inline std::vector<std::vector<LongNodeID>>* graph_io_stream::loadLinesFromStreamToBinary(PartitionConfig & partition_config, LongNodeID num_lines) {
	std::vector<std::vector<LongNodeID>>* input;                                                          
	input = new std::vector<std::vector<LongNodeID>>(num_lines);                                                                                          
	std::vector<std::string>* lines;                                                                      
	lines = new std::vector<std::string>(1);                                                              
	LongNodeID node_counter = 0;                                                                          
	buffered_input *ss2 = NULL;                                                                           
	while( node_counter < num_lines) {                                                                    
		std::getline(*(partition_config.stream_in),(*lines)[0]);                                      
		if ((*lines)[0][0] == '%') { // a comment in the file                                         
			continue;                                                                             
		}                                                                                             
		ss2 = new buffered_input(lines);                                                              
		ss2->simple_scan_line((*input)[node_counter++]);                                              
		(*lines)[0].clear(); delete ss2;                                                              
	}                                                                                                     
	delete lines;                                                                                         
	return input;                                                                                         
}                                                                                                             



#endif /*GRAPHIOSTREAM_H_*/
