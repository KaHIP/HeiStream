/******************************************************************************
 * graph_io_stream.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <math.h>
#include <sstream>
#include "graph_io_stream.h"
#include "timer.h"

#define MIN(A,B) ((A)<(B))?(A):(B)
#define MAX(A,B) ((A)>(B))?(A):(B)


graph_io_stream::graph_io_stream() {

}

graph_io_stream::~graph_io_stream() {

}


NodeID graph_io_stream::createModel (PartitionConfig & config, graph_access & G, std::vector<std::vector<LongNodeID>>* &input) {
	NodeWeight total_nodeweight = 0;                                                                      
	NodeID node_counter = 0;                                                                              
	EdgeID edge_counter = 0;                                                                              
	LongEdgeID used_edges = 0;                                                                            
	bool read_ew = false;                                                                                 
	bool read_nw = false;                                                                                 
	LongEdgeID nmbEdges = 2*config.remaining_stream_edges;                                                
	LongNodeID target;                                                                                    
	NodeWeight weight;                                                                                    
	std::vector<std::vector<std::pair<NodeID,EdgeWeight>>> all_edges;                                     
	std::vector<NodeWeight> all_nodes;                                                                    

	std::vector<NodeWeight> all_assigned_ghost_nodes(config.nmbNodes + config.quotient_nodes,0);          
	all_edges.resize(config.nmbNodes + config.quotient_nodes);                                            
	all_nodes.resize(config.nmbNodes + config.quotient_nodes);                                            
	config.lower_global_node = config.total_stream_nodecounter + 1; // Bounds below start from 1 instead of 0
	config.upper_global_node = config.total_stream_nodecounter + config.nmbNodes;                         
	LongNodeID cursor = 0;                                                                                
	NodeID node = 0;                                                                                      

	config.curr_batch++;                                                                                  

	if (config.ram_stream) {                                                                              
		cursor = input->size() - config.remaining_stream_nodes;                                       
	}                                                                                                     


	if( nmbEdges > std::numeric_limits<EdgeWeight>::max() || config.nmbNodes > std::numeric_limits<LongNodeID>::max()) {
#ifdef MODE64BITEDGES                                                                                         
		std::cerr <<  "The graph is too large. Currently only 64bits supported!"  << std::endl;       
#else                                                                                                         
		std::cerr <<  "The graph is too large. Currently only 32bits supported!"  << std::endl;       
#endif                                                                                                        
		exit(0);                                                                                      
	}                                                                                                     

	switch(config.remaining_stream_ew) {                                                                  
		case 1:                                                                                       
			read_ew = true;                                                                       
			break;                                                                                
		case 10:                                                                                      
			read_nw = true;                                                                       
			break;                                                                                
		case 11:                                                                                      
			read_ew = true;                                                                       
			read_nw = true;                                                                       
			break;                                                                                
	}                                                                                                     

	/* config.degree_nodeBlock = new std::vector<std::vector<EdgeWeight>> (config.nmbNodes, std::vector<EdgeWeight>(config.k,0)); */
	config.edge_block_nodes = new std::vector<std::vector<NodeID>> (config.k, std::vector<NodeID>());

	setupForGhostNeighbors(config);                                                                       

	for (node_counter=0; node_counter < config.nmbNodes; node_counter++) {                                
		std::vector<LongNodeID> &line_numbers = (*input)[cursor];                                     
		LongNodeID col_counter = 0;                                                                   
		node = (NodeID) node_counter;                                                                 
		weight = 1;                                                                                   
		if( read_nw ) {                                                                               
			weight = line_numbers[col_counter++];                                                 
			if( total_nodeweight > std::numeric_limits<NodeWeight>::max()) {                      
				std::cerr <<  "The sum of the node weights is too large (it exceeds the node weight type)."  << std::endl;
				std::cerr <<  "Currently not supported. Please scale your node weights."  << std::endl;
				exit(0);                                                                      
			}                                                                                     
		}                                                                                             
		total_nodeweight += weight;                                                                   
		processNodeWeight(config, all_nodes, node, weight);                                           

		while (col_counter < line_numbers.size()) {                                                   
			target = line_numbers[col_counter++];                                                 
			EdgeWeight edge_weight = 1;                                                           
			if( read_ew ) {                                                                       
				edge_weight = line_numbers[col_counter++];                                    
			}                                                                                     

			if(target > config.upper_global_node) { // edge to future batch                       
				processGhostNeighborInBatch(config, node, target, edge_weight);               
			} else if(target < config.lower_global_node ) { // edge to previous batch             
				used_edges++;                                                                 
				processQuotientEdgeInBatch(config, node, target, edge_weight);                
			} else { // edge to current batch                                                     
				used_edges += ((NodeID)(target - config.lower_global_node) < node); // used_edges only counts arcs to previus nodes
				edge_counter += insertRegularEdgeInBatch(config, all_edges, node, target, edge_weight);
			}                                                                                     
		}                                                                                             

		cursor++;                                                                                     
	}                                                                                                     
	if (!config.ram_stream) {                                                                             
		delete input;                                                                                 
	}                                                                                                     

	NodeID uncontracted_ghost_nodes = mapGhostKeysToNodesInBatch(config, all_edges, all_nodes, all_assigned_ghost_nodes, node_counter);
	insertQuotientNodesInBatch(config, all_nodes, uncontracted_ghost_nodes, node_counter);                
	edge_counter += insertGhostEdgesInBatch(config, all_edges);                                           
	edge_counter += insertQuotientEdgesInBatch(config, all_edges, uncontracted_ghost_nodes);              

	createGraphForBatch(config, G, node_counter, edge_counter, all_edges, all_nodes, all_assigned_ghost_nodes);

	/* delete config.degree_nodeBlock; */
	delete config.edge_block_nodes;


	config.total_stream_nodecounter += config.nmbNodes;                                                   
	config.total_stream_nodeweight  += total_nodeweight;                                                  
	config.remaining_stream_nodes   -= config.nmbNodes;                                                   
	config.remaining_stream_edges   -= used_edges;                                                        

	if( node_counter != (NodeID) config.nmbNodes + uncontracted_ghost_nodes + config.quotient_nodes) {    
		std::cerr <<  "number of specified nodes mismatch"  << std::endl;                             
		std::cerr <<  (config.nmbNodes + uncontracted_ghost_nodes + config.quotient_nodes) <<  " " <<  node_counter  << std::endl;
		exit(0);                                                                                      
	}                                                                                                     

	return node_counter;                                                                                  
}                                                                                                             




void graph_io_stream::insertQuotientNodesInBatch(PartitionConfig & config, std::vector<NodeWeight>& all_nodes, NodeID uncontracted_ghost_nodes, NodeID& node_counter) {
        while (node_counter < config.nmbNodes + uncontracted_ghost_nodes + config.quotient_nodes) { // Add artificial nodes representing current blocks
                NodeID node = node_counter++;
		NodeID targetPar = node - (NodeID)config.nmbNodes - uncontracted_ghost_nodes;
                NodeWeight weight = (*config.stream_blocks_weight)[targetPar];
		all_nodes[node] = weight;
        }
}


/* EdgeID graph_io_stream::insertQuotientEdgesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, NodeID uncontracted_ghost_nodes) { */
/* 	EdgeID inserted_edges = 0; */
/* 	for (NodeID node=0; node < config.nmbNodes; node++) { */ 
/* 		for (PartitionID block=0; block < config.k; block++) { */ 
/* 			NodeID target = config.nmbNodes + uncontracted_ghost_nodes + block; */
/* 			EdgeWeight edge_weight = (*config.degree_nodeBlock)[node][block]; */
/* 			if (edge_weight > 0) { */
/* 				inserted_edges += includeEdgeInBatch(all_edges, node, target, (1+config.double_non_ghost_edges)*edge_weight); */
/* 			} */
/* 		} */
/* 	} */ 
/* 	return inserted_edges; */
/* } */
EdgeID graph_io_stream::insertQuotientEdgesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, NodeID uncontracted_ghost_nodes) {
	EdgeID inserted_edges = 0;
	NodeID prev_node;
	EdgeWeight edge_weight;
	for (PartitionID block=0; block < config.k; block++) { 
		NodeID target = config.nmbNodes + uncontracted_ghost_nodes + block;
		if ((*config.edge_block_nodes)[block].size() >= 1) {
			prev_node = (*config.edge_block_nodes)[block][0];
			edge_weight = 0;
		} else continue; // if block has no neighbors in batch, continue outer for loop
		for (auto& node : (*config.edge_block_nodes)[block]) {
			if (node == prev_node) {
				edge_weight++;
			} else {
				inserted_edges += includeEdgeInBatch(all_edges, prev_node, target, (1+config.double_non_ghost_edges)*edge_weight);
				prev_node = node;
				edge_weight = 1;
			}
		}
		// necessary because last neighbor of each block is not directly included inside the inner for loop
		inserted_edges += includeEdgeInBatch(all_edges, prev_node, target, (1+config.double_non_ghost_edges)*edge_weight);
	} 
	return inserted_edges;
}


EdgeID graph_io_stream::insertGhostEdgesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges) { 
	EdgeID inserted_edges = 0;
	if (config.stream_allow_ghostnodes || config.restream_number) {
		for (LongNodeID ghost_key=0; ghost_key < config.ghostkey_to_edges->size(); ghost_key++) {
			NodeID ghost_node = (*config.ghostkey_to_node)[ghost_key];
			auto& neighbors_list = (*config.ghostkey_to_edges)[ghost_key];
			for (auto& edge : neighbors_list) {
				inserted_edges += includeEdgeInBatch(all_edges, ghost_node, edge.first, edge.second);
			}
		}
	}
	return inserted_edges;
}

NodeID graph_io_stream::mapGhostKeysToNodesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, 
					std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes, NodeID& node_counter) {
	NodeID inserted_nodes = 0;
	if (config.restream_number) {
		inserted_nodes = restreamMapGhostKeysToNodes(config);
		config.ghostglobal_to_ghostkey->clear();
	} else if (config.stream_allow_ghostnodes) {
		inserted_nodes = greedyMapGhostKeysToNodes(config, all_edges, all_nodes, all_assigned_ghost_nodes, node_counter);
		config.ghostglobal_to_ghostkey->clear();
	}
	return inserted_nodes;
}

NodeID graph_io_stream::restreamMapGhostKeysToNodes(PartitionConfig & config) {
	if (config.ghostkey_to_node != NULL) {
		delete config.ghostkey_to_node;
	}
	config.ghostkey_to_node = new std::vector<NodeID>(config.ghost_nodes,0);
	for(PartitionID targetPar = 0; targetPar<config.ghost_nodes; targetPar++) {
		NodeID node = targetPar + config.nmbNodes;
		(*config.ghostkey_to_node)[targetPar] = node; // PS: weight of quotient nodes is already up-to-date 
	} 
	return 0;
}

NodeID graph_io_stream::greedyMapGhostKeysToNodes(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, 
						std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes, NodeID& node_counter) {
	if (config.ghostkey_to_node != NULL) 
		delete config.ghostkey_to_node;
	config.ghostkey_to_node = new std::vector<NodeID>(config.ghost_nodes,0);
	NodeID inserted_nodes =  MIN(config.ghost_nodes,config.ghost_nodes_threshold);
	all_nodes.resize(config.nmbNodes + inserted_nodes + config.quotient_nodes);
	all_edges.resize(config.nmbNodes + inserted_nodes + config.quotient_nodes);
	LongNodeID ghost_key;
	for (ghost_key=0; ghost_key<config.ghost_nodes; ghost_key++) {
		if (ghost_key < (LongNodeID) inserted_nodes) { // uncontracted ghost node
			NodeID node = node_counter++;
			(*config.ghostkey_to_node)[ghost_key] = node;
			all_nodes[node] = 1;
		} else { // contracted ghost node
			auto & list_edges = (*config.ghostkey_to_edges)[ghost_key];
//			int contr_pos = random_functions::nextInt(0, list_edges.size()-1);
			LongNodeID contr_pos = ((LongNodeID)ghost_key*config.curr_batch) % list_edges.size();
			NodeID node = return_and_delete_element(list_edges, contr_pos).first; // remove self edge 
			(*config.ghostkey_to_node)[ghost_key] = node;
			all_nodes[node]++; 
			all_assigned_ghost_nodes[node]++; 
		}
	}
	return inserted_nodes;
}

void graph_io_stream::processGhostNeighborInBatch(PartitionConfig & config, NodeID node, LongNodeID ghost_target, EdgeWeight edge_weight) {
	LongNodeID ghost_key;
	PartitionID targetGlobalPar = (*config.stream_nodes_assign)[ghost_target-1];
	if (config.stream_allow_ghostnodes || (config.restream_number && targetGlobalPar != INVALID_PARTITION) ) {
		if (!config.ghostglobal_to_ghostkey->has_key(ghost_target-1)) {
			ghost_key = config.ghost_nodes++;
			config.ghostglobal_to_ghostkey->push_back(ghost_target-1, ghost_key);
			config.ghostkey_to_edges->push_back(std::vector<std::pair<NodeID,ShortEdgeWeight>>());
		} else {
			ghost_key = (*config.ghostglobal_to_ghostkey)[ghost_target-1]; 
		}
		(*config.ghostkey_to_edges)[ghost_key].push_back(std::make_pair(node, (ShortEdgeWeight)edge_weight));
	}
}


/* void graph_io_stream::processQuotientEdgeInBatch(PartitionConfig & config, NodeID node, LongNodeID global_target, EdgeWeight edge_weight) { */
/* 	PartitionID targetGlobalPar = (*config.stream_nodes_assign)[global_target-1]; */
/* 	if ( targetGlobalPar == config.k ) { // SIGNAL: neighbor not yet assigned */
/* 		processGhostNeighborInBatch(config, node, global_target, edge_weight); */ 
/* 		return; */	
/* 	} */
/* 	if ( targetGlobalPar > config.k ) { */
/* 		std::cerr << "ERROR regarding number of parts.\n"; */
/* 		exit(0); */
/* 	} */
/* 	(*config.degree_nodeBlock)[node][targetGlobalPar] += edge_weight; */
/* } */
void graph_io_stream::processQuotientEdgeInBatch(PartitionConfig & config, NodeID node, LongNodeID global_target, EdgeWeight edge_weight) {
	PartitionID targetGlobalPar = (*config.stream_nodes_assign)[global_target-1];
	if ( targetGlobalPar == config.k ) { // SIGNAL: neighbor not yet assigned
		processGhostNeighborInBatch(config, node, global_target, edge_weight); 
		return;	
	}
	if ( targetGlobalPar > config.k ) {
		std::cerr << "ERROR regarding number of parts.\n";
		exit(0);
	}
	for (EdgeWeight i=0; i<edge_weight; i++) (*config.edge_block_nodes)[targetGlobalPar].push_back(node);
}


void graph_io_stream::processNodeWeight(PartitionConfig & config, std::vector<NodeWeight>& all_nodes, NodeID node, NodeWeight weight) {
	LongNodeID global_node = config.lower_global_node + (LongNodeID) node - 1;
	PartitionID nodeGlobalPar = (*config.stream_nodes_assign)[global_node];
	if (config.restream_number) {
		(*config.stream_blocks_weight)[nodeGlobalPar] -= weight;
	}
	all_nodes[node] = weight;
}


EdgeID graph_io_stream::insertRegularEdgeInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, 
										NodeID node, LongNodeID global_target, EdgeWeight edge_weight) {
	NodeID target = (NodeID) (global_target - config.lower_global_node);
	if(target == node) {
		std::cerr <<  "The graph file contains self-loops, which are not supported. Please remove them from the file."  << std::endl;
		exit(0);
	}
	if (target > node) {
		return 0; 
	} 
	 
	return includeEdgeInBatch(all_edges, node, target, (1+config.double_non_ghost_edges)*edge_weight);
}


EdgeID graph_io_stream::includeEdgeInBatch(std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, NodeID node, NodeID target, EdgeWeight edge_weight) {
	all_edges[node].push_back(std::make_pair(target,edge_weight));
	all_edges[target].push_back(std::make_pair(node,edge_weight));

	return 2;
}


void graph_io_stream::setupForGhostNeighbors(PartitionConfig & config) {
	if (config.stream_allow_ghostnodes || config.restream_number) {
		if (config.ghostglobal_to_ghostkey == NULL) {
			config.ghostglobal_to_ghostkey = new buffered_map(config.stream_nodes_assign, config.restream_number);
		}
		if (config.ghostkey_to_edges != NULL)  {
			config.ghostkey_to_edges->clear();
		} else {
			config.ghostkey_to_edges = new std::vector<std::vector<std::pair<NodeID,ShortEdgeWeight>>>();
		}
		config.ghost_nodes = 0 + (config.restream_number>0)*(config.k);
		if (config.restream_number > 0) {
			config.ghost_nodes = config.k;
			config.ghostkey_to_edges->resize(config.ghost_nodes);
		}
	}
}


void graph_io_stream::createGraphForBatch(PartitionConfig & config, graph_access & G, NodeID node_counter, EdgeID edge_counter,
		std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes) {
        G.stream_repeat_construction(node_counter, edge_counter);
	G.resizeImplicitGhostNodesList(node_counter);
	std::vector<EdgeWeight>  neighbor_edgeweight(node_counter,0);
	std::vector<EdgeID> neighbor_edgeid(node_counter,0);
	std::vector<NodeID> neighbors;
	EdgeID e;
	for(NodeID i = 0 ; i < node_counter ; i++) { // actual insertion of nodes and edges
                NodeID node = G.new_node(); 
                NodeWeight weight = all_nodes[node];
                NodeWeight ghost_weight = all_assigned_ghost_nodes[node];
                G.setNodeWeight(node, weight);
                G.setImplicitGhostNodes(node, ghost_weight);
		G.setPartitionIndex(node, 0);
		recoverBlockAssignedToNode(config, G, node, node_counter);
		for (auto& [target,edge_weight] : all_edges[node]) {
			if (neighbor_edgeweight[target] == 0) {
				e = G.new_edge(node, target);
				neighbor_edgeid[target] = e;
				neighbors.push_back(target);
			} else {
				e = neighbor_edgeid[target];
			}
			neighbor_edgeweight[target] += edge_weight;
                        G.setEdgeWeight(e, neighbor_edgeweight[target]);
                }
		for (auto& target : neighbors) {
			neighbor_edgeweight[target] = 0;
		}
		neighbors.clear();
        }
        G.stream_finish_construction();
}


void graph_io_stream::recoverBlockAssignedToNode(PartitionConfig & config, graph_access & G, NodeID node, NodeID node_counter) {
	if (node >= node_counter - config.quotient_nodes) { // quotient nodes
		NodeID targetPar = node - (node_counter - config.quotient_nodes);
		G.setPartitionIndex(node, targetPar);
	} else if (config.restream_number && (config.restream_vcycle || config.initial_partitioning_type == INITIAL_PARTITIONING_FENNEL) ) { 
		if (node < config.nmbNodes) { // regular nodes
			LongNodeID global_node = config.lower_global_node + (LongNodeID) node - 1;
			PartitionID targetPar = (*config.stream_nodes_assign)[global_node];
			G.setPartitionIndex(node, targetPar);
		} else { // ghost node. PS: there are no isolated ghost nodes in restream
			std::cerr << "Unexpected branch.\n";
			exit(0);
		}
	}
}



void graph_io_stream::generalizeStreamPartition(PartitionConfig & config, graph_access & G_local) {
	for(NodeID node = 0, end = config.nmbNodes; node < end; node++) {
		PartitionID block = G_local.getPartitionIndex(node);
		LongNodeID global_node = (LongNodeID) node + config.lower_global_node - 1;
		(*config.stream_nodes_assign)[global_node] = block;
		(*config.stream_blocks_weight)[block] += G_local.getNodeWeight(node) - G_local.getImplicitGhostNodes(node);
	} 
}

void graph_io_stream::countAssignedNodes(PartitionConfig & config) {
	config.stream_assigned_nodes = 0;
	for(int i=0; i<(int)config.stream_blocks_weight->size(); i++) { 
		config.stream_assigned_nodes += (*config.stream_blocks_weight)[i];
	}
}

// This method modifies the matrix edges_virtualReal 
void graph_io_stream::onePassPartition(PartitionConfig & config, std::vector<std::vector<EdgeWeight>> & edges_virtualReal,
				std::vector<PartitionID> & blockVirtualToReal, std::vector<NodeWeight> & weight_VirtualBlocks) {
	int len = edges_virtualReal.size();
	for (int i=0; i<len; i++) {
		PartitionID decision = onePassDecide(config, i, edges_virtualReal[i]);
		blockVirtualToReal[i] = decision;
		(*config.stream_blocks_weight)[decision] += weight_VirtualBlocks[i];
	}
}

int graph_io_stream::onePassDecide(PartitionConfig & config, NodeID node, std::vector<EdgeWeight> & edges_i_real) {
	PartitionID decision;
	double best = std::numeric_limits<double>::lowest();
	LongNodeID global_node = (LongNodeID) node + config.lower_global_node - 1;
	double score = 0;
	double fennel_weight = getFennelWeight(config);
	int blocks = config.k;
	EdgeWeight block_weight;
	switch(config.one_pass_algorithm) {
		case ONEPASS_HASHING:
			decision = fnv1a(global_node) % config.k;
			break;
		case ONEPASS_GREEDY:
			for (int j=0; j<blocks; j++) {
				if ((double)edges_i_real[j] > best) {
					decision = j;
					best =(double) edges_i_real[j];
				}
			}
			break;
		case ONEPASS_LDG:
			for (int j=0; j<blocks; j++) {
				block_weight = (*config.stream_blocks_weight)[j] + (*config.add_blocks_weight)[j];
				if (block_weight  >= (EdgeWeight)config.stream_total_upperbound) {
					continue;
				}
				score = (0.1+edges_i_real[j])*(config.stream_total_upperbound - block_weight)/(double)config.stream_total_upperbound;
				if (score > best) {
					decision = j;
					best = score;
				}
			}
			break;
		case ONEPASS_FENNEL:
			for (int j=0; j<blocks; j++) {
				block_weight = (*config.stream_blocks_weight)[j] + (*config.add_blocks_weight)[j];
				if (block_weight  >= (EdgeWeight)config.stream_total_upperbound) {
					continue;
				}
				score = (0.1+edges_i_real[j])-fennel_weight*(config.fennel_alpha_gamma*	std::pow(block_weight, config.fennel_gamma-1) );
				if (score > best) {
					decision = j;
					best = score;
				}
			}
			break;
		case ONEPASS_CHUNKING:
			decision = global_node % config.k;
			break;
		case ONEPASS_FRACTIONAL_GREEDY:
			for (int j=0; j<blocks; j++) {
				score = (0.1+edges_i_real[j]) - (double)config.stream_total_upperbound / (config.stream_total_upperbound -
					(*config.stream_blocks_weight)[j] - (*config.add_blocks_weight)[j]);
				if (score > best) {
					decision = j;
					best = score;
				}
			}
			break;
	}
	return decision;
}


double graph_io_stream::getFennelWeight(PartitionConfig & partition_config) {
	double fennel_weight = 2;
	double fennel_tmp = 0;
	switch(partition_config.fennel_dynamics) {
		case FENNELADP_ORIGINAL:
			fennel_weight = 1;
			break;
		case FENNELADP_DOUBLE:
			fennel_weight = 2;
			break;
		case FENNELADP_LINEAR:
			fennel_weight = 2*partition_config.remaining_stream_nodes/(double)partition_config.stream_nodes_assign->size();
			break;
		case FENNELADP_MID_LINEAR:
			fennel_tmp = 2*partition_config.remaining_stream_nodes/(double)partition_config.stream_nodes_assign->size();
			if (fennel_tmp <= 1) {
				fennel_weight = 2 * (fennel_tmp);
			}
			break;
		case FENNELADP_QUADRATIC:
			fennel_tmp = partition_config.remaining_stream_nodes/(double)partition_config.stream_nodes_assign->size();
			fennel_weight = 2*fennel_tmp*fennel_tmp;
			break;
		case FENNELADP_MID_QUADRATIC:
			fennel_tmp = 2*partition_config.remaining_stream_nodes/(double)partition_config.stream_nodes_assign->size();
			if (fennel_tmp <= 1) {
				fennel_weight = 2*fennel_tmp*fennel_tmp;
			}
			break;
		case FENNELADP_MID_CONSTANT:
			fennel_tmp = partition_config.remaining_stream_nodes/(double)partition_config.stream_nodes_assign->size();
			if (fennel_tmp <= 1.5) {
				fennel_weight = 0.5;
			}
			break;
		case FENNELADP_EDGE_CUT:
			fennel_weight = 0;
			break;
	}
	return fennel_weight;
}


void graph_io_stream::writePartitionStream(PartitionConfig & config, const std::string & filename) {
        std::ofstream f(filename.c_str());
        std::cout << "writing partition to " << filename << " ... " << std::endl;
	
	for (LongNodeID node = 0; node < config.stream_nodes_assign->size(); node++) {
                f << (*config.stream_nodes_assign)[node] <<  "\n";
        } 

        f.close();
}


void graph_io_stream::readFirstLineStream(PartitionConfig & partition_config, std::string graph_filename, EdgeWeight& total_edge_cut) {
	if (partition_config.stream_in != NULL) {
		delete partition_config.stream_in;
	}
	partition_config.stream_in = new std::ifstream(graph_filename.c_str());
	if (!(*(partition_config.stream_in))) {
		std::cerr << "Error opening " << graph_filename << std::endl;
		exit(1);
	}
	std::vector<std::string>* lines;

	lines = new std::vector<std::string>(1);
	std::getline(*(partition_config.stream_in),(*lines)[0]);

	//skip comments
	while( (*lines)[0][0] == '%' ) {
		std::getline(*(partition_config.stream_in),(*lines)[0]);
	}

	std::stringstream ss((*lines)[0]);
	ss >> partition_config.remaining_stream_nodes;
	ss >> partition_config.remaining_stream_edges;
	ss >> partition_config.remaining_stream_ew;

	if (partition_config.stream_nodes_assign == NULL) {
		partition_config.stream_nodes_assign  = new std::vector<PartitionID>(partition_config.remaining_stream_nodes, INVALID_PARTITION);
	}
	if (partition_config.stream_blocks_weight == NULL) {
		partition_config.stream_blocks_weight = new std::vector<NodeWeight>(partition_config.k, 0);
	}
	if (partition_config.add_blocks_weight == NULL) {
		partition_config.add_blocks_weight = new std::vector<NodeWeight>(partition_config.k, 0);
	}
	partition_config.total_stream_nodeweight = 0;
	partition_config.total_stream_nodecounter = 0;
	partition_config.stream_n_nodes = partition_config.remaining_stream_nodes;

	if (partition_config.num_streams_passes > 1 + partition_config.restream_number) {
		partition_config.stream_total_upperbound = ceil(((100+1.5*partition_config.imbalance)/100.)*
					(partition_config.remaining_stream_nodes/(double)partition_config.k));
	} else {
		partition_config.stream_total_upperbound = ceil(((100+partition_config.imbalance)/100.)*
					(partition_config.remaining_stream_nodes/(double)partition_config.k));
	}

	partition_config.fennel_alpha = partition_config.remaining_stream_edges * 
				std::pow(partition_config.k,partition_config.fennel_gamma-1) / 
				(std::pow(partition_config.remaining_stream_nodes,partition_config.fennel_gamma));

        partition_config.fennel_alpha_gamma = partition_config.fennel_alpha * partition_config.fennel_gamma;

	partition_config.quotient_nodes = partition_config.k;
	 
	total_edge_cut = 0;
	if (partition_config.stream_buffer_len == 0) { // signal of partial restream standard buffer size
		partition_config.stream_buffer_len = (LongNodeID) ceil(partition_config.remaining_stream_nodes/(double)partition_config.k);
	}
	partition_config.nmbNodes = MIN(partition_config.stream_buffer_len, partition_config.remaining_stream_nodes);
	partition_config.n_batches = ceil(partition_config.remaining_stream_nodes / (double)partition_config.nmbNodes);
	partition_config.curr_batch = 0;
//	partition_config.stream_global_epsilon = (partition_config.imbalance)/100.;

	delete lines;
}



void graph_io_stream::loadRemainingLines(PartitionConfig & partition_config, LINE_BUFFER &lines) {
	if (partition_config.ram_stream) {
		lines = graph_io_stream::loadLinesFromStream(partition_config, partition_config.remaining_stream_nodes);
	}
}


void graph_io_stream::loadBufferLines(PartitionConfig & partition_config, LINE_BUFFER &lines, LongNodeID num_lines) {
	if (!partition_config.ram_stream) {
		lines = graph_io_stream::loadLinesFromStream(partition_config, num_lines);
	}
}


std::vector<std::string>* graph_io_stream::loadLinesFromStream(PartitionConfig & partition_config, LongNodeID num_lines) {
	std::vector<std::string>* lines;
	lines = new std::vector<std::string>(num_lines);
	LongNodeID node_counter = 0;
	while( node_counter < num_lines) {
		std::getline(*(partition_config.stream_in),(*lines)[node_counter]);
		if ((*lines)[node_counter][0] == '%') { // a comment in the file
			continue;
		}
		node_counter++;
	}
	return lines;
}

void graph_io_stream::prescribeBufferInbalance(PartitionConfig & partition_config) {
	double &global_epsilon = partition_config.stream_global_epsilon;
	int &passes = partition_config.num_streams_passes;
	if (partition_config.restream_number && partition_config.restream_number >= partition_config.num_streams_passes -1) {
		partition_config.imbalance = 100*global_epsilon;
	} else {
		double current_nodes = (double)partition_config.stream_assigned_nodes + partition_config.nmbNodes + partition_config.ghost_nodes;
		partition_config.imbalance = 100 * partition_config.stream_n_nodes * (1+global_epsilon) / current_nodes - 100;
		if (passes > 1) {
			partition_config.imbalance = MIN(partition_config.imbalance,partition_config.batch_inbalance);
		} else {
			partition_config.imbalance = MIN(MAX(100*global_epsilon,0.75*partition_config.imbalance),partition_config.batch_inbalance);
		}
	}
}

void graph_io_stream::streamEvaluatePartition(PartitionConfig & config, const std::string & filename, EdgeWeight& edgeCut) {
	std::vector<std::vector<LongNodeID>>* input;
	std::vector<std::string>* lines;
	lines = new std::vector<std::string>(1);
	LongNodeID node_counter = 0;
	buffered_input *ss2 = NULL;
	std::string line;
	std::ifstream in(filename.c_str());
	if (!in) {
		std::cerr << "Error opening " << filename << std::endl;
		return 1;
	}
	long nmbNodes;
	long nmbEdges;
	int ew = 0;
	std::getline(in,(*lines)[0]);
	while ((*lines)[0][0] == '%') std::getline(in,(*lines)[0]); // a comment in the file

	std::stringstream ss((*lines)[0]);
	ss >> nmbNodes;
	ss >> nmbEdges;
	ss >> ew;
	bool read_ew = false;
	bool read_nw = false;
	if(ew == 1) {
		read_ew = true;
	} else if (ew == 11) {
		read_ew = true;
		read_nw = true;
	} else if (ew == 10) {
		read_nw = true;
	}
	NodeID target;
	NodeWeight total_nodeweight = 0;
	EdgeWeight total_edgeweight = 0;
	edgeCut = 0;

	while(  std::getline(in, (*lines)[0])) {
		if ((*lines)[0][0] == '%') continue; // a comment in the file
		NodeID node = node_counter++;
		PartitionID partitionIDSource = (*config.stream_nodes_assign)[node];

		input = new std::vector<std::vector<LongNodeID>>(1);
		ss2 = new buffered_input(lines);
		ss2->simple_scan_line((*input)[0]);
		std::vector<LongNodeID> &line_numbers = (*input)[0];
		LongNodeID col_counter = 0;

		NodeWeight weight = 1;
		if( read_nw ) {
			weight = line_numbers[col_counter++];
			total_nodeweight += weight;
		}
		while (col_counter < line_numbers.size()) {
			target = line_numbers[col_counter++];
			target = target-1;
			EdgeWeight edge_weight = 1;
			if( read_ew ) {
				edge_weight = line_numbers[col_counter++];
			}
			total_edgeweight += edge_weight;
			PartitionID partitionIDTarget = (*config.stream_nodes_assign)[target];
			if (partitionIDSource != partitionIDTarget) {
				edgeCut += edge_weight;
			}
		}
		(*lines)[0].clear(); delete ss2;
		delete input;
		if(in.eof()) {
			break;
		}
	}
	edgeCut = edgeCut/2; // Since every edge is counted twice
	delete lines;
}


template< typename T>
T graph_io_stream::return_and_delete_element(std::vector<T> & vec, LongNodeID pos) {
	T val = vec[pos];
	vec[pos] = vec[vec.size()-1];
	vec.erase(vec.begin()+vec.size()-1); 
	return val;
}


