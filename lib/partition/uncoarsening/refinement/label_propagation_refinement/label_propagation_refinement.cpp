/******************************************************************************
 * label_propagation_refinement.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#include "label_propagation_refinement.h"
#include "partition/coarsening/clustering/node_ordering.h"
#include "tools/random_functions.h"
#include <map>
#include <set>



label_propagation_refinement::label_propagation_refinement() {
                
}

label_propagation_refinement::~label_propagation_refinement() {
                
}

EdgeWeight label_propagation_refinement::perform_refinement(PartitionConfig & partition_config, 
                                                            graph_access & G, 
                                                            complete_boundary & boundary) {
	return (EdgeWeight) 0;
}

EdgeWeight label_propagation_refinement::perform_refinement_fennel(PartitionConfig & partition_config, graph_access & G, 
										    complete_boundary & boundary) {
	random_functions::fastRandBool<uint64_t> random_obj;
        std::vector<PartitionID> hash_map(partition_config.k,0);
//        std::vector<NodeID> permutation(G.number_of_nodes());
        std::vector<NodeWeight> cluster_sizes(partition_config.k, 0);
        std::vector<NodeWeight> cluster_ghost_nodes(partition_config.k, 0);

        node_ordering n_ordering;
//        n_ordering.order_nodes(partition_config, G, permutation);

        std::queue< NodeID > * Q             = new std::queue< NodeID >();
        std::queue< NodeID > * next_Q        = new std::queue< NodeID >();
        std::vector<bool> * Q_contained      = new std::vector<bool>(G.number_of_nodes(), false);
        std::vector<bool> * next_Q_contained = new std::vector<bool> (G.number_of_nodes(), false);
        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
		cluster_ghost_nodes[G.getPartitionIndex(node)] += G.getImplicitGhostNodes(node);
//                Q->push(permutation[node]);
                Q->push(node);
        } endfor

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
			fennel_weight = 2*partition_config.remaining_stream_nodes/(double)partition_config.total_stream_edges;
			break;
		case FENNELADP_MID_LINEAR:
			fennel_tmp = 2*partition_config.remaining_stream_nodes/(double)partition_config.total_stream_edges;
			if (fennel_tmp <= 1) {
				fennel_weight = 2 * (fennel_tmp);
			}
			break;
		case FENNELADP_QUADRATIC:
			fennel_tmp = partition_config.remaining_stream_nodes/(double)partition_config.total_stream_edges;
			fennel_weight = 2*fennel_tmp*fennel_tmp;
			break;
		case FENNELADP_MID_QUADRATIC:
			fennel_tmp = 2*partition_config.remaining_stream_nodes/(double)partition_config.total_stream_edges;
			if (fennel_tmp <= 1) {
				fennel_weight = 2*fennel_tmp*fennel_tmp;
			}
			break;
		case FENNELADP_MID_CONSTANT:
			fennel_tmp = partition_config.remaining_stream_nodes/(double)partition_config.total_stream_edges;
			if (fennel_tmp <= 1.5) {
				fennel_weight = 0.5;
			}
			break;
		case FENNELADP_EDGE_CUT:
			fennel_weight = 0;
			break;
	}

	if (partition_config.remaining_stream_nodes == 0 && partition_config.fennel_dynamics != FENNELADP_ORIGINAL) {
	//if (partition_config.remaining_stream_nodes == 0) {
		fennel_weight = 0;
	}

        for( int j = 0; j < partition_config.label_iterations_refinement; j++) {
                unsigned int change_counter = 0;
                while( !Q->empty() ) {
                        NodeID node = Q->front();
                        Q->pop();
                        (*Q_contained)[node] = false;

			if (node >= G.number_of_nodes() - partition_config.quotient_nodes) {
				continue;
			}

                        //now move the node to the cluster that is most common in the neighborhood
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                hash_map[G.getPartitionIndex(target)]+=G.getEdgeWeight(e);
                        } endfor

                        //second sweep for finding max and resetting array
                        PartitionID max_block = G.getPartitionIndex(node);
                        PartitionID my_block  = G.getPartitionIndex(node);

//                        double max_value = std::numeric_limits<double>::lowest();
			double max_value = hash_map[my_block];
			max_value -= fennel_weight*(G.getNodeWeight(node)*partition_config.fennel_alpha_gamma*random_functions::approx_sqrt(cluster_sizes[my_block] -G.getNodeWeight(node)) ); 
			hash_map[my_block] = 0;

			double cur_value=0;
			NodeID target;
			PartitionID cur_block;
                        forall_out_edges(G, e, node) {
                                target	   = G.getEdgeTarget(e);
                                cur_block  = G.getPartitionIndex(target);
				if (hash_map[cur_block] == 0) {
					continue;
				}
                                cur_value  = hash_map[cur_block];
				cur_value -= fennel_weight * (G.getNodeWeight(node) * partition_config.fennel_alpha_gamma * random_functions::approx_sqrt(cluster_sizes[cur_block] - (my_block==cur_block)*G.getNodeWeight(node)) ); 
                                if((cur_value > max_value  || (cur_value == max_value && random_obj.nextBool()))
				    && (cluster_sizes[cur_block] -cluster_ghost_nodes[cur_block] + G.getNodeWeight(node) -G.getImplicitGhostNodes(node) <= partition_config.stream_total_upperbound
				    || cur_block == my_block
				    ||  cluster_sizes[cur_block] -cluster_ghost_nodes[cur_block] + G.getNodeWeight(node) -G.getImplicitGhostNodes(node) <= (cluster_sizes[my_block] - cluster_ghost_nodes[my_block] ) ) ) { // if solution already violates imbalance
                                        max_value = cur_value;
                                        max_block = cur_block;
                                }

                                hash_map[cur_block] = 0;
                        } endfor

                        cluster_sizes[G.getPartitionIndex(node)]	-= G.getNodeWeight(node);
			cluster_ghost_nodes[G.getPartitionIndex(node)]	-= G.getImplicitGhostNodes(node);
                        cluster_sizes[max_block]		+= G.getNodeWeight(node);
			cluster_ghost_nodes[max_block]		+= G.getImplicitGhostNodes(node);
                        bool changed_label                = G.getPartitionIndex(node) != max_block; 
                        change_counter                   += changed_label;
                        G.setPartitionIndex(node, max_block);


                        if(changed_label) {
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        if(!(*next_Q_contained)[target]) {
                                                next_Q->push(target);
                                                (*next_Q_contained)[target] = true;
                                        } 
                                } endfor

                                G.setPartitionIndex(node, max_block);        
                        }
                } 

                std::swap( Q, next_Q);
                std::swap( Q_contained, next_Q_contained);

        }
        

        delete Q;
        delete next_Q;
        delete Q_contained;
        delete next_Q_contained;
        
        return 0;

}

bool label_propagation_refinement::is_boundary(NodeID node, graph_access & G) {

        PartitionID my_part   = G.getPartitionIndex(node);

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if (G.getPartitionIndex(target) == my_part) {
                        return true;
                }
        } endfor
        return false;
}  
