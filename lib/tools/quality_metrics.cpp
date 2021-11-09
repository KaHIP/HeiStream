/******************************************************************************
 * quality_metrics.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include <cmath>

#include "quality_metrics.h"
#include "data_structure/union_find.h"
#include "random_functions.h"

#include <unordered_map>
#include <vector>
#include <list>
#include <numeric>

quality_metrics::quality_metrics() {
}

quality_metrics::~quality_metrics () {
}

float quality_metrics::getFennelWeight(PartitionConfig & partition_config) {
	float fennel_weight = 2;
	float fennel_tmp = 0;
	switch(partition_config.fennel_dynamics) {
		case FENNELADP_ORIGINAL:
			fennel_weight = 1;
			break;
		case FENNELADP_DOUBLE:
			fennel_weight = 2;
			break;
		case FENNELADP_LINEAR:
			fennel_weight = 2*partition_config.remaining_stream_nodes/(float)partition_config.stream_nodes_assign->size();
			break;
		case FENNELADP_MID_LINEAR:
			fennel_tmp = 2*partition_config.remaining_stream_nodes/(float)partition_config.stream_nodes_assign->size();
			if (fennel_tmp <= 1) {
				fennel_weight = 2 * (fennel_tmp);
			}
			break;
		case FENNELADP_QUADRATIC:
			fennel_tmp = partition_config.remaining_stream_nodes/(float)partition_config.stream_nodes_assign->size();
			fennel_weight = 2*fennel_tmp*fennel_tmp;
			break;
		case FENNELADP_MID_QUADRATIC:
			fennel_tmp = 2*partition_config.remaining_stream_nodes/(float)partition_config.stream_nodes_assign->size();
			if (fennel_tmp <= 1) {
				fennel_weight = 2*fennel_tmp*fennel_tmp;
			}
			break;
		case FENNELADP_MID_CONSTANT:
			fennel_tmp = partition_config.remaining_stream_nodes/(float)partition_config.stream_nodes_assign->size();
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

EdgeWeight quality_metrics::fennel_objective(PartitionConfig & partition_config, graph_access & G, int * partition_map, double fennel_gamma, double fennel_alpha) {
        EdgeWeight interiorEdges = 0;
	int k = G.get_partition_count();
	std::vector<int> weights(k,0);
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = partition_map[n];
		weights[partitionIDSource] = G.getNodeWeight(n);
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = partition_map[targetNode];

                        if (partitionIDSource == partitionIDTarget) {
                                interiorEdges += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
	double alpha_gamma = fennel_gamma * fennel_alpha;
        double objective = 0;
        for(int i = 0; i < k; i++) {
		objective += alpha_gamma * weights[i] * random_functions::approx_sqrt(weights[i]); // exponent is gamma +1 because each node contributes
        }
        objective *= getFennelWeight(partition_config);
        objective -= interiorEdges;
	return objective*0.5;	// because edges and nonneighbors are computed twice
}

EdgeWeight quality_metrics::fennel_objective(PartitionConfig & partition_config, graph_access & G, double fennel_gamma, double fennel_alpha) {
        EdgeWeight interiorEdges = 0;
	int k = G.get_partition_count();
	std::vector<int> weights(k,0);
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = G.getPartitionIndex(n);
		weights[partitionIDSource] = G.getNodeWeight(n);
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                        if (partitionIDSource == partitionIDTarget) {
                                interiorEdges += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
	double alpha_gamma = fennel_gamma * fennel_alpha;
        double objective = 0;
        for(int i = 0; i < k; i++) {
		objective += alpha_gamma * weights[i] * random_functions::approx_sqrt(weights[i]); // exponent is gamma +1 because each node contributes
        }
        objective *= getFennelWeight(partition_config);
        objective -= interiorEdges;
	return objective*0.5;	// because edges and nonneighbors are computed twice
}

EdgeWeight quality_metrics::ghost_edge_cut(const PartitionConfig & config, graph_access & G) {
        EdgeWeight edgeCut = 0;
	for (LongNodeID ghostkey=0; ghostkey < config.ghostkey_to_edges->size(); ghostkey++) {
		NodeID node = (*config.ghostkey_to_node)[ghostkey];
		auto& node_list = (*config.ghostkey_to_edges)[ghostkey];
		for (auto& edge : node_list) {
			if(G.getPartitionIndex(node) != G.getPartitionIndex(edge.first)) {
                                edgeCut += edge.second;
			}
		}
	}
        return edgeCut;
}

EdgeWeight quality_metrics::edge_cut(graph_access & G) {
        EdgeWeight edgeCut = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = G.getPartitionIndex(n);
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                        if (partitionIDSource != partitionIDTarget) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
        return edgeCut/2;
}

EdgeWeight quality_metrics::edge_cut_full_stream(const PartitionConfig & config, graph_access & G, 
						std::vector<std::vector<EdgeWeight>> & edges_virtualReal) {
	PartitionID partitionIDSource = 0;
        EdgeWeight edgeCut = 0;
	int blocks = edges_virtualReal[0].size();
        forall_nodes(G, n) { 
		if (n >= config.nmbNodes) {
			break;
		}
                partitionIDSource = G.getPartitionIndex(n);
		for (PartitionID i=0; i<blocks; i++) {
			if (i != partitionIDSource) {
				edgeCut += edges_virtualReal[n][i];
			}
		}
        } endfor
        return edgeCut;
}

EdgeWeight quality_metrics::edge_cut(graph_access & G, int * partition_map) {
        EdgeWeight edgeCut = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = partition_map[n];
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = partition_map[targetNode];

                        if (partitionIDSource != partitionIDTarget) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
        return edgeCut/2;
}

EdgeWeight quality_metrics::edge_cut(graph_access & G, PartitionID lhs, PartitionID rhs) {
        EdgeWeight edgeCut = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = G.getPartitionIndex(n);
                if(partitionIDSource != lhs) continue;
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                        if(partitionIDTarget == rhs) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
        return edgeCut;
}

EdgeWeight quality_metrics::edge_cut_connected(graph_access & G, int * partition_map) {
        EdgeWeight edgeCut = 0;
        EdgeWeight sumEW = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = partition_map[n];
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = partition_map[targetNode];

                        if (partitionIDSource != partitionIDTarget) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                        sumEW+=G.getEdgeWeight(e);
                } endfor 
        } endfor
        union_find uf(G.number_of_nodes());
        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] == partition_map[target]) {
                                uf.Union(node, target); 
                        }
                } endfor
        } endfor

        std::unordered_map<NodeID, NodeID> size_right;
        forall_nodes(G, node) {
                size_right[uf.Find(node)] = 1;
        } endfor


        std::cout <<  "number of connected comp " <<  size_right.size()  << std::endl;
        if( size_right.size() == G.get_partition_count()) {
                return edgeCut/2;
        } else {
                return edgeCut/2+sumEW*size_right.size();
        }

}


EdgeWeight quality_metrics::max_communication_volume(graph_access & G, int * partition_map) {
    std::vector<EdgeWeight> block_volume(G.get_partition_count(),0);
    forall_nodes(G, node) {
        PartitionID block = partition_map[node];
        std::vector<bool> block_incident(G.get_partition_count(), false);
        block_incident[block] = true;

        int num_incident_blocks = 0;

        forall_out_edges(G, e, node) {
            NodeID target = G.getEdgeTarget(e);
            PartitionID target_block = partition_map[target];
            if(!block_incident[target_block]) {
                block_incident[target_block] = true;
                num_incident_blocks++;
            }
        } endfor
        block_volume[block] += num_incident_blocks;
    } endfor

    EdgeWeight max_comm_volume = *(std::max_element(block_volume.begin(), block_volume.end()));
    return max_comm_volume;
}

EdgeWeight quality_metrics::min_communication_volume(graph_access & G) {
    std::vector<EdgeWeight> block_volume(G.get_partition_count(),0);
    forall_nodes(G, node) {
        PartitionID block = G.getPartitionIndex(node);
        std::vector<bool> block_incident(G.get_partition_count(), false);
        block_incident[block] = true;
        int num_incident_blocks = 0;

        forall_out_edges(G, e, node) {
            NodeID target = G.getEdgeTarget(e);
            PartitionID target_block = G.getPartitionIndex(target);
            if(!block_incident[target_block]) {
                block_incident[target_block] = true;
                num_incident_blocks++;
            }
        } endfor
        block_volume[block] += num_incident_blocks;
    } endfor

    EdgeWeight min_comm_volume = *(std::min_element(block_volume.begin(), block_volume.end()));
    return min_comm_volume;
}

EdgeWeight quality_metrics::max_communication_volume(graph_access & G) {
    std::vector<EdgeWeight> block_volume(G.get_partition_count(),0);
    forall_nodes(G, node) {
        PartitionID block = G.getPartitionIndex(node);
        std::vector<bool> block_incident(G.get_partition_count(), false);
        block_incident[block] = true;
        int num_incident_blocks = 0;

        forall_out_edges(G, e, node) {
            NodeID target = G.getEdgeTarget(e);
            PartitionID target_block = G.getPartitionIndex(target);
            if(!block_incident[target_block]) {
                block_incident[target_block] = true;
                num_incident_blocks++;
            }
        } endfor
        block_volume[block] += num_incident_blocks;
    } endfor

    EdgeWeight max_comm_volume = *(std::max_element(block_volume.begin(), block_volume.end()));
    return max_comm_volume;
}

EdgeWeight quality_metrics::total_communication_volume(graph_access & G) {
    std::vector<EdgeWeight> block_volume(G.get_partition_count(),0);
    forall_nodes(G, node) {
        PartitionID block = G.getPartitionIndex(node);
        std::vector<bool> block_incident(G.get_partition_count(), false);
        block_incident[block] = true;
        int num_incident_blocks = 0;

        forall_out_edges(G, e, node) {
            NodeID target = G.getEdgeTarget(e);
            PartitionID target_block = G.getPartitionIndex(target);
            if(!block_incident[target_block]) {
                block_incident[target_block] = true;
                num_incident_blocks++;
            }
        } endfor
        block_volume[block] += num_incident_blocks;
    } endfor

    EdgeWeight total_comm_volume = std::accumulate(block_volume.begin(), block_volume.end(),0);
    return total_comm_volume;
}



int quality_metrics::boundary_nodes(graph_access& G) {
        int no_of_boundary_nodes = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = G.getPartitionIndex(n);

                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                        if (partitionIDSource != partitionIDTarget) {
                                no_of_boundary_nodes++;
                                break; 
                        }
                } endfor 
        }       endfor
        return no_of_boundary_nodes;
}

double quality_metrics::balance_separator(graph_access& G) {
        std::vector<PartitionID> part_weights(G.get_partition_count(), 0);

        double overallWeight = 0;

        forall_nodes(G, n) {
                PartitionID curPartition = G.getPartitionIndex(n);
                part_weights[curPartition] += G.getNodeWeight(n);
                overallWeight += G.getNodeWeight(n);
        } endfor

        double balance_part_weight = ceil(overallWeight / (double)(G.get_partition_count()-1));
        double cur_max             = -1;

        PartitionID separator_block = G.getSeparatorBlock();
        forall_blocks(G, p) {
                if( p == separator_block ) continue;
                double cur = part_weights[p];
                if (cur > cur_max) {
                        cur_max = cur;
                }
        } endfor

        double percentage = cur_max/balance_part_weight;
        return percentage;
}

NodeWeight quality_metrics::separator_weight(graph_access& G) {
        NodeWeight separator_size = 0;
        PartitionID separator_ID = G.getSeparatorBlock();
        forall_nodes(G, node) {
                if( G.getPartitionIndex(node) == separator_ID) {
                        separator_size += G.getNodeWeight(node);
		}
        } endfor

        return separator_size;
}



double quality_metrics::balance_full_stream(std::vector<NodeWeight> &stream_blocks_weight) {
	double total_weight = 0;
	double max_block_weight = 0;
	for (int i=0; i<(int)stream_blocks_weight.size(); i++) {
		total_weight += (double) stream_blocks_weight[i];
		if (stream_blocks_weight[i] > max_block_weight) {
			max_block_weight = (double) stream_blocks_weight[i];
		}
	}
        double balance_part_weight = ceil(total_weight / (double)stream_blocks_weight.size());
        double percentage = max_block_weight/balance_part_weight;
        return percentage;
}


double quality_metrics::balance(graph_access& G) {
        std::vector<PartitionID> part_weights(G.get_partition_count(), 0);

        double overallWeight = 0;

        forall_nodes(G, n) {
                PartitionID curPartition = G.getPartitionIndex(n);
                part_weights[curPartition] += G.getNodeWeight(n);
                overallWeight += G.getNodeWeight(n);
        } endfor

        double balance_part_weight = ceil(overallWeight / (double)G.get_partition_count());
        double cur_max             = -1;

        forall_blocks(G, p) {
                double cur = part_weights[p];
                if (cur > cur_max) {
                        cur_max = cur;
                }
        } endfor

        double percentage = cur_max/balance_part_weight;
        return percentage;
}

double quality_metrics::edge_balance(graph_access &G, const std::vector<PartitionID> &edge_partition) {
    std::vector<PartitionID> part_weights(G.get_partition_count(), 0);

    double overallWeight = 0;

    forall_edges(G, e) {
        PartitionID curPartition = edge_partition[e];
        ++part_weights[curPartition];
        ++overallWeight;
    } endfor

    double balance_part_weight = ceil(overallWeight / (double)G.get_partition_count());
    double cur_max             = -1;

    forall_blocks(G, p) {
        double cur = part_weights[p];
        if (cur > cur_max) {
            cur_max = cur;
        }
    } endfor

    double percentage = cur_max/balance_part_weight;
    return percentage;
}

double quality_metrics::balance_edges(graph_access& G) {
        std::vector<PartitionID> part_weights(G.get_partition_count(), 0);

        double overallWeight = 0;

        forall_nodes(G, n) {
                PartitionID curPartition = G.getPartitionIndex(n);
                part_weights[curPartition] += G.getNodeDegree(n);
                overallWeight += G.getNodeDegree(n);
        } endfor

        double balance_part_weight = ceil(overallWeight / (double)G.get_partition_count());
        double cur_max             = -1;

        forall_blocks(G, p) {
                double cur = part_weights[p];
                if (cur > cur_max) {
                        cur_max = cur;
                }
        } endfor

        double percentage = cur_max/balance_part_weight;
        return percentage;
}

EdgeWeight quality_metrics::objective(const PartitionConfig & config, graph_access & G, int* partition_map) {
        if(config.mh_optimize_communication_volume) {
                return max_communication_volume(G, partition_map);
        } else if(config.mh_penalty_for_unconnected) {
                return edge_cut_connected(G, partition_map);
        } else {
                return edge_cut(G, partition_map);
        }
}

