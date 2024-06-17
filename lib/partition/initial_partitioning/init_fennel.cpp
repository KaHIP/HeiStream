/******************************************************************************
 * init_fennel.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "init_fennel.h"
#include "data_structure/priority_queues/maxNodeHeap.h"


init_fennel::init_fennel() {

}

init_fennel::~init_fennel() {

}

void init_fennel::initial_partition(PartitionConfig &config,
                                    const unsigned int seed,
                                    graph_access &G,
                                    int *partition_map,
                                    int ismultisec/*=0*/) {

    timer t;
    t.restart();
    unsigned iterations = 1;
    EdgeWeight best_cut = std::numeric_limits<EdgeWeight>::max();

    for (unsigned i = 0; i < iterations; i++) {
        fennel(config, G);
        G.set_partition_count(config.quotient_nodes);

        quality_metrics qm;
        EdgeWeight curcut = 0;

        if (config.use_fennel_objective) {
            curcut = qm.fennel_objective(config, G, config.fennel_gamma, config.fennel_alpha);
        } else {
            curcut = qm.edge_cut(G);
        }

        if (curcut < best_cut) {
            best_cut = curcut;
            forall_nodes(G, n)
                    {
                        partition_map[n] = G.getPartitionIndex(n);
                    }
            endfor
        }
    }
    PRINT(std::cout << "init_fennel took " << t.elapsed() << std::endl;)
}

void init_fennel::initial_partition(PartitionConfig &config,
                                    const unsigned int seed,
                                    graph_access &G,
                                    int *xadj,
                                    int *adjncy,
                                    int *vwgt,
                                    int *adjwgt,
                                    int *partition_map,
                                    int ismultisec/*=0*/) {

    std::cout << "not implemented yet" << std::endl;

}


EdgeWeight init_fennel::fennel(PartitionConfig &partition_config, graph_access &G) {
    random_functions::fastRandBool<uint64_t> random_obj;
    double cur_value = 0;
    bool node_too_large = false;
    std::vector <PartitionID> hash_map(partition_config.k, 0);
    std::vector <NodeWeight> cluster_sizes(partition_config.k, 0);
    std::vector <NodeWeight> cluster_ghost_nodes(partition_config.k, 0);
    maxNodeHeap *queue = new maxNodeHeap();

    if (partition_config.restream_number) {
        for (NodeID node = 0; node < G.number_of_nodes() - partition_config.quotient_nodes; node++) {
            cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
        }
    }
    for (NodeID node = G.number_of_nodes() - partition_config.quotient_nodes, end = G.number_of_nodes();
         node < end; node++) {
        cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
        if (partition_config.use_queue) {
            queue->insert(G.getPartitionIndex(node), -G.getNodeWeight(node));
        }
        cluster_ghost_nodes[G.getPartitionIndex(node)] += G.getImplicitGhostNodes(node);
    }

    double fennel_weight = 2;
    double fennel_tmp = 0;
    switch (partition_config.fennel_dynamics) {
        case FENNELADP_ORIGINAL:
            fennel_weight = 1;
            break;
        case FENNELADP_DOUBLE:
            fennel_weight = 2;
            break;
        case FENNELADP_LINEAR:
            fennel_weight = 2 * partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            break;
        case FENNELADP_MID_LINEAR:
            fennel_tmp = 2 * partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            if (fennel_tmp <= 1) {
                fennel_weight = 2 * (fennel_tmp);
            }
            break;
        case FENNELADP_QUADRATIC:
            fennel_tmp = partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            fennel_weight = 2 * fennel_tmp * fennel_tmp;
            break;
        case FENNELADP_MID_QUADRATIC:
            fennel_tmp = 2 * partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            if (fennel_tmp <= 1) {
                fennel_weight = 2 * fennel_tmp * fennel_tmp;
            }
            break;
        case FENNELADP_MID_CONSTANT:
            fennel_tmp = partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            if (fennel_tmp <= 1.5) {
                fennel_weight = 0.5;
            }
            break;
        case FENNELADP_EDGE_CUT:
            fennel_weight = 0;
            break;
    }

//	if (partition_config.remaining_stream_nodes == 0 && partition_config.fennel_dynamics != FENNELADP_ORIGINAL) {
//		fennel_weight = 0;
//	}
    fennel_weight = 1;

    for (int j = 0; j < partition_config.initial_part_fennel_tries; j++) {
        bool preliminary_sol = j || partition_config.restream_number;
        forall_nodes(G, node)
                {
                    if (node >= G.number_of_nodes() - partition_config.quotient_nodes) {
                        break;
                    }
                    node_too_large = true;

                    //now move the node to the cluster that is most common in the neighborhood
                    forall_out_edges(G, e, node)
                            {
                                NodeID target = G.getEdgeTarget(e);
                                if (target >= G.number_of_nodes() - partition_config.quotient_nodes || target < node ||
                                    preliminary_sol) {
                                    hash_map[G.getPartitionIndex(target)] += G.getEdgeWeight(e);
                                }
                            }
                    endfor

                    //second sweep for finding max and resetting array
                    PartitionID my_block = 0;
                    if (preliminary_sol) {
                        my_block = G.getPartitionIndex(node);
                        cluster_sizes[my_block] -= G.getNodeWeight(node);
                        cluster_ghost_nodes[my_block] -= G.getImplicitGhostNodes(node);
                    }
                    PartitionID max_block = my_block;
                    double max_value = std::numeric_limits<double>::lowest();

                    if (partition_config.use_queue) {
                        forall_out_edges(G, e, node) {
                                    NodeID target = G.getEdgeTarget(e);
                                    PartitionID cur_block = G.getPartitionIndex(target);
                                    cur_value = hash_map[cur_block];
                                    cur_value -=
                                            fennel_weight *
                                            (G.getNodeWeight(node) * partition_config.fennel_alpha_gamma *
                                             random_functions::approx_sqrt(cluster_sizes[cur_block]));

                                    if ((cur_value > max_value ||
                                         (cur_value == max_value && random_obj.nextBool())) &&
                                        (cluster_sizes[cur_block] - cluster_ghost_nodes[cur_block] +
                                         G.getNodeWeight(node) - G.getImplicitGhostNodes(node) <=
                                         partition_config.stream_total_upperbound)) {
                                        node_too_large = false;
                                        max_value = cur_value;
                                        max_block = cur_block;
                                    }
                                    hash_map[cur_block] = 0;
                                }
                        endfor PartitionID cur_block = queue->maxElement();
                        cur_value -=
                                fennel_weight *
                                (G.getNodeWeight(node) * partition_config.fennel_alpha_gamma *
                                 random_functions::approx_sqrt(cluster_sizes[cur_block]));

                        if ((cur_value > max_value ||
                             (cur_value == max_value && random_obj.nextBool())) &&
                            (cluster_sizes[cur_block] - cluster_ghost_nodes[cur_block] +
                             G.getNodeWeight(node) - G.getImplicitGhostNodes(node) <=
                             partition_config.stream_total_upperbound)) {
                            node_too_large = false;
                            max_value = cur_value;
                            max_block = cur_block;
                        }
                    } else {
                        for (PartitionID cur_block = 0; cur_block < hash_map.size(); cur_block++) {
                            cur_value = hash_map[cur_block];
                            cur_value -= fennel_weight * (G.getNodeWeight(node) * partition_config.fennel_alpha_gamma *
                                                          random_functions::approx_sqrt(cluster_sizes[cur_block]));

                            if ((cur_value > max_value || (cur_value == max_value && random_obj.nextBool()))
                                && (cluster_sizes[cur_block] - cluster_ghost_nodes[cur_block] + G.getNodeWeight(node) -
                                    G.getImplicitGhostNodes(node) <= partition_config.stream_total_upperbound)) {
                                node_too_large = false;
                                max_value = cur_value;
                                max_block = cur_block;
                            }

                            hash_map[cur_block] = 0;
                        }
                    }

                    if (node_too_large) {
                        max_block =
                                fnv2a(partition_config.lower_global_node + node) % partition_config.k; // Random choice
                    }

                    cluster_sizes[max_block] += G.getNodeWeight(node);
                    if (partition_config.use_queue) {
                        queue->subtractKey(max_block, -G.getNodeWeight(node));
                    }
                    cluster_ghost_nodes[max_block] += G.getImplicitGhostNodes(node);
                    G.setPartitionIndex(node, max_block);

                }
        endfor

    }
    delete queue;

    return 0;
}

