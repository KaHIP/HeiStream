/******************************************************************************
 * init_fennel.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "init_fennel.h"

#include <algorithm>
#include <vector>

#include "core/timing/scoped_stage_timer.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "partition/heuristics/fennel_scoring.h"

init_fennel::init_fennel() {}

init_fennel::~init_fennel() {}

void init_fennel::initial_partition(Config& config, const unsigned int seed, graph_access& G,
                                    int* partition_map, int ismultisec /*=0*/,
                                    int stream_pass_index /*=0*/) {
    double t = 0.0;
    TIMED_SCOPE(t, "init_fennel_total");
    unsigned iterations = 1;
    EdgeWeight best_cut = std::numeric_limits<EdgeWeight>::max();

    for (unsigned i = 0; i < iterations; i++) {
        fennel(config, G, stream_pass_index > 0);
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
            forall_nodes(G, n) {
                partition_map[n] = G.getPartitionIndex(n);
            }
            endfor
        }
    }
    PRINT(std::cout << "init_fennel took " << t << '\n';)
}

void init_fennel::initial_partition(Config& config, const unsigned int seed, graph_access& G,
                                    int* xadj, int* adjncy, int* vwgt, int* adjwgt,
                                    int* partition_map, int ismultisec /*=0*/,
                                    int stream_pass_index /*=0*/) {
    (void)stream_pass_index;
    std::cout << "not implemented yet" << '\n';
}

EdgeWeight init_fennel::fennel(Config& partition_config, graph_access& G, bool is_restream_pass) {
    bool node_too_large = false;
    std::vector<PartitionID> hash_map(partition_config.k, 0);
    std::vector<NodeWeight> cluster_sizes(partition_config.k, 0);
    std::vector<NodeWeight> cluster_ghost_nodes(partition_config.k, 0);
    std::vector<int> touched_stamp(partition_config.k, 0);
    int current_stamp = 0;
    std::vector<PartitionID> touched_blocks;
    touched_blocks.reserve(32);
    std::unique_ptr<maxNodeHeap> queue;
    if (partition_config.use_queue) {
        queue = std::make_unique<maxNodeHeap>();
    }

    if (is_restream_pass) {
        for (NodeID node = 0; node < G.number_of_nodes() - partition_config.quotient_nodes;
             node++) {
            cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
        }
    }
    for (NodeID node = G.number_of_nodes() - partition_config.quotient_nodes,
                end = G.number_of_nodes();
         node < end; node++) {
        cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
        cluster_ghost_nodes[G.getPartitionIndex(node)] += G.getImplicitGhostNodes(node);
    }
    if (partition_config.use_queue) {
        for (PartitionID block = 0; block < partition_config.k; ++block) {
            queue->insert(block, -cluster_sizes[block]);
        }
    }

    double fennel_weight = 2;
    double fennel_tmp = 0;
    const double remaining_stream_nodes = std::max<LongNodeID>(
        0, partition_config.total_nodes - partition_config.stream_assigned_nodes);
    switch (partition_config.fennel_dynamics) {
        case FENNELADP_ORIGINAL:
            fennel_weight = 1;
            break;
        case FENNELADP_DOUBLE:
            fennel_weight = 2;
            break;
        case FENNELADP_LINEAR:
            fennel_weight =
                2 * remaining_stream_nodes / (double)partition_config.total_stream_edges;
            break;
        case FENNELADP_MID_LINEAR:
            fennel_tmp = 2 * remaining_stream_nodes / (double)partition_config.total_stream_edges;
            if (fennel_tmp <= 1) {
                fennel_weight = 2 * (fennel_tmp);
            }
            break;
        case FENNELADP_QUADRATIC:
            fennel_tmp = remaining_stream_nodes / (double)partition_config.total_stream_edges;
            fennel_weight = 2 * fennel_tmp * fennel_tmp;
            break;
        case FENNELADP_MID_QUADRATIC:
            fennel_tmp = 2 * remaining_stream_nodes / (double)partition_config.total_stream_edges;
            if (fennel_tmp <= 1) {
                fennel_weight = 2 * fennel_tmp * fennel_tmp;
            }
            break;
        case FENNELADP_MID_CONSTANT:
            fennel_tmp = remaining_stream_nodes / (double)partition_config.total_stream_edges;
            if (fennel_tmp <= 1.5) {
                fennel_weight = 0.5;
            }
            break;
        case FENNELADP_EDGE_CUT:
            fennel_weight = 0;
            break;
    }

    //	if (partition_config.remaining_stream_nodes == 0 && partition_config.fennel_dynamics !=
    // FENNELADP_ORIGINAL) { 		fennel_weight = 0;
    //	}
    fennel_weight = 1;

    for (int j = 0; j < partition_config.initial_part_fennel_tries; j++) {
        bool preliminary_sol = j || is_restream_pass;
        forall_nodes(G, node) {
            if (node >= G.number_of_nodes() - partition_config.quotient_nodes) {
                break;
            }
            node_too_large = true;

            // now move the node to the cluster that is most common in the neighborhood
            forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if (target >= G.number_of_nodes() - partition_config.quotient_nodes ||
                    target < node || preliminary_sol) {
                    hash_map[G.getPartitionIndex(target)] += G.getEdgeWeight(e);
                }
            }
            endfor

                // second sweep for finding max and resetting array
            PartitionID my_block = 0;
            if (preliminary_sol) {
                my_block = G.getPartitionIndex(node);
                cluster_sizes[my_block] -= G.getNodeWeight(node);
                cluster_ghost_nodes[my_block] -= G.getImplicitGhostNodes(node);
                if (partition_config.use_queue) {
                    queue->changeKey(my_block, -cluster_sizes[my_block]);
                }
            }
            PartitionID max_block = my_block;
            double max_value = std::numeric_limits<double>::lowest();
            NodeWeight max_block_load = std::numeric_limits<NodeWeight>::max();
            touched_blocks.clear();

            auto is_feasible = [&](PartitionID cur_block) {
                return (cluster_sizes[cur_block] - cluster_ghost_nodes[cur_block] +
                            G.getNodeWeight(node) - G.getImplicitGhostNodes(node) <=
                        partition_config.stream_total_upperbound);
            };

            auto maybe_take_block = [&](PartitionID cur_block, bool enforce_feasible) {
                if (enforce_feasible && !is_feasible(cur_block)) {
                    return;
                }
                const NodeWeight cur_load = cluster_sizes[cur_block];
                const double cur_value = partition::heuristics::fennel::assignment_score(
                    hash_map[cur_block], partition_config.fennel_alpha_gamma,
                    cluster_sizes[cur_block], G.getNodeWeight(node), fennel_weight);
                if (cur_value > max_value ||
                    (cur_value == max_value && cur_load < max_block_load) ||
                    (cur_value == max_value && cur_load == max_block_load && cur_block < max_block)) {
                    node_too_large = false;
                    max_value = cur_value;
                    max_block = cur_block;
                    max_block_load = cur_load;
                }
            };
            if (partition_config.use_queue) {
                ++current_stamp;
                if (current_stamp == std::numeric_limits<int>::max()) {
                    std::fill(touched_stamp.begin(), touched_stamp.end(), 0);
                    current_stamp = 1;
                }
            }

            forall_out_edges(G, e, node) {
                const NodeID target = G.getEdgeTarget(e);
                if (target >= G.number_of_nodes() - partition_config.quotient_nodes || target < node ||
                    preliminary_sol) {
                    const PartitionID cur_block = G.getPartitionIndex(target);
                    if (partition_config.use_queue && touched_stamp[cur_block] != current_stamp) {
                        touched_stamp[cur_block] = current_stamp;
                        touched_blocks.push_back(cur_block);
                    }
                }
            }
            endfor

            if (partition_config.use_queue) {
                for (const PartitionID cur_block : touched_blocks) {
                    maybe_take_block(cur_block, true);
                    hash_map[cur_block] = 0;
                }

                if (!queue->empty()) {
                    const PartitionID min_load_block = queue->maxElement();
                    if (touched_stamp[min_load_block] != current_stamp) {
                        maybe_take_block(min_load_block, false);
                    }
                }
            } else {
                for (PartitionID cur_block = 0; cur_block < hash_map.size(); cur_block++) {
                    maybe_take_block(cur_block, true);
                    hash_map[cur_block] = 0;
                }
            }

            if (node_too_large) {
                max_block = fnv2a(partition_config.lower_global_store + node) %
                            partition_config.k;  // Random choice
            }

            cluster_sizes[max_block] += G.getNodeWeight(node);
            if (partition_config.use_queue) {
                queue->changeKey(max_block, -cluster_sizes[max_block]);
            }
            cluster_ghost_nodes[max_block] += G.getImplicitGhostNodes(node);
            G.setPartitionIndex(node, max_block);
        }
        endfor
    }
    return 0;
}
