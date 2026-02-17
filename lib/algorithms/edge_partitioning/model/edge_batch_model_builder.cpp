/****************************************************************************
 * edge_batch_model_builder.cpp
 *****************************************************************************/

#include "algorithms/edge_partitioning/model/edge_batch_model_builder.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <random>

#include "io/edge_stream_reader.h"
#include "partition/model/model_graph_ops.h"
#include "partition/model/model_partition_mapper.h"

namespace {
void process_quotient_edge_in_batch_edge(partition::state::EdgePartitionerPassState& pass_state,
                                         NodeID node, LongNodeID global_target,
                                         EdgeWeight edge_weight) {
    (void)edge_weight;
    PartitionID target_global_par = global_target;
    if ((*pass_state.edge_block_nodes)[target_global_par].size() >= 1) {
        auto& curr_element = (*pass_state.edge_block_nodes)[target_global_par].back();
        if (curr_element.first == node) {
            curr_element.second += edge_weight;
            return;
        }
    }
    (*pass_state.edge_block_nodes)[target_global_par].emplace_back(node, edge_weight);
}

void add_edges_to_incident_block(Config& partition_config,
                                 partition::state::EdgePartitionerPassState& pass_state,
                                 NodeID node, LongEdgeID& used_edges) {
    constexpr NodeID kUnsetBlock = std::numeric_limits<NodeID>::max();
    if (partition_config.minimal_mode) {
        const NodeID incident_source = (*pass_state.nodes_on_edge_conv)[node][0];
        const NodeID target = (*pass_state.blocks_on_node_minimal)[incident_source];
        if (target != kUnsetBlock) {
            used_edges++;
            partition_config.quotient_edges_count++;
            process_quotient_edge_in_batch_edge(pass_state, node, target, 1);
        }
    } else {
        const NodeID incident_source = (*pass_state.nodes_on_edge_conv)[node][0];
        if (!(*pass_state.blocks_on_node)[incident_source].empty()) {
            if (partition_config.past_subset_size == -1 ||
                (partition_config.past_subset_size != 0 &&
                 (*pass_state.blocks_on_node)[incident_source].size() <= static_cast<size_t>(
                     partition_config.past_subset_size))) {
                for (google::dense_hash_set<PartitionID>::const_iterator it =
                         (*pass_state.blocks_on_node)[incident_source].begin();
                     it != (*pass_state.blocks_on_node)[incident_source].end();
                     ++it) {
                    used_edges++;
                    const PartitionID target = *it;
                    partition_config.quotient_edges_count++;
                    process_quotient_edge_in_batch_edge(pass_state, node, target, 1);
                }
            } else if (partition_config.past_subset_size != 0 &&
                       partition_config.past_subset_size != -1) {
                std::vector<PartitionID> in;
                for (google::dense_hash_set<PartitionID>::const_iterator it =
                         (*pass_state.blocks_on_node)[incident_source].begin();
                     it != (*pass_state.blocks_on_node)[incident_source].end();
                     ++it) {
                    in.push_back(*it);
                }
                std::vector<unsigned int> reservoir;
                reservoir.reserve(partition_config.past_subset_size);
                std::random_device rd;
                std::mt19937_64 gen(rd());
                for (size_t i = 0; i < in.size(); ++i) {
                    if (reservoir.size() < static_cast<size_t>(partition_config.past_subset_size)) {
                        reservoir.push_back(in[i]);
                    } else {
                        std::uniform_int_distribution<size_t> dis(0, i);
                        size_t j = dis(gen);
                        if (j < static_cast<size_t>(partition_config.past_subset_size)) {
                            reservoir[j] = in[i];
                        }
                    }
                }
                for (auto i : reservoir) {
                    used_edges++;
                    const PartitionID target = i;
                    partition_config.quotient_edges_count++;
                    process_quotient_edge_in_batch_edge(pass_state, node, target, 1);
                }
            }
        }
    }
}

void build_graph_model(graph_access& graph, Config& partition_config,
                       partition::state::EdgePartitionerPassState& pass_state,
                       std::vector<EdgeID>& contracted_edge_graph_id, std::vector<EdgeID>& rev_edge,
                       EdgeID split_nodes, NodeID split_edges, graph_access& split_graph) {
    // Construct the CSPAC graph model from the contracted split graph.
    partition_config.nmbNodes = split_nodes / 2;
    pass_state.prev_batch_edge_id = split_nodes / 2;
    const LongEdgeID remaining_batch_edges = split_edges / 2;

    if (partition_config.batch_alpha) {
        partition_config.fennel_alpha =
            remaining_batch_edges *
            std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
            (std::pow(pass_state.remaining_stream_partition_nodes, partition_config.fennel_gamma));

        partition_config.fennel_alpha_gamma =
            partition_config.fennel_alpha * partition_config.fennel_gamma;
    }

    NodeID node_counter = 0;
    EdgeID edge_counter = 0;
    NodeID node = 0;
    LongNodeID target;
    LongEdgeID used_edges = 0;
    LongEdgeID nmbEdges = 2 * remaining_batch_edges;
    const NodeWeight weight = 1;
    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>> all_edges;
    std::vector<NodeWeight> all_nodes;
    all_edges.resize(partition_config.nmbNodes + partition_config.quotient_nodes);
    all_nodes.resize(partition_config.nmbNodes + partition_config.quotient_nodes);
    pass_state.lower_global_node =
        pass_state.total_stream_nodecounter + 1;  // Bounds below start from 1 instead of 0
    pass_state.upper_global_node = pass_state.total_stream_nodecounter + partition_config.nmbNodes;
    pass_state.curr_batch++;

    if (nmbEdges > std::numeric_limits<EdgeWeight>::max() ||
        partition_config.nmbNodes > std::numeric_limits<LongNodeID>::max()) {
#ifdef MODE64BITEDGES
        std::cerr << "The graph is too large. Currently only 64bits supported!" << '\n';
#else
        std::cerr << "The graph is too large. Currently only 32bits supported!" << '\n';
#endif
        exit(0);
    }

    pass_state.edge_block_nodes =
        std::make_unique<std::vector<std::vector<std::pair<NodeID, NodeWeight>>>>(
            partition_config.k, std::vector<std::pair<NodeID, NodeWeight>>());

    NodeID in_deg = 0;
    NodeID out_deg = 0;
    NodeID current_edge_node;

    for (NodeID u = 0; u < split_graph.number_of_nodes(); ++u) {
        NodeID first_edge_u = split_graph.get_first_edge(u);
        for (EdgeID e = split_graph.get_first_edge(u); e < split_graph.get_first_invalid_edge(u);
             ++e) {
            NodeID v = split_graph.getEdgeTarget(e);
            if (u < v) {
                // Process each undirected edge once and connect it to its
                // predecessor/successor around both endpoints in the split graph.
                NodeID first_edge_v = split_graph.get_first_edge(v);
                in_deg = split_graph.getNodeDegree(u);
                out_deg = split_graph.getNodeDegree(v);
                current_edge_node = contracted_edge_graph_id[e];
                const EdgeWeight edge_weight = 1;
                node_counter++;
                node = current_edge_node;

                if (pass_state.restream_number > 0 && pass_state.stream_nodes_assign &&
                    pass_state.stream_blocks_weight) {
                    // Edge restreaming (minimal mode): this edge-node was already assigned in a
                    // previous pass. Remove its old block load before re-partitioning so quotient
                    // weights remain consistent and are not double-counted across passes.
                    const LongNodeID global_edge_node =
                        pass_state.lower_global_node + static_cast<LongNodeID>(node) - 1;
                    const PartitionID prev_block =
                        (*pass_state.stream_nodes_assign)[global_edge_node];
                    if (prev_block < partition_config.k) {
                        (*pass_state.stream_blocks_weight)[prev_block] -= 1;
                    }
                }

                all_nodes[node] = weight;

                if (in_deg > 1) {
                    EdgeID u_next = (in_deg + e - first_edge_u + 1) % in_deg + first_edge_u;
                    EdgeID u_prev = (in_deg + e - first_edge_u - 1) % in_deg + first_edge_u;

                    target = contracted_edge_graph_id[u_next] + 1 + pass_state.last_edge_count;
                    used_edges += ((NodeID)(target - pass_state.lower_global_node) < node);
                    edge_counter += partition_model::graph_ops::insert_regular_edge(
                        partition_config, all_edges, node, target, pass_state.lower_global_node,
                        edge_weight);
                    if (u_next != u_prev) {
                        target = contracted_edge_graph_id[u_prev] + 1 + pass_state.last_edge_count;
                        used_edges += ((NodeID)(target - pass_state.lower_global_node) < node);
                        edge_counter += partition_model::graph_ops::insert_regular_edge(
                            partition_config, all_edges, node, target, pass_state.lower_global_node,
                            edge_weight);
                    }
                }

                if (out_deg > 1) {
                    EdgeID u_next =
                        (out_deg + rev_edge[e] - first_edge_v + 1) % out_deg + first_edge_v;
                    EdgeID u_prev =
                        (out_deg + rev_edge[e] - first_edge_v - 1) % out_deg + first_edge_v;

                    target = contracted_edge_graph_id[u_next] + 1 + pass_state.last_edge_count;
                    used_edges += ((NodeID)(target - pass_state.lower_global_node) < node);
                    edge_counter += partition_model::graph_ops::insert_regular_edge(
                        partition_config, all_edges, node, target, pass_state.lower_global_node,
                        edge_weight);
                    if (u_next != u_prev) {
                        target = contracted_edge_graph_id[u_prev] + 1 + pass_state.last_edge_count;
                        used_edges += ((NodeID)(target - pass_state.lower_global_node) < node);
                        edge_counter += partition_model::graph_ops::insert_regular_edge(
                            partition_config, all_edges, node, target, pass_state.lower_global_node,
                            edge_weight);
                    }
                }

                add_edges_to_incident_block(partition_config, pass_state, node, used_edges);
                partition_config.quotient_edges_count = 0;
            }
        }
    }

    assert(pass_state.stream_blocks_weight);
    partition_model::graph_ops::insert_quotient_nodes(partition_config,
                                                      *pass_state.stream_blocks_weight, all_nodes,
                                                      partition_config.nmbNodes, node_counter);
    assert(pass_state.edge_block_nodes);
    edge_counter += partition_model::graph_ops::insert_quotient_edges(
        partition_config, *pass_state.edge_block_nodes, all_edges, partition_config.nmbNodes);

    std::vector<NodeWeight> implicit_ghost_nodes(node_counter, 0);
    partition_model::graph_ops::materialize_graph(
        graph, node_counter, edge_counter, all_edges, all_nodes, implicit_ghost_nodes,
        [&](NodeID node) {
            return partition_model::recover_partition_for_model_node(
                partition_config, pass_state.stream_nodes_assign.get(),
                pass_state.lower_global_node, node, node_counter, pass_state.restream_number);
        });

    pass_state.edge_block_nodes.reset();

    pass_state.total_stream_nodecounter += partition_config.nmbNodes;
    pass_state.remaining_stream_partition_nodes -= partition_config.nmbNodes;
    (void)used_edges;

    if (node_counter != (NodeID)partition_config.nmbNodes + partition_config.quotient_nodes) {
        std::cerr << "number of specified nodes mismatch" << '\n';
        std::cerr << (partition_config.nmbNodes + partition_config.quotient_nodes) << " "
                  << node_counter << '\n';
        exit(0);
    }
}

}  // namespace


namespace edge_partitioning::model {

void prepare_edge_batch_window(Config& partition_config,
                               partition::state::EdgePartitionerPassState& pass_state) {
    if (pass_state.total_stream_nodecounter == 0) {
        pass_state.lower_global_node_conv = 0;
    } else {
        pass_state.lower_global_node_conv = pass_state.upper_global_node_conv;
    }
    pass_state.upper_global_node_conv =
        pass_state.lower_global_node_conv + partition_config.nmbNodes;
    pass_state.last_edge_count = pass_state.incremental_edge_id;
    pass_state.prev_batch_edge_id = pass_state.last_edge_count;
    pass_state.remaining_stream_graph_nodes -= partition_config.nmbNodes;
}

void build_edge_split_batch_graph(Config& partition_config,
                                  partition::state::EdgePartitionerPassState& pass_state,
                                  std::vector<std::vector<Edge>>& subgraph_edges,
                                  std::vector<NodeID>& mapping, std::vector<EdgeID>& rev_edge,
                                  std::vector<EdgeID>& contracted_edge_graph_id,
                                  NodeID& number_of_deg1_vertices, NodeID& number_of_deg2_vertices,
                                  NodeID g_temp_no_nodes, NodeID g_temp_no_edges,
                                  graph_access& g_temp) {
    number_of_deg1_vertices = 0;
    number_of_deg2_vertices = 0;
    NodeID edge_to_node_key = 0;
    NodeID u_og;
    NodeID v_og;
    contracted_edge_graph_id.resize(g_temp_no_edges, -1);

    pass_state.nodes_on_edge_conv = std::make_unique<std::vector<std::vector<NodeID>>>(
        g_temp_no_edges / 2, std::vector<NodeID>(2));

    rev_edge.resize(g_temp_no_edges);

    g_temp.start_construction_light(g_temp_no_nodes, g_temp_no_edges);
    for (NodeID u = 0; u < g_temp_no_nodes; ++u) {
        if (subgraph_edges[u].size() == 1) {
            ++number_of_deg1_vertices;
        } else if (subgraph_edges[u].size() == 2) {
            ++number_of_deg2_vertices;
        }
        u_og = mapping[u];

        g_temp.new_node();
        g_temp.setNodeWeight(u, 1);
        for (auto& j : subgraph_edges[u]) {
            NodeID target_idx = j.target;
            v_og = mapping[target_idx];
            const EdgeID e = g_temp.new_edge(u, target_idx);
            if (u < target_idx) {
                j.weight = e;
                edge_to_node_key = pass_state.incremental_edge_id++;
                contracted_edge_graph_id[e] = edge_to_node_key - pass_state.last_edge_count;
            } else {
                EdgeWeight rev_edge_id = subgraph_edges[target_idx][j.weight].weight;
                rev_edge[rev_edge_id] = e;
                rev_edge[e] = rev_edge_id;
                contracted_edge_graph_id[e] = contracted_edge_graph_id[rev_edge[e]];
            }

            if (u_og < v_og) {
                (*pass_state.nodes_on_edge_conv)[contracted_edge_graph_id[e]][0] = u_og;
                (*pass_state.nodes_on_edge_conv)[contracted_edge_graph_id[e]][1] = v_og;
            }
            g_temp.setEdgeWeight(e, 1);
        }
    }
    g_temp.finish_construction_light();
}

void build_batch_model(Config& partition_config,
                       partition::state::EdgePartitionerPassState& pass_state, graph_access& G,
                       io::stream::StreamReader& stream_reader) {
    // Build the batch graph and then construct the model.
    std::unique_ptr<graph_access> G_temp(new graph_access());
    std::vector<EdgeID> rev_edge;
    std::vector<EdgeID> contracted_edge_graph_id;
    NodeID number_of_deg1_vertices, number_of_deg2_vertices;

    if (pass_state.stream_input_binary) {
        // Binary and text readers share downstream model construction.
        io::edge_stream::read_input_as_graph_binary(
            partition_config, pass_state, stream_reader, rev_edge, contracted_edge_graph_id,
            number_of_deg1_vertices, number_of_deg2_vertices, *G_temp);
    } else {
        io::edge_stream::read_input_as_graph(partition_config, pass_state, stream_reader, rev_edge,
                                             contracted_edge_graph_id, number_of_deg1_vertices,
                                             number_of_deg2_vertices, *G_temp);
    }

    NodeID split_n = G_temp->number_of_edges();
    EdgeID split_m =
        3 * G_temp->number_of_edges() - 2 * (number_of_deg1_vertices + number_of_deg2_vertices);

    if (partition_config.dynamic_alpha) {
        partition_config.fennel_edges -=
            (2 * (number_of_deg2_vertices + number_of_deg1_vertices) + split_n);
        partition_config.fennel_alpha =
            partition_config.fennel_edges *
            std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
            (std::pow(pass_state.remaining_stream_partition_nodes, partition_config.fennel_gamma));
        partition_config.fennel_alpha_gamma =
            partition_config.fennel_alpha * partition_config.fennel_gamma;
    }
    build_graph_model(G, partition_config, pass_state, contracted_edge_graph_id, rev_edge, split_n,
                      split_m - split_n, *G_temp);
}

} // namespace edge_partitioning::model
