/****************************************************************************
 * node_batch_model_builder.cpp
 *****************************************************************************/

#include "algorithms/node_partitioning/model/node_batch_model_builder.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <memory>

#include "partition/assignment/node_fennel_assign.h"
#include "partition/model/model_graph_ops.h"
#include "partition/model/model_partition_mapper.h"

namespace {

template <typename T>
T return_and_delete_element(std::vector<T>& vec, LongNodeID pos) {
    T val = vec[pos];
    vec[pos] = vec[vec.size() - 1];
    vec.erase(vec.begin() + vec.size() - 1);
    return val;
}

EdgeID insert_regular_local_edge_in_batch(
    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges, NodeID node,
    NodeID local_target, EdgeWeight edge_weight, EdgeWeight non_ghost_multiplier) {
    if (local_target == node) {
        std::cerr << "The graph file contains self-loops, which are not supported. "
                     "Please remove them from the file."
                  << '\n';
        exit(0);
    }
    return partition_model::graph_ops::include_edge(all_edges, node, local_target,
                                                    non_ghost_multiplier * edge_weight);
}

void process_node_weight(Config& config, partition::state::NodePartitionerPassState& pass_state,
                         std::vector<NodeWeight>& all_nodes, NodeID node, NodeWeight weight) {
    all_nodes[node] = weight;
    if (config.edge_partition) {
        return;
    }
    LongNodeID global_node = pass_state.lower_global_node + static_cast<LongNodeID>(node) - 1;
    PartitionID node_global_par = (*pass_state.stream_nodes_assign)[global_node];
    if (pass_state.restream_number) {
        node_fennel_assign::sub_block_weight(config, pass_state, node_global_par, weight);
    }
}

void process_node_weight(Config& config, partition::state::NodePartitionerPassState& pass_state,
                         std::vector<NodeWeight>& all_nodes, NodeID node, NodeWeight weight,
                         LongNodeID global_node) {
    all_nodes[node] = weight;
    if (config.edge_partition) {
        return;
    }
    PartitionID node_global_par = (*pass_state.stream_nodes_assign)[global_node - 1];
    if (pass_state.restream_number) {
        node_fennel_assign::sub_block_weight(config, pass_state, node_global_par, weight);
    }
}

void setup_for_ghost_neighbors(Config& config,
                               partition::state::NodePartitionerPassState& pass_state) {
    if (config.stream_allow_ghostnodes || pass_state.restream_number) {
        if (!pass_state.ghostglobal_to_ghostkey) {
            pass_state.ghostglobal_to_ghostkey = std::make_unique<buffered_map>(
                pass_state.stream_nodes_assign.get(), pass_state.restream_number);
        }
        if (pass_state.ghostkey_to_edges) {
            pass_state.ghostkey_to_edges->clear();
        } else {
            pass_state.ghostkey_to_edges =
                std::make_unique<std::vector<std::vector<std::pair<NodeID, ShortEdgeWeight>>>>();
        }
        config.ghost_nodes = (pass_state.restream_number > 0) ? config.k : 0;
        if (pass_state.restream_number > 0) {
            pass_state.ghostkey_to_edges->resize(config.ghost_nodes);
        }
    }
}

void process_ghost_neighbor_in_batch(Config& config,
                                     partition::state::NodePartitionerPassState& pass_state,
                                     NodeID node, LongNodeID ghost_target, EdgeWeight edge_weight) {
    LongNodeID ghost_key;
    PartitionID target_global_par = (*pass_state.stream_nodes_assign)[ghost_target - 1];
    if (config.stream_allow_ghostnodes ||
        (pass_state.restream_number && target_global_par != INVALID_PARTITION)) {
        if (!pass_state.ghostglobal_to_ghostkey->has_key(ghost_target - 1)) {
            ghost_key = config.ghost_nodes++;
            pass_state.ghostglobal_to_ghostkey->push_back(ghost_target - 1, ghost_key);
            pass_state.ghostkey_to_edges->push_back(
                std::vector<std::pair<NodeID, ShortEdgeWeight>>());
        } else {
            ghost_key = (*pass_state.ghostglobal_to_ghostkey)[ghost_target - 1];
        }
        (*pass_state.ghostkey_to_edges)[ghost_key].push_back(
            std::make_pair(node, static_cast<ShortEdgeWeight>(edge_weight)));
    }
}

bool process_quotient_edge_in_batch(Config& config,
                                    partition::state::NodePartitionerPassState& pass_state,
                                    NodeID node, LongNodeID global_target, EdgeWeight edge_weight) {
    PartitionID target_global_par = (*pass_state.stream_nodes_assign)[global_target - 1];
    if (config.k <= target_global_par && target_global_par < 2 * config.k) {
        config.num_ghost_nodes++;
        target_global_par -= config.k;
        edge_weight = 1;
    }
    if (target_global_par > config.k) {
        return false;
    }
    if ((*pass_state.edge_block_nodes)[target_global_par].size() >= 1) {
        auto& curr_element = (*pass_state.edge_block_nodes)[target_global_par].back();
        if (curr_element.first == node) {
            curr_element.second += edge_weight;
            return true;
        }
    }
    std::pair<NodeID, NodeWeight> element = {node, edge_weight};
    (*pass_state.edge_block_nodes)[target_global_par].push_back(element);
    return true;
}

NodeID restream_map_ghost_keys_to_nodes(Config& config,
                                        partition::state::NodePartitionerPassState& pass_state) {
    pass_state.ghostkey_to_node = std::make_unique<std::vector<NodeID>>(config.ghost_nodes, 0);
    for (PartitionID target_par = 0; target_par < config.ghost_nodes; target_par++) {
        NodeID node = target_par + config.nmbNodes;
        (*pass_state.ghostkey_to_node)[target_par] = node;
    }
    return 0;
}

NodeID greedy_map_ghost_keys_to_nodes(
    Config& config, partition::state::NodePartitionerPassState& pass_state,
    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges,
    std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes,
    NodeID& node_counter) {
    pass_state.ghostkey_to_node = std::make_unique<std::vector<NodeID>>(config.ghost_nodes, 0);
    NodeID inserted_nodes = std::min(config.ghost_nodes, config.ghost_nodes_threshold);
    all_nodes.resize(config.nmbNodes + inserted_nodes + config.quotient_nodes);
    all_edges.resize(config.nmbNodes + inserted_nodes + config.quotient_nodes);
    for (LongNodeID ghost_key = 0; ghost_key < config.ghost_nodes; ghost_key++) {
        if (ghost_key < static_cast<LongNodeID>(inserted_nodes)) {
            NodeID node = node_counter++;
            (*pass_state.ghostkey_to_node)[ghost_key] = node;
            all_nodes[node] = 1;
        } else {
            auto& list_edges = (*pass_state.ghostkey_to_edges)[ghost_key];
            LongNodeID contr_pos =
                (static_cast<LongNodeID>(ghost_key) * pass_state.curr_batch) % list_edges.size();
            NodeID node = return_and_delete_element(list_edges, contr_pos).first;
            (*pass_state.ghostkey_to_node)[ghost_key] = node;
            all_nodes[node]++;
            all_assigned_ghost_nodes[node]++;
        }
    }
    return inserted_nodes;
}

NodeID map_ghost_keys_to_nodes_in_batch(
    Config& config, partition::state::NodePartitionerPassState& pass_state,
    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges,
    std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes,
    NodeID& node_counter) {
    NodeID inserted_nodes = 0;
    if (pass_state.restream_number) {
        inserted_nodes = restream_map_ghost_keys_to_nodes(config, pass_state);
        pass_state.ghostglobal_to_ghostkey->clear();
    } else if (config.stream_allow_ghostnodes) {
        inserted_nodes = greedy_map_ghost_keys_to_nodes(config, pass_state, all_edges, all_nodes,
                                                        all_assigned_ghost_nodes, node_counter);
        pass_state.ghostglobal_to_ghostkey->clear();
    }
    return inserted_nodes;
}

EdgeID insert_ghost_edges_in_batch(
    Config& config, partition::state::NodePartitionerPassState& pass_state,
    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>>& all_edges) {
    EdgeID inserted_edges = 0;
    if (config.stream_allow_ghostnodes || pass_state.restream_number) {
        for (LongNodeID ghost_key = 0; ghost_key < pass_state.ghostkey_to_edges->size();
             ghost_key++) {
            NodeID ghost_node = (*pass_state.ghostkey_to_node)[ghost_key];
            auto& neighbors_list = (*pass_state.ghostkey_to_edges)[ghost_key];
            for (auto& edge : neighbors_list) {
                inserted_edges += partition_model::graph_ops::include_edge(all_edges, ghost_node,
                                                                           edge.first, edge.second);
            }
        }
    }
    return inserted_edges;
}

}  // namespace


namespace node_partitioning::model {

NodeID build_batch_model_from_stream_input(Config& config,
                                           partition::state::NodePartitionerPassState& pass_state,
                                           graph_access& G,
                                           std::vector<std::vector<LongNodeID>>& input) {
    // Build a per-batch model graph with quotient and ghost nodes.
    NodeWeight total_nodeweight = 0;
    NodeID node_counter = 0;
    EdgeID edge_counter = 0;
    LongEdgeID used_edges = 0;
    bool read_ew = false;
    bool read_nw = false;
    LongEdgeID nmbEdges = 2 * pass_state.remaining_stream_edges;
    LongNodeID target;
    NodeWeight weight;
    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>> all_edges;
    std::vector<NodeWeight> all_nodes;

    std::vector<NodeWeight> all_assigned_ghost_nodes(config.nmbNodes + config.quotient_nodes, 0);
    all_edges.resize(config.nmbNodes + config.quotient_nodes);
    all_nodes.resize(config.nmbNodes + config.quotient_nodes);
    pass_state.lower_global_node =
        pass_state.total_stream_nodecounter + 1;  // Bounds below start from 1
    pass_state.upper_global_node = pass_state.total_stream_nodecounter + config.nmbNodes;
    LongNodeID cursor = 0;
    NodeID node = 0;

    pass_state.curr_batch++;

    if (config.ram_stream) {
        cursor = input.size() - pass_state.remaining_stream_nodes;
    }

    if (nmbEdges > std::numeric_limits<EdgeWeight>::max() ||
        config.nmbNodes > std::numeric_limits<LongNodeID>::max()) {
#ifdef MODE64BITEDGES
        std::cerr << "The graph is too large. Currently only 64bits supported!" << '\n';
#else
        std::cerr << "The graph is too large. Currently only 32bits supported!" << '\n';
#endif
        exit(0);
    }
    switch (pass_state.remaining_stream_ew) {
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

    // Track incident block edges for quotient-node construction.
    pass_state.edge_block_nodes =
        std::make_unique<std::vector<std::vector<std::pair<NodeID, NodeWeight>>>>(
            config.k, std::vector<std::pair<NodeID, NodeWeight>>());

    setup_for_ghost_neighbors(config, pass_state);

    for (node_counter = 0; node_counter < config.nmbNodes; node_counter++) {
        std::vector<LongNodeID>& line_numbers = input[cursor];
        LongNodeID col_counter = 0;
        node = (NodeID)node_counter;
        weight = 1;
        if (read_nw) {
            weight = line_numbers[col_counter++];
            if (total_nodeweight > std::numeric_limits<NodeWeight>::max()) {
                std::cerr
                    << "The sum of the node weights is too large (it exceeds the node weight type)."
                    << '\n';
                std::cerr << "Currently not supported. Please scale your node weights." << '\n';
                exit(0);
            }
        }
        total_nodeweight += weight;
        process_node_weight(config, pass_state, all_nodes, node, weight);

        // Classify edges into current/past/future batches.
        while (col_counter < line_numbers.size()) {
            target = line_numbers[col_counter++];
            EdgeWeight edge_weight = 1;
            if (read_ew) {
                edge_weight = line_numbers[col_counter++];
            }

            if (target > pass_state.upper_global_node) {
                process_ghost_neighbor_in_batch(config, pass_state, node, target, edge_weight);
            } else if (target < pass_state.lower_global_node) {
                used_edges++;
                process_quotient_edge_in_batch(config, pass_state, node, target, edge_weight);
            } else {
                used_edges += ((NodeID)(target - pass_state.lower_global_node) < node);
                edge_counter += partition_model::graph_ops::insert_regular_edge(
                    config, all_edges, node, target, pass_state.lower_global_node, edge_weight);
            }
        }

        cursor++;
    }
    // Map ghost keys to nodes and insert quotient/ghost edges.
    NodeID uncontracted_ghost_nodes = map_ghost_keys_to_nodes_in_batch(
        config, pass_state, all_edges, all_nodes, all_assigned_ghost_nodes, node_counter);
    assert(pass_state.stream_blocks_weight);
    partition_model::graph_ops::insert_quotient_nodes(
        config, *pass_state.stream_blocks_weight, all_nodes,
        config.nmbNodes + uncontracted_ghost_nodes, node_counter);
    edge_counter += insert_ghost_edges_in_batch(config, pass_state, all_edges);
    assert(pass_state.edge_block_nodes);
    edge_counter += partition_model::graph_ops::insert_quotient_edges(
        config, *pass_state.edge_block_nodes, all_edges,
        config.nmbNodes + uncontracted_ghost_nodes);

    partition_model::graph_ops::materialize_graph(
        G, node_counter, edge_counter, all_edges, all_nodes, all_assigned_ghost_nodes,
        [&](NodeID node) {
            return partition_model::recover_partition_for_model_node(
                config, pass_state.stream_nodes_assign.get(), pass_state.lower_global_node, node,
                node_counter, pass_state.restream_number);
        });

    pass_state.edge_block_nodes.reset();

    pass_state.total_stream_nodecounter += config.nmbNodes;
    pass_state.total_stream_nodeweight += total_nodeweight;
    pass_state.remaining_stream_nodes -= config.nmbNodes;
    pass_state.remaining_stream_edges -= used_edges;

    if (node_counter !=
        (NodeID)config.nmbNodes + uncontracted_ghost_nodes + config.quotient_nodes) {
        std::cerr << "number of specified nodes mismatch" << '\n';
        std::cerr << (config.nmbNodes + uncontracted_ghost_nodes + config.quotient_nodes) << " "
                  << node_counter << '\n';
        exit(0);
    }

    return node_counter;
}

NodeID build_batch_model_from_reordered_nodes(
    Config& config, partition::state::NodePartitionerPassState& pass_state, graph_access& G,
    std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>& batch_nodes, size_t batch_id) {
    NodeWeight total_nodeweight = 0;
    NodeID node_counter = 0;
    EdgeID edge_counter = 0;
    LongEdgeID used_edges = 0;
    LongNodeID target;
    NodeWeight weight;
    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>> all_edges;
    std::vector<NodeWeight> all_nodes;

    std::vector<NodeWeight> all_assigned_ghost_nodes(config.nmbNodes + config.quotient_nodes, 0);
    all_edges.resize(config.nmbNodes + config.quotient_nodes);
    all_nodes.resize(config.nmbNodes + config.quotient_nodes);

    PartitionID batch_marker = config.batch_manager->get_batch_marker(batch_id);

    pass_state.curr_batch++;

    pass_state.edge_block_nodes =
        std::make_unique<std::vector<std::vector<std::pair<NodeID, NodeWeight>>>>(
            config.k, std::vector<std::pair<NodeID, NodeWeight>>());

    setup_for_ghost_neighbors(config, pass_state);

    std::vector<NodeID> global_to_local_map(config.number_of_nodes, UNDEFINED_NODE);

    node_counter = 0;
    config.num_ghost_nodes = 0;

    std::vector<PartitionID>& node_to_batch_marker = config.sep_batch_marker
                                                         ? (*pass_state.stream_nodes_batch_marker)
                                                         : (*pass_state.stream_nodes_assign);

    for (auto& entry : batch_nodes) {
        LongNodeID global_node_id = entry.first;
        std::vector<LongNodeID>& line_numbers = entry.second;
        LongNodeID col_counter = 0;
        NodeID node = node_counter;

        if (config.local_to_global_map) {
            (*config.local_to_global_map)[node] = global_node_id;
        }
        global_to_local_map[global_node_id - 1] = node;

        weight = 1;
        total_nodeweight += weight;
        process_node_weight(config, pass_state, all_nodes, node, weight, global_node_id);

        if (pass_state.restream_number) {
            while (col_counter < static_cast<LongNodeID>(line_numbers.size())) {
                target = line_numbers[col_counter++];
                EdgeWeight edge_weight = 1;
                NodeID local_target = global_to_local_map[target - 1];
                if (local_target != UNDEFINED_NODE) {
                    used_edges++;
                    edge_counter += insert_regular_local_edge_in_batch(
                        all_edges, node, local_target, edge_weight,
                        1 + config.double_non_ghost_edges);
                } else {
                    used_edges++;
                    process_quotient_edge_in_batch(config, pass_state, node, target, edge_weight);
                }
            }
        } else {
            while (col_counter < static_cast<LongNodeID>(line_numbers.size())) {
                target = line_numbers[col_counter++];
                EdgeWeight edge_weight = 1;
                if (config.ghost_neighbors_enabled) {
                    edge_weight = config.default_weight_non_ghost;
                }

                if (node_to_batch_marker[target - 1] == batch_marker) {
                    NodeID local_target = global_to_local_map[target - 1];
                    if (local_target != UNDEFINED_NODE) {
                        used_edges++;
                        edge_counter += insert_regular_local_edge_in_batch(
                            all_edges, node, local_target, edge_weight,
                            1 + config.double_non_ghost_edges);
                    }
                } else if ((*pass_state.stream_nodes_assign)[target - 1] < 2 * config.k) {
                    bool inserted = process_quotient_edge_in_batch(config, pass_state, node, target,
                                                                   edge_weight);
                    if (inserted) {
                        used_edges++;
                    }

                    if (config.ghost_neighbors_enabled &&
                        (*pass_state.stream_nodes_assign)[target - 1] >= config.k) {
                        (*config.batch_unpartitioned_neighbors)[node].push_back(target);
                    }
                } else {
                    if (config.ghost_neighbors_enabled) {
                        (*config.batch_unpartitioned_neighbors)[node].push_back(target);
                    }
                }
            }
        }

        node_counter++;
    }

    NodeID uncontracted_ghost_nodes = map_ghost_keys_to_nodes_in_batch(
        config, pass_state, all_edges, all_nodes, all_assigned_ghost_nodes, node_counter);
    assert(pass_state.stream_blocks_weight);
    partition_model::graph_ops::insert_quotient_nodes(
        config, *pass_state.stream_blocks_weight, all_nodes,
        config.nmbNodes + uncontracted_ghost_nodes, node_counter);
    edge_counter += insert_ghost_edges_in_batch(config, pass_state, all_edges);
    assert(pass_state.edge_block_nodes);
    edge_counter += partition_model::graph_ops::insert_quotient_edges(
        config, *pass_state.edge_block_nodes, all_edges,
        config.nmbNodes + uncontracted_ghost_nodes);

    pass_state.edge_block_nodes.reset();

    partition_model::graph_ops::materialize_graph(
        G, node_counter, edge_counter, all_edges, all_nodes, all_assigned_ghost_nodes,
        [&](NodeID node) {
            return partition_model::recover_partition_for_model_node(
                config, pass_state.stream_nodes_assign.get(), pass_state.lower_global_node, node,
                node_counter, pass_state.restream_number);
        });

    pass_state.total_stream_nodecounter += config.nmbNodes;
    pass_state.total_stream_nodeweight += total_nodeweight;

    if (node_counter !=
        (NodeID)(config.nmbNodes + uncontracted_ghost_nodes + config.quotient_nodes)) {
        std::cerr << "number of specified nodes mismatch" << '\n';
        std::cerr << (config.nmbNodes + uncontracted_ghost_nodes + config.quotient_nodes) << " "
                  << node_counter << '\n';
        exit(0);
    }

    return node_counter;
}

} // namespace node_partitioning::model
