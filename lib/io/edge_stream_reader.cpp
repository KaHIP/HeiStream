/****************************************************************************
 * edge_stream_reader.cpp
 *****************************************************************************/

#include "io/edge_stream_reader.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include <sstream>
#include <vector>

#include "algorithms/edge_partitioning/model/edge_batch_model_builder.h"
#include "core/timing/scoped_stage_timer.h"
#include "io/stream_reader.h"


namespace io::edge_stream {

namespace {

void add_split_edge(std::vector<std::vector<Edge>>& g_temp_edges, std::vector<NodeID>& mapping,
                    google::dense_hash_map<NodeID, NodeID>& prev_batch_mapping, NodeID node,
                    NodeID& target_idx, NodeID& map_index, NodeID& g_temp_no_nodes,
                    NodeID& g_temp_no_edges, LongNodeID target_global,
                    LongNodeID lower_global_node_conv, LongNodeID upper_global_node_conv) {
    if (target_global >= upper_global_node_conv) {
        return;
    }

    if (target_global < lower_global_node_conv) {
        if (prev_batch_mapping.find(target_global) == prev_batch_mapping.end()) {
            prev_batch_mapping[target_global] = map_index;
            target_idx = map_index;
            g_temp_edges.emplace_back();
            mapping.push_back(static_cast<NodeID>(target_global));
            map_index++;
            g_temp_no_nodes++;
        } else {
            target_idx = prev_batch_mapping[target_global];
        }

        Edge forward_edge{};
        forward_edge.target = target_idx;
        forward_edge.weight = -1;
        g_temp_edges[node].emplace_back(forward_edge);

        Edge backward_edge{};
        backward_edge.target = node;
        backward_edge.weight = g_temp_edges[node].size() - 1;
        g_temp_edges[target_idx].emplace_back(backward_edge);
        g_temp_no_edges += 2;
        return;
    }

    target_idx = static_cast<NodeID>(target_global - lower_global_node_conv);
    if (node < target_idx) {
        Edge forward_edge{};
        forward_edge.target = target_idx;
        forward_edge.weight = -1;
        g_temp_edges[node].emplace_back(forward_edge);

        Edge backward_edge{};
        backward_edge.target = node;
        backward_edge.weight = g_temp_edges[node].size() - 1;
        g_temp_edges[target_idx].emplace_back(backward_edge);
        g_temp_no_edges += 2;
    }
}

}  // namespace

void read_input_as_graph(Config& partition_config,
                         partition::state::EdgePartitionerPassState& pass_state,
                         io::stream::StreamReader& stream_reader, std::vector<EdgeID>& rev_edge,
                         std::vector<EdgeID>& contracted_edge_graph_id,
                         NodeID& number_of_deg1_vertices, NodeID& number_of_deg2_vertices,
                         graph_access& G_temp) {
    edge_partitioning::model::prepare_edge_batch_window(partition_config, pass_state);
    std::vector<std::vector<Edge>> G_temp_edges;
    G_temp_edges.resize(partition_config.nmbNodes);
    std::vector<NodeID> mapping(partition_config.nmbNodes, -1);
    NodeID target_idx = -1;
    LongNodeID node_counter = 0;
    NodeID G_temp_no_nodes = partition_config.nmbNodes;
    NodeID G_temp_no_edges = 0;
    NodeID map_index = partition_config.nmbNodes;
    NodeID node = 0;
    google::dense_hash_map<NodeID, NodeID> prev_batch_mapping;
    prev_batch_mapping.set_empty_key(-1);

    while (node_counter < partition_config.nmbNodes) {
        std::vector<LongNodeID> input;
        if (!stream_reader.read_next_parsed_line(input)) {
            break;
        }

        LongNodeID col_counter = 0;
        node = static_cast<NodeID>(node_counter);
        mapping[node] = node + pass_state.lower_global_node_conv;

        while (col_counter < input.size()) {
            const LongNodeID target = input[col_counter++] - 1;
            add_split_edge(G_temp_edges, mapping, prev_batch_mapping, node, target_idx, map_index,
                           G_temp_no_nodes, G_temp_no_edges, target,
                           pass_state.lower_global_node_conv, pass_state.upper_global_node_conv);
        }
        node_counter++;
    }
    edge_partitioning::model::build_edge_split_batch_graph(
        partition_config, pass_state, G_temp_edges, mapping, rev_edge, contracted_edge_graph_id,
        number_of_deg1_vertices, number_of_deg2_vertices, G_temp_no_nodes, G_temp_no_edges, G_temp);
}

void read_input_as_graph_binary(Config& partition_config,
                                partition::state::EdgePartitionerPassState& pass_state,
                                io::stream::StreamReader& stream_reader,
                                std::vector<EdgeID>& rev_edge,
                                std::vector<EdgeID>& contracted_edge_graph_id,
                                NodeID& number_of_deg1_vertices, NodeID& number_of_deg2_vertices,
                                graph_access& G_temp) {
    edge_partitioning::model::prepare_edge_batch_window(partition_config, pass_state);
    std::vector<std::vector<Edge>> G_temp_edges;
    G_temp_edges.resize(partition_config.nmbNodes);
    std::vector<NodeID> mapping(partition_config.nmbNodes, -1);

    NodeID target_idx = -1;
    NodeID G_temp_no_nodes = partition_config.nmbNodes;
    NodeID G_temp_no_edges = 0;
    NodeID map_index = partition_config.nmbNodes;
    NodeID node = 0;
    google::dense_hash_map<NodeID, NodeID> prev_batch_mapping;
    prev_batch_mapping.set_empty_key(-1);

    unsigned long long nodes_in_batch = static_cast<unsigned long long>(partition_config.nmbNodes);
    std::vector<unsigned long long> vertex_offsets(nodes_in_batch + 1);
    stream_reader.stream().seekg(stream_reader.binary_start_pos());
    {
        double io_seconds = 0.0;
        {
            TIMED_SCOPE(io_seconds, "edge_binary_read_offsets");
            stream_reader.stream().read(reinterpret_cast<char*>(vertex_offsets.data()),
                                        (nodes_in_batch + 1) * sizeof(unsigned long long));
        }
        stream_reader.add_disk_read_time(io_seconds);
    }
    unsigned long long next_pos =
        stream_reader.binary_start_pos() + (nodes_in_batch) * sizeof(unsigned long long);

    unsigned long long edge_start_pos = vertex_offsets[0];
    unsigned long long num_reads = vertex_offsets[nodes_in_batch] - vertex_offsets[0];
    unsigned long long num_edges_to_read = num_reads / sizeof(unsigned long long);
    std::vector<unsigned long long> edges(num_edges_to_read);
    stream_reader.stream().seekg(edge_start_pos);
    {
        double io_seconds = 0.0;
        {
            TIMED_SCOPE(io_seconds, "edge_binary_read_edges");
            stream_reader.stream().read(reinterpret_cast<char*>(edges.data()),
                                        (num_edges_to_read) * sizeof(unsigned long long));
        }
        stream_reader.add_disk_read_time(io_seconds);
    }
    stream_reader.set_binary_start_pos(next_pos);

    unsigned long long pos = 0;
    for (unsigned long long i = 0; i < nodes_in_batch; ++i) {
        node = static_cast<NodeID>(i);
        unsigned long long degree =
            (vertex_offsets[i + 1] - vertex_offsets[i]) / sizeof(unsigned long long);
        G_temp_edges[node].reserve(degree);
        mapping[node] = node + pass_state.lower_global_node_conv;
        for (unsigned long long j = 0; j < degree; j++, pos++) {
            const LongNodeID target = static_cast<LongNodeID>(edges[pos]);
            add_split_edge(G_temp_edges, mapping, prev_batch_mapping, node, target_idx, map_index,
                           G_temp_no_nodes, G_temp_no_edges, target,
                           pass_state.lower_global_node_conv, pass_state.upper_global_node_conv);
        }
    }
    edge_partitioning::model::build_edge_split_batch_graph(
        partition_config, pass_state, G_temp_edges, mapping, rev_edge, contracted_edge_graph_id,
        number_of_deg1_vertices, number_of_deg2_vertices, G_temp_no_nodes, G_temp_no_edges, G_temp);
}

} // namespace io::edge_stream

