/******************************************************************************
 * edge_partitioning_quality_eval.cpp
 *****************************************************************************/

#include "partition/metrics/edge_partitioning_quality_eval.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "data_structure/buffered_map.h"
#include "data_structure/graph_access.h"



namespace partition::metrics::edge {

static bool hasEnding(std::string const& string, std::string const& ending) {
    // Local suffix helper for text/binary input dispatch.
    if (string.length() >= ending.length()) {
        return (0 == string.compare(string.length() - ending.length(), ending.length(), ending));
    }
    return false;
}

static void find_reverse_edges(graph_access& G, std::vector<EdgeID>& m_reverse_edge) {
    if (!m_reverse_edge.empty()) {
        return;
    }

    const EdgeID guard = std::numeric_limits<EdgeID>::max();
    m_reverse_edge.resize(G.number_of_edges(), guard);

    for (NodeID u = 0; u < G.number_of_nodes(); ++u) {
        for (EdgeID e_uv = G.get_first_edge(u); e_uv < G.get_first_invalid_edge(u); ++e_uv) {
            NodeID v = G.getEdgeTarget(e_uv);
            if (u < v) {
                assert(m_reverse_edge[e_uv] == guard);
            } else {
                assert(m_reverse_edge[e_uv] != guard);
                continue;
            }

            for (EdgeID e_vu = G.get_first_edge(v); e_vu < G.get_first_invalid_edge(v); ++e_vu) {
                if (G.getEdgeTarget(e_vu) == u) {
                    m_reverse_edge[e_uv] = e_vu;
                    m_reverse_edge[e_vu] = e_uv;
                    break;
                }
            }
            assert(m_reverse_edge[e_uv] != guard);
        }
    }
}

void evaluate_partition_batch(Config& config, const std::string& filename,
                              std::vector<PartitionID>* stream_nodes_assign_in, NodeID& vertex_cut,
                              NodeID& replicas, double& replication_factor, double& balance) {
    // Load the entire graph and compute per-node block incidence.
    graph_access G_temp;
    long nmbNodes;
    long nmbEdges;
    if (hasEnding(config.graph_filename, ".bin") || hasEnding(config.graph_filename, ".parhip")) {
        // Reconstruct graph structure from the binary CSR-like format.
        std::vector<unsigned long long> buffer(3, 0);
        std::ifstream filebin(config.graph_filename.c_str(), std::ios::binary | std::ios::in);
        if (filebin) {
            filebin.read((char*)(&buffer[0]), 3 * sizeof(unsigned long long));
        }

        unsigned long long version = buffer[0];
        (void)version;
        nmbNodes = static_cast<NodeID>(buffer[1]);
        nmbEdges = static_cast<NodeID>(buffer[2]) / 2;

        unsigned long long nodes_in_batch = static_cast<unsigned long long>(nmbNodes);
        std::vector<unsigned long long> vertex_offsets(nodes_in_batch + 1);
        filebin.seekg(3 * (sizeof(unsigned long long)));
        filebin.read(reinterpret_cast<char*>(vertex_offsets.data()),
                     (nodes_in_batch + 1) * sizeof(unsigned long long));

        unsigned long long edge_start_pos = vertex_offsets[0];
        unsigned long long num_reads = vertex_offsets[nodes_in_batch] - vertex_offsets[0];
        unsigned long long num_edges_to_read = num_reads / sizeof(unsigned long long);
        std::vector<unsigned long long> edges(num_edges_to_read);
        filebin.seekg(edge_start_pos);
        filebin.read(reinterpret_cast<char*>(edges.data()),
                     (num_edges_to_read) * sizeof(unsigned long long));

        NodeID node;
        G_temp.start_construction(nmbNodes, nmbEdges * 2);
        unsigned long long pos = 0;
        for (unsigned long long i = 0; i < nodes_in_batch; ++i) {
            node = static_cast<NodeID>(i);
            G_temp.new_node();
            G_temp.setPartitionIndex(node, 0);
            G_temp.setNodeWeight(node, 1);
            unsigned long long degree =
                (vertex_offsets[i + 1] - vertex_offsets[i]) / sizeof(unsigned long long);
            for (unsigned long long j = 0; j < degree; j++, pos++) {
                auto target = static_cast<NodeID>(edges[pos]);
                const EdgeID e = G_temp.new_edge(node, target);
                G_temp.setEdgeWeight(e, 1);
            }
        }
        G_temp.finish_construction();
        // vectors free themselves
        filebin.close();
    } else {
        // Text input path for METIS-style adjacency lists.
        std::vector<std::vector<LongNodeID>> input(1);
        std::vector<std::string> lines(1);
        LongNodeID node_counter = 0;
        std::ifstream in(filename.c_str());
        if (!in) {
            std::cerr << "Error opening " << filename << '\n';
        }
        int ew = 0;
        std::getline(in, lines[0]);
        while (lines[0][0] == '%')
            std::getline(in, lines[0]);

        std::stringstream ss(lines[0]);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;
        bool read_ew = false;
        bool read_nw = false;
        if (ew == 1) {
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

        G_temp.start_construction(nmbNodes, nmbEdges * 2);

        while (std::getline(in, lines[0])) {
            if (lines[0][0] == '%')
                continue;
            NodeID node = node_counter++;
            G_temp.new_node();
            G_temp.setPartitionIndex(node, 0);
            G_temp.setNodeWeight(node, 1);
            buffered_input ss2(&lines);
            ss2.simple_scan_line(input[0]);
            std::vector<LongNodeID>& line_numbers = input[0];
            LongNodeID col_counter = 0;

            NodeWeight weight = 1;
            if (read_nw) {
                weight = line_numbers[col_counter++];
                total_nodeweight += weight;
            }
            while (col_counter < line_numbers.size()) {
                target = line_numbers[col_counter++];
                target = target - 1;
                EdgeWeight edge_weight = 1;
                if (read_ew) {
                    edge_weight = line_numbers[col_counter++];
                }
                total_edgeweight += edge_weight;
                const EdgeID e = G_temp.new_edge(node, target);
                G_temp.setEdgeWeight(e, 1);
            }
            lines[0].clear();
            if (in.eof()) {
                break;
            }
        }
        G_temp.finish_construction();
    }

    // Use in-memory assignment state when available; otherwise read the
    // streamed output file (progress mode without O(m) assignment vector).
    std::vector<PartitionID> stream_nodes_assign_eval;
    std::vector<PartitionID>* stream_nodes_assign = stream_nodes_assign_in;
    if (stream_nodes_assign == nullptr) {
        stream_nodes_assign_eval = std::vector<PartitionID>(nmbEdges, INVALID_PARTITION);
        stream_nodes_assign = &stream_nodes_assign_eval;
        std::string line;
        std::ifstream part_file(config.filename_output);
        if (!part_file) {
            std::cerr << "Error opening partition ID file." << filename << '\n';
            exit(1);
        }
        for (LongEdgeID i = 0; i < nmbEdges; i++) {
            std::getline(part_file, line);
            if (line[0] == '%') {
                i--;
                continue;
            }
            (*stream_nodes_assign)[i] = (PartitionID)atol(line.c_str());
        }
    }

    // Evaluate replication and balance.
    std::vector<EdgeID> m_reverse_edge;
    find_reverse_edges(G_temp, m_reverse_edge);

    std::vector<PartitionID> edge_partition(G_temp.number_of_edges(), -1);
    NodeID edgeCount = 0;
    NodeID remaining_nodes = nmbNodes;
    NodeID batchSize = std::min(static_cast<NodeID>(config.batch_size), remaining_nodes);
    NodeID lower_node = 0;
    NodeID upper_node = lower_node + batchSize;

    while (remaining_nodes != 0) {
        // Rebuild edge-partition IDs in the same stream-batch order as writer.
        for (NodeID u = lower_node; u < upper_node; u++) {
            for (EdgeID e = G_temp.get_first_edge(u); e < G_temp.get_first_invalid_edge(u); ++e) {
                NodeID v = G_temp.getEdgeTarget(e);
                if (v >= upper_node) {
                    continue;
                }
                if (v < lower_node) {
                    edge_partition[e] = (*stream_nodes_assign)[edgeCount];
                    edge_partition[m_reverse_edge[e]] = (*stream_nodes_assign)[edgeCount];
                    edgeCount++;
                } else {
                    if (u < v) {
                        edge_partition[e] = (*stream_nodes_assign)[edgeCount];
                        edge_partition[m_reverse_edge[e]] = (*stream_nodes_assign)[edgeCount];
                        edgeCount++;
                    }
                }
            }
        }
        lower_node = upper_node;
        remaining_nodes = remaining_nodes - batchSize;
        batchSize = std::min(static_cast<NodeID>(config.batch_size), remaining_nodes);
        upper_node = lower_node + batchSize;
    }

    std::vector<NodeID> v_ei(config.k, 0);
    std::vector<NodeID> block_weights(config.k, 0);
    vertex_cut = 0;
    replicas = 0;

    for (NodeID u = 0; u < G_temp.number_of_nodes(); ++u) {
        if (G_temp.getNodeDegree(u) == 0) {
            continue;
        }

        std::vector<bool> counted(config.k);
        NodeID incident_part_count = 0;

        for (EdgeID e = G_temp.get_first_edge(u); e < G_temp.get_first_invalid_edge(u); ++e) {
            PartitionID part = edge_partition[e];
            if (u < G_temp.getEdgeTarget(e)) {
                block_weights[part]++;
            }
            if (!counted[part]) {
                counted[part] = true;
                v_ei[part]++;
                ++incident_part_count;
                ++replicas;
            }
        }
        --replicas;
        if (incident_part_count > 1) {
            vertex_cut++;
        }
    }

    double sum = 0;
    double max = -1;
    double total_weight = 0;
    for (int i = 0; i < v_ei.size(); i++) {
        sum += v_ei[i];
        if (block_weights[i] > max) {
            max = block_weights[i];
        }
        if (block_weights[i] > config.stream_total_upperbound) {
            std::cout << "Partition is imbalanced" << '\n';
        }
        total_weight += block_weights[i];
    }
    replication_factor = sum / (double)G_temp.number_of_nodes();
    double balance_part_weight = ceil(total_weight / (double)config.k);
    balance = max / balance_part_weight;
}

void evaluate_partition(Config& config, const std::string& filename,
                        std::vector<PartitionID>* stream_nodes_assign_in, NodeID& vertex_cut,
                        NodeID& replicas, double& replication_factor, double& balance) {
    // Stream through the graph and compute replication factor.
    std::unique_ptr<std::ifstream> in;
    std::string bin_ending(".bin");
    std::string parhip_ending(".parhip");
    if (hasEnding(config.graph_filename, bin_ending) ||
        hasEnding(config.graph_filename, parhip_ending)) {
        in = std::unique_ptr<std::ifstream>(
            new std::ifstream(config.graph_filename.c_str(), std::ios::binary | std::ios::in));
    } else {
        in = std::unique_ptr<std::ifstream>(new std::ifstream(config.graph_filename.c_str()));
    }

    if (!(*in)) {
        std::cerr << "Error opening " << filename << '\n';
    }

    std::vector<std::string> lines(1);
    std::string line;
    LongNodeID node_counter = 0;
    LongEdgeID nmbEdges;
    LongNodeID nmbNodes;
    int ew = 0;
    if (hasEnding(config.graph_filename, bin_ending) ||
        hasEnding(config.graph_filename, parhip_ending)) {
        std::vector<unsigned long long> buffer(3, 0);
        if ((*in)) {
            (*in).read((char*)(&buffer[0]), 3 * sizeof(unsigned long long));
        }
        unsigned long long version = buffer[0];
        (void)version;
        nmbNodes = static_cast<NodeID>(buffer[1]);
        nmbEdges = static_cast<NodeID>(buffer[2]) / 2;
    } else {
        std::getline(*in, lines[0]);
        while (lines[0][0] == '%')
            std::getline(*in, lines[0]);

        std::stringstream ss(lines[0]);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;
    }

    std::vector<PartitionID> stream_nodes_assign_eval(config.total_stream_edges, INVALID_PARTITION);
    std::vector<PartitionID>* stream_nodes_assign = stream_nodes_assign_in;
    if (stream_nodes_assign == nullptr) {
        stream_nodes_assign = &stream_nodes_assign_eval;
    }
    std::ifstream part_file(config.filename_output);
    if (!part_file) {
        std::cerr << "Error opening partition ID file." << filename << '\n';
        exit(1);
    }
    for (LongEdgeID i = 0; i < config.total_stream_edges; i++) {
        std::getline(part_file, line);
        if (line[0] == '%') {
            i--;
            continue;
        }
        stream_nodes_assign_eval[i] = (PartitionID)atol(line.c_str());
    }

    vertex_cut = 0;
    replicas = 0;
    double max = -1;
    double total_weight = 0;
    NodeID remaining_nodes = nmbNodes;
    NodeID batch_size = std::min(static_cast<NodeID>(config.batch_size), remaining_nodes);
    NodeID lower_node = 0;
    NodeID upper_node = lower_node + batch_size;
    NodeID node = 0;
    NodeID target = 0;
    EdgeID incremental_edge_counter = 0;
    unsigned long long start_pos = 3 * sizeof(unsigned long long);
    std::vector<std::vector<PartitionID>> part_ids_on_nodes(nmbNodes,
                                                            std::vector<PartitionID>(config.k, 0));

    if (hasEnding(config.graph_filename, bin_ending) ||
        hasEnding(config.graph_filename, parhip_ending)) {
        // Binary stream path: read one batch of offset/edge arrays at a time.
        while (remaining_nodes != 0) {
            unsigned long long nodes_in_batch = static_cast<unsigned long long>(batch_size);
            std::vector<unsigned long long> vertex_offsets(nodes_in_batch + 1);
            (*in).seekg(start_pos);
            (*in).read(reinterpret_cast<char*>(vertex_offsets.data()),
                       (nodes_in_batch + 1) * sizeof(unsigned long long));
            unsigned long long next_pos = start_pos + (nodes_in_batch) * sizeof(unsigned long long);

            unsigned long long edge_start_pos = vertex_offsets[0];
            unsigned long long num_reads = vertex_offsets[nodes_in_batch] - vertex_offsets[0];
            unsigned long long num_edges_to_read = num_reads / sizeof(unsigned long long);
            std::vector<unsigned long long> edges(num_edges_to_read);
            (*in).seekg(edge_start_pos);
            (*in).read(reinterpret_cast<char*>(edges.data()),
                       (num_edges_to_read) * sizeof(unsigned long long));
            start_pos = next_pos;

            unsigned long long pos = 0;
            for (unsigned long long i = 0; i < nodes_in_batch; ++i) {
                node = static_cast<NodeID>(i);
                unsigned long long degree =
                    (vertex_offsets[i + 1] - vertex_offsets[i]) / sizeof(unsigned long long);
                for (unsigned long long j = 0; j < degree; j++, pos++) {
                    target = static_cast<NodeID>(edges[pos]);
                    if (target >= upper_node) {
                        // edge to future batch
                        continue;
                    } else if (target < lower_node) {
                        // edge to past batch
                        part_ids_on_nodes[node + lower_node]
                                         [(*stream_nodes_assign)[incremental_edge_counter]]++;
                        part_ids_on_nodes[target]
                                         [(*stream_nodes_assign)[incremental_edge_counter]]++;
                        incremental_edge_counter++;
                    } else {
                        // edge to current batch
                        if (node + lower_node < target) {
                            part_ids_on_nodes[node + lower_node]
                                             [(*stream_nodes_assign)[incremental_edge_counter]]++;
                            part_ids_on_nodes[target]
                                             [(*stream_nodes_assign)[incremental_edge_counter]]++;
                            incremental_edge_counter++;
                        }
                    }
                }
            }
            if (incremental_edge_counter != 0) {
                lower_node = upper_node;
                remaining_nodes = remaining_nodes - batch_size;
                batch_size = std::min(static_cast<NodeID>(config.batch_size), remaining_nodes);
                upper_node = lower_node + batch_size;
            }
        }
    } else {
        // Text stream path: consume the same batch windows as partitioning.
        NodeID node_counter_text = 0;
        while (remaining_nodes != 0) {
            std::vector<LongNodeID> input;
            while (node_counter_text < batch_size) {
                std::getline(*in, lines[0]);
                if (lines[0][0] == '%') {
                    continue;
                }
                buffered_input ss2(&lines);
                ss2.simple_scan_line(input);
                lines[0].clear();

                LongNodeID col_counter = 0;
                node = node_counter_text;
                while (col_counter < input.size()) {
                    target = input[col_counter++] - 1;
                    if (target >= upper_node) {
                        // edge to future batch
                        continue;
                    } else if (target < lower_node) {
                        // edge to past batch
                        part_ids_on_nodes[node + lower_node]
                                         [(*stream_nodes_assign)[incremental_edge_counter]]++;
                        part_ids_on_nodes[target]
                                         [(*stream_nodes_assign)[incremental_edge_counter]]++;
                        incremental_edge_counter++;
                    } else {
                        // edge to current batch
                        if (node + lower_node < target) {
                            part_ids_on_nodes[node + lower_node]
                                             [(*stream_nodes_assign)[incremental_edge_counter]]++;
                            part_ids_on_nodes[target]
                                             [(*stream_nodes_assign)[incremental_edge_counter]]++;
                            incremental_edge_counter++;
                        }
                    }
                }
                node_counter_text++;
            }
            node_counter_text = 0;
            lower_node = upper_node;
            remaining_nodes = remaining_nodes - batch_size;
            batch_size = std::min(static_cast<NodeID>(config.batch_size), remaining_nodes);
            upper_node = lower_node + batch_size;
        }
    }

    std::vector<NodeID> v_ei(config.k, 0);
    std::vector<NodeID> block_weights(config.k, 0);

    for (auto& part_ids_on_node : part_ids_on_nodes) {
        // Count unique incident partitions for each node (replicas/vertex-cut).
        bool set = false;
        NodeID blocks_on_node_count = 0;
        PartitionID cur_part = 0;
        for (auto& j : part_ids_on_node) {
            if (j > 0) {
                block_weights[cur_part] += j;
                v_ei[cur_part]++;
                replicas++;
                set = true;
                blocks_on_node_count++;
            }
            cur_part++;
        }
        if (set) {
            --replicas;
        }
        if (blocks_on_node_count > 1) {
            ++vertex_cut;
        }
    }

    double sum = 0;
    for (int i = 0; i < v_ei.size(); i++) {
        sum += v_ei[i];
        if (block_weights[i] / 2 > max) {
            max = block_weights[i] / 2;
        }
        if (block_weights[i] / 2 > config.stream_total_upperbound) {
            std::cout << "Partition is imbalanced" << '\n';
        }
        total_weight += block_weights[i] / 2;
    }
    replication_factor = sum / (double)nmbNodes;
    double balance_part_weight = ceil(total_weight / (double)config.k);
    balance = max / balance_part_weight;
    // streams free themselves
}

} // namespace partition::metrics::edge

