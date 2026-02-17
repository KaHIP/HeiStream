/****************************************************************************
 * node_partitioning_quality_eval.cpp
 *****************************************************************************/

#include "partition/metrics/node_partitioning_quality_eval.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "data_structure/buffered_map.h"



namespace partition::metrics::node {

void evaluate_partition(const Config& config, const std::vector<PartitionID>& stream_nodes_assign,
                        const std::string& filename, EdgeWeight& edgeCut) {
    // Stream through the input graph and count cross-partition edges.
    std::vector<std::vector<LongNodeID>> input(1);
    std::vector<std::string> lines(1);
    LongNodeID node_counter = 0;
    std::ifstream in(filename.c_str());
    if (!in) {
        std::cerr << "Error opening " << filename << '\n';
        return;
    }
    long nmbNodes;
    long nmbEdges;
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
    edgeCut = 0;

    while (std::getline(in, lines[0])) {
        if (lines[0][0] == '%')
            continue;
        NodeID node = node_counter++;
        PartitionID partitionIDSource = stream_nodes_assign[node];

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
            PartitionID partitionIDTarget = stream_nodes_assign[target];
            if (partitionIDSource != partitionIDTarget) {
                edgeCut += edge_weight;
            }
        }
        lines[0].clear();
        if (in.eof()) {
            break;
        }
    }
    edgeCut = edgeCut / 2;  // Edges are counted twice.
}

} // namespace partition::metrics::node


