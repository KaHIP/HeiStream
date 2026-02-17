/******************************************************************************
 * quality_metrics.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

#ifndef QUALITY_METRICS_10HC2I5M
#define QUALITY_METRICS_10HC2I5M

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"

class quality_metrics {
   public:
    quality_metrics();
    virtual ~quality_metrics();

    float getFennelWeight(Config& partition_config);
    EdgeWeight fennel_objective(Config& partition_config, graph_access& G, int* partition_map,
                                double fennel_gamma, double fennel_alpha);
    EdgeWeight fennel_objective(Config& partition_config, graph_access& G, double fennel_gamma,
                                double fennel_alpha);
    EdgeWeight ghost_edge_cut(const Config& config, graph_access& G);
    EdgeWeight edge_cut(graph_access& G);
    EdgeWeight edge_cut(graph_access& G, int* partition_map);
    EdgeWeight edge_cut(graph_access& G, PartitionID lhs, PartitionID rhs);
    EdgeWeight max_communication_volume(graph_access& G);
    EdgeWeight min_communication_volume(graph_access& G);
    EdgeWeight max_communication_volume(graph_access& G, int* partition_map);
    EdgeWeight total_communication_volume(graph_access& G);
    EdgeWeight objective(const Config& config, graph_access& G, int* partition_map);
    EdgeWeight edge_cut_connected(graph_access& G, int* partition_map);
    int boundary_nodes(graph_access& G);
    NodeWeight separator_weight(graph_access& G);
    double balance(graph_access& G);
    double balance_edges(graph_access& G);
    double balance_separator(graph_access& G);
    double edge_balance(graph_access& G, const std::vector<PartitionID>& edge_partition);

    EdgeWeight edge_cut_full_stream(const Config& config, graph_access& G,
                                    std::vector<std::vector<EdgeWeight>>& edges_virtualReal);
    double balance_full_stream(std::vector<NodeWeight>& stream_blocks_weight);
};

#endif /* end of include guard: QUALITY_METRICS_10HC2I5M */
