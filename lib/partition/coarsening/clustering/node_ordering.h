/******************************************************************************
 * node_ordering.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef NODE_ORDERING_HM1YMLB1
#define NODE_ORDERING_HM1YMLB1

#include <algorithm>

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "partition/partition_config.h"
#include "tools/random_functions.h"

class node_ordering {
   public:
    node_ordering();
    virtual ~node_ordering();

    void order_nodes(const Config& config, graph_access& G, std::vector<NodeID>& ordered_nodes) {
        for (unsigned int i = 0; i < ordered_nodes.size(); i++) {
            ordered_nodes[i] = i;
        }

        switch (config.node_ordering) {
            case RANDOM_NODEORDERING:
                order_nodes_random(config, G, ordered_nodes);
                break;
            case DEGREE_NODEORDERING:
                order_nodes_degree(config, G, ordered_nodes);
                break;
            case NATURAL_NODEORDERING:
                order_nodes_natural(config, G, ordered_nodes);
                break;
        }
    }

    void order_nodes_random(const Config& config, graph_access& G,
                            std::vector<NodeID>& ordered_nodes) {
        random_functions::permutate_vector_fast(ordered_nodes, false);
    }

    void order_nodes_degree(const Config& config, graph_access& G,
                            std::vector<NodeID>& ordered_nodes) {
        std::sort(ordered_nodes.begin(), ordered_nodes.end(),
                  [&](const NodeID& lhs, const NodeID& rhs) -> bool {
                      return (G.getNodeDegree(lhs) < G.getNodeDegree(rhs));
                  });
    }

    void order_nodes_natural(const Config& config, graph_access& G,
                             std::vector<NodeID>& ordered_nodes) {}
};

#endif /* end of include guard: NODE_ORDERING_HM1YMLB1 */
