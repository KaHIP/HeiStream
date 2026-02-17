/****************************************************************************
 * model_partition_mapper.cpp
 *****************************************************************************/

#include "partition/model/model_partition_mapper.h"

#include <cstdlib>
#include <iostream>

namespace partition_model {

PartitionID recover_partition_for_model_node(const Config& config,
                                             const std::vector<PartitionID>* stream_nodes_assign,
                                             LongNodeID lower_global_node, NodeID node,
                                             NodeID node_counter, int restream_number) {
    if (node >= node_counter - config.quotient_nodes) {
        return node - (node_counter - config.quotient_nodes);
    }

    if (restream_number &&
        (config.restream_vcycle || config.initial_partitioning_type == INITIAL_PARTITIONING_FENNEL)) {
        if (stream_nodes_assign == nullptr) {
            std::cerr << "Missing stream assignment state for restream model recovery." << '\n';
            std::exit(1);
        }
        if (node < config.nmbNodes) {
            const LongNodeID global_node_id =
                config.local_to_global_map ? (*config.local_to_global_map)[node]
                                           : lower_global_node + static_cast<LongNodeID>(node);
            return (*stream_nodes_assign)[global_node_id - 1];
        }
        std::cerr << "Unexpected model-node partition recovery branch." << '\n';
        std::exit(1);
    }

    return 0;
}

}  // namespace partition_model
