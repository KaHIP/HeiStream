/******************************************************************************
 * node_fennel_assign.h
 *****************************************************************************/

#ifndef PARTITION_ASSIGNMENT_NODE_FENNEL_ASSIGN_H_
#define PARTITION_ASSIGNMENT_NODE_FENNEL_ASSIGN_H_

#include <algorithm>
#include <limits>
#include <vector>

#include "definitions.h"
#include "partition/heuristics/fennel_scoring.h"
#include "partition/partition_config.h"
#include "partition/state/node_partitioner_pass_state.h"

namespace node_fennel_assign {

inline Gain queue_key_from_weight(NodeWeight weight) {
    const NodeWeight capped = std::min<NodeWeight>(
        weight, static_cast<NodeWeight>(std::numeric_limits<Gain>::max()));
    return -static_cast<Gain>(capped);
}

inline void rebuild_block_min_queue(const Config& partition_config,
                                    partition::state::NodePartitionerPassState& pass_state) {
    auto& queue_ptr = pass_state.stream_block_min_queue;
    queue_ptr = std::make_unique<maxNodeHeap>();

    if (!pass_state.stream_blocks_weight) {
        return;
    }
    for (PartitionID p = 0; p < partition_config.k; ++p) {
        queue_ptr->insert(p, queue_key_from_weight((*pass_state.stream_blocks_weight)[p]));
    }
}

inline void ensure_block_min_queue_ready(const Config& partition_config,
                                         partition::state::NodePartitionerPassState& pass_state) {
    if (!partition_config.use_queue) {
        return;
    }
    if (!pass_state.stream_block_min_queue ||
        pass_state.stream_block_min_queue->size() != partition_config.k) {
        rebuild_block_min_queue(partition_config, pass_state);
    }
}

inline void set_block_weight_and_queue_key(Config& partition_config,
                                           partition::state::NodePartitionerPassState& pass_state,
                                           PartitionID block) {
    if (!partition_config.use_queue || !pass_state.stream_blocks_weight) {
        return;
    }
    ensure_block_min_queue_ready(partition_config, pass_state);
    const Gain key = queue_key_from_weight((*pass_state.stream_blocks_weight)[block]);
    if (!pass_state.stream_block_min_queue->contains(block)) {
        pass_state.stream_block_min_queue->insert(block, key);
    } else {
        pass_state.stream_block_min_queue->changeKey(block, key);
    }
}

inline void add_block_weight(Config& partition_config,
                             partition::state::NodePartitionerPassState& pass_state,
                             PartitionID block, NodeWeight delta) {
    (*pass_state.stream_blocks_weight)[block] += delta;
    set_block_weight_and_queue_key(partition_config, pass_state, block);
}

inline void sub_block_weight(Config& partition_config,
                             partition::state::NodePartitionerPassState& pass_state,
                             PartitionID block, NodeWeight delta) {
    (*pass_state.stream_blocks_weight)[block] -= delta;
    set_block_weight_and_queue_key(partition_config, pass_state, block);
}

}  // namespace node_fennel_assign

// Assigns one stream node to a block using the shared Fennel heuristic.
inline void partition_single_node(Config& partition_config,
                                  partition::state::NodePartitionerPassState& pass_state,
                                  LongNodeID global_node_id, std::vector<LongNodeID>& adjacents) {
    // Reuse scratch storage to avoid per-node allocations in the hot path.
    static thread_local std::vector<float> hash_map;
    static thread_local std::vector<PartitionID> touched_partitions;
    const PartitionID k = partition_config.k;
    if (hash_map.size() < static_cast<size_t>(k)) {
        hash_map.resize(static_cast<size_t>(k), 0.0f);
    }
    touched_partitions.clear();

    auto& stream_nodes_assign = *pass_state.stream_nodes_assign;
    for (LongNodeID adj_id : adjacents) {
        PartitionID adj_part = stream_nodes_assign[adj_id - 1];
        if (adj_part < k) {
            if (hash_map[adj_part] == 0.0f) {
                touched_partitions.push_back(adj_part);
            }
            hash_map[adj_part] += 1.0f;
        } else if (adj_part < 2 * k) {
            PartitionID ghost_part = adj_part - k;
            if (hash_map[ghost_part] == 0.0f) {
                touched_partitions.push_back(ghost_part);
            }
            // Ghost neighbors contribute with ghost_weight.
            hash_map[ghost_part] += partition_config.ghost_weight;
        }
    }

    auto& stream_blocks_weight = *pass_state.stream_blocks_weight;
    PartitionID best_partition = 0;
    if (partition_config.use_queue) {
        node_fennel_assign::ensure_block_min_queue_ready(partition_config, pass_state);
        const PartitionID min_block = pass_state.stream_block_min_queue->maxElement();
        best_partition = partition::heuristics::fennel::choose_stream_block(
            hash_map, touched_partitions, stream_blocks_weight, partition_config.max_block_weight,
            partition_config.fennel_alpha_gamma, true, min_block, false);
    } else {
        best_partition = partition::heuristics::fennel::choose_stream_block(
            hash_map, touched_partitions, stream_blocks_weight, partition_config.max_block_weight,
            partition_config.fennel_alpha_gamma, false);
    }

    // Assign node and account load.
    stream_nodes_assign[global_node_id - 1] = best_partition;
    node_fennel_assign::add_block_weight(partition_config, pass_state, best_partition, 1);

    if (partition_config.ghost_neighbors_enabled) {
        for (LongNodeID adj_id : adjacents) {
            PartitionID adj_part = stream_nodes_assign[adj_id - 1];
            if (adj_part == INVALID_PARTITION || (k <= adj_part && adj_part < 2 * k)) {
                stream_nodes_assign[adj_id - 1] = best_partition + k;
            }
        }
    }

    // Reset only touched entries for the next call.
    for (PartitionID part : touched_partitions) {
        hash_map[part] = 0.0f;
    }
}

#endif  // PARTITION_ASSIGNMENT_NODE_FENNEL_ASSIGN_H_
