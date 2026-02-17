#ifndef BUFFER_NKJSAF9
#define BUFFER_NKJSAF9

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include "algorithms/node_partitioning/execution/batch_id_manager.h"
#include "core/timing/scoped_stage_timer.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "definitions.h"
#include "macros_assertions.h"
#include "partition/assignment/node_fennel_assign.h"
#include "partition/partition_config.h"
#include "partition/state/node_partitioner_pass_state.h"
#include "random_functions.h"

class Buffer {
   private:
    struct NeighborScanStats {
        unsigned assigned = 0;
        unsigned ghost_assigned = 0;
        unsigned in_buffer = 0;
    };

    Config& config;
    partition::state::NodePartitionerPassState& pass_state;
    bucket_pq pq;

    LongNodeID total_degree_sum;
    LongNodeID node_counter;
    double progress;
    float current_beta;

    TIMING_DECLARE(double update_adj_time) = 0.0;
    TIMING_DECLARE(double insert_time) = 0.0;
    TIMING_DECLARE(double extract_time) = 0.0;

    std::function<void(PartitionTask&&)> push_task_callback;

    float b_score_dmax;  /// Maximum buffer score dmax, degree factor of all degrees higher than
                         /// this is counted as 1.0 in the score calculation

    NeighborScanStats scan_neighbors(const std::vector<LongNodeID>& adjacents,
                                     bool count_sep_batch_marker, bool count_buffer_membership,
                                     bool split_ghost_neighbors) {
        NeighborScanStats stats;
        for (const LongNodeID global_adj_id : adjacents) {
            const PartitionID pid = (*pass_state.stream_nodes_assign)[global_adj_id - 1];
            if (pid != INVALID_PARTITION) {
                if (split_ghost_neighbors && pid >= config.k && pid < 2 * config.k) {
                    ++stats.ghost_assigned;
                } else {
                    ++stats.assigned;
                }
                continue;
            }

            if (count_sep_batch_marker && config.sep_batch_marker &&
                (*pass_state.stream_nodes_batch_marker)[global_adj_id - 1] != INVALID_PARTITION) {
                ++stats.assigned;
                continue;
            }

            if (count_buffer_membership && pq.contains(global_adj_id)) {
                ++stats.in_buffer;
            }
        }
        return stats;
    }

    void update_neighbours_priority_impl(std::vector<LongNodeID>& adjacents, bool part_adj_directly,
                                         bool is_active_ghost_neighbor, bool emit_partition_tasks,
                                         const char* timing_scope_label) {
        TIMED_SCOPE(update_adj_time, timing_scope_label);
        if (part_adj_directly) {
            part_adj_directly = config.part_adj_directly;
        }
        if (adjacents.empty()) {
            return;
        }

        for (LongNodeID adj_id : adjacents) {
            PQItem* adj_buffer_item_ptr = pq.findItem(adj_id);
            if (adj_buffer_item_ptr == nullptr) {
                continue;
            }

            PQItem& adj_buffer_item = *adj_buffer_item_ptr;
            auto& adj_adjacents = adj_buffer_item.get_adjacents();
            const unsigned adj_degree = adj_adjacents.size();
            adj_buffer_item.num_adj_partitioned++;

            if (part_adj_directly && adj_degree > 3 &&
                adj_degree == adj_buffer_item.num_adj_partitioned) {
                if (!emit_partition_tasks) {
                    pq.deleteNode(adj_id);
                    partition_single_node(config, pass_state, adj_id, adj_adjacents);
                    update_neighbours_priority_impl(adj_adjacents, false, false, false,
                                                    "buffer_update_neighbors");
                    completely_remove_node(adj_id);
                } else {
                    std::vector<LongNodeID> adj_adjacents_copy = std::move(adj_adjacents);
                    pq.deleteNode(adj_id);
                    completely_remove_node(adj_id);

                    bool adj_is_active_ghost_neighbor = false;
                    if (config.ghost_neighbors_enabled) {
                        adj_is_active_ghost_neighbor =
                            config.k <= (*pass_state.stream_nodes_assign)[adj_id - 1] &&
                            (*pass_state.stream_nodes_assign)[adj_id - 1] < 2 * config.k;
                    }

                    update_neighbours_priority_impl(adj_adjacents_copy, true,
                                                    adj_is_active_ghost_neighbor, true,
                                                    "buffer_update_neighbors_parallel");
                    (*pass_state.stream_nodes_assign)[adj_id - 1] = TO_BE_PARTITIONED;
                    PartitionTask task(
                        -1, std::vector<BatchNode>{{adj_id, std::move(adj_adjacents_copy)}});
                    if (push_task_callback) {
                        push_task_callback(std::move(task));
                    }
                }
            } else {
                const float updated_buffer_score =
                    calc_updated_buffer_score(adj_id, adj_buffer_item, is_active_ghost_neighbor);
                pq.increaseKey(adj_id, updated_buffer_score);
            }
        }
    }

   public:
    Buffer(Config& partition_config, partition::state::NodePartitionerPassState& pass_state_in,
           LongNodeID max_buffer_size)
        : config(partition_config),
          pass_state(pass_state_in),
          pq(partition_config,
             static_cast<unsigned>(std::floor(get_max_buffer_score(partition_config) *
                                              partition_config.bq_disc_factor)) +
                 1,
             partition_config.number_of_nodes, max_buffer_size, partition_config.bq_disc_factor) {
        b_score_dmax = std::min(static_cast<float>(config.d_max),
                                10000.0f);  // Cap d_max to avoid too large values
        current_beta = config.haa_beta;

        total_degree_sum = 0;
        node_counter = 0;
        progress = 0.0;
    }

    void set_push_task_callback(std::function<void(PartitionTask&&)> callback) {
        push_task_callback = std::move(callback);
    }

    static float get_max_buffer_score(const Config& cfg) {
        switch (cfg.buffer_score_type) {
            case BUFFER_SCORE_ANR:
            case BUFFER_SCORE_HAA:
                return std::max(1.0f, cfg.haa_theta);
            case BUFFER_SCORE_CMS:
            case BUFFER_SCORE_NSS:
            case BUFFER_SCORE_CBSQ:
            case BUFFER_SCORE_CBS:
            default:
                return 1.0f + cfg.cbs_theta;
        }
    }

    // Berechnet den Score für einen Knoten
    float calc_buffer_score(LongNodeID global_node_id, const std::vector<LongNodeID>& adjacents,
                            unsigned& cnt_adj_partitioned) {
        (void)global_node_id;
        if (adjacents.size() == 0) {
            return 0;
        }

        float degree = static_cast<float>(adjacents.size());
        float buffer_score;
        assert(cnt_adj_partitioned == 0);

        switch (config.buffer_score_type) {
            case BUFFER_SCORE_ANR:  // Assigned Neighbor Count
            {
                const NeighborScanStats stats =
                    scan_neighbors(adjacents, false, false, true /* split ghosts */);
                cnt_adj_partitioned = stats.assigned;
                buffer_score = cnt_adj_partitioned / degree;
                break;
            }
            case BUFFER_SCORE_HAA: {
                if (!config.ghost_neighbors_enabled) {
                    const NeighborScanStats stats =
                        scan_neighbors(adjacents, true, false, false /* split ghosts */);
                    cnt_adj_partitioned = stats.assigned;
                    float degree_factor = std::min(degree / b_score_dmax, 1.0f);

                    /// OPTIMIZATION for BETA = 2
                    // buffer_score = degree_factor * degree_factor
                    buffer_score = std::pow(degree_factor, current_beta) +
                                   (config.haa_theta * (1 - degree_factor)) *
                                       (static_cast<float>(cnt_adj_partitioned) / degree);

                } else {
                    const NeighborScanStats stats =
                        scan_neighbors(adjacents, true, false, true /* split ghosts */);
                    cnt_adj_partitioned = stats.assigned;
                    const unsigned cnt_adj_ghost = stats.ghost_assigned;

                    float degree_factor = std::min(degree / b_score_dmax, 1.0f);

                    /// OPTIMIZATION for BETA = 2
                    // buffer_score = degree_factor * degree_factor
                    buffer_score =
                        std::pow(degree_factor, current_beta) +
                        (config.haa_theta * (1 - degree_factor)) *
                            (static_cast<float>(cnt_adj_partitioned) / degree +
                             config.ghost_weight * (static_cast<float>(cnt_adj_ghost) /
                                                    degree));  /// Add ghost neighbors to the score
                }

                break;
            }

            case BUFFER_SCORE_CMS:  // Community - Majority Score
            {
                static thread_local std::vector<int> hash_map;
                static thread_local std::vector<PartitionID> touched_partitions;
                if (hash_map.size() < static_cast<size_t>(config.k)) {
                    hash_map.resize(static_cast<size_t>(config.k), 0);
                }
                touched_partitions.clear();

                int most_common_partition_count = 0;
                for (LongNodeID adj_id : adjacents) {
                    PartitionID adj_part = (*pass_state.stream_nodes_assign)[adj_id - 1];
                    if (adj_part < TO_BE_PARTITIONED - 10000) {
                        cnt_adj_partitioned++;
                        if (hash_map[adj_part] == 0) {
                            touched_partitions.push_back(adj_part);
                        }
                        hash_map[adj_part]++;
                        if (hash_map[adj_part] > most_common_partition_count) {
                            most_common_partition_count = hash_map[adj_part];
                            if (most_common_partition_count > degree / 2) {
                                break;
                            }
                        }
                    }
                }
                for (PartitionID p : touched_partitions) {
                    hash_map[p] = 0;
                }
                buffer_score = (float)most_common_partition_count /
                               degree;  // +  (float) degree / b_score_dmax;
                break;
            }

            case BUFFER_SCORE_NSS:  // Neighborhood Seen Score
            {
                const NeighborScanStats stats =
                    scan_neighbors(adjacents, false, true, false /* split ghosts */);
                cnt_adj_partitioned = stats.assigned;
                const unsigned cnt_adj_in_buffer = stats.in_buffer;
                buffer_score = (float)(cnt_adj_partitioned + cnt_adj_in_buffer) / degree;
                break;
            }

            case BUFFER_SCORE_CBSQ:  // Cuttana buffer score 2
            {
                const NeighborScanStats stats =
                    scan_neighbors(adjacents, false, false, true /* split ghosts */);
                cnt_adj_partitioned = stats.assigned;

                buffer_score = std::min(std::pow(degree / b_score_dmax, current_beta), 1.0f) +
                               config.cbs_theta * (float)cnt_adj_partitioned / degree;
                break;
            }

            case BUFFER_SCORE_CBS:  // Cuttana buffer score
            {
                const NeighborScanStats stats =
                    scan_neighbors(adjacents, false, false, true /* split ghosts */);
                cnt_adj_partitioned = stats.assigned;

                buffer_score = std::min(degree / b_score_dmax, 1.0f) +
                               config.cbs_theta * (float)cnt_adj_partitioned /
                                   degree;  // Range: [0, 3] (first term: [0, 1], second term: [0,
                                            // config.cbs_theta])
                // if (config.ghost_neighbors_enabled) {
                //     // Add ghost neighbors to the score
                //     buffer_score += config.cbs_theta * config.ghost_weight * (float)
                //     cnt_adj_ghost / degree;
                // }
                break;
            }
            default:  // Default to classic buffer score
                std::cerr << "Unknown buffer score type: " << config.buffer_score_type << '\n';
                exit(1);
        }

        return buffer_score;
    }

    float calc_updated_buffer_score(LongNodeID node_id, PQItem& buffer_item,
                                    bool is_active_ghost_neighbor = false) {
        const std::vector<LongNodeID>& adjacents = buffer_item.get_adjacents();
        float degree = static_cast<float>(adjacents.size());

        switch (config.buffer_score_type) {
            case BUFFER_SCORE_ANR:
                return buffer_item.buffer_score + 1.0 / degree;

            case BUFFER_SCORE_HAA: {
                if (!config.ghost_neighbors_enabled) {
                    float degree_factor = std::min(degree / b_score_dmax, 1.0f);
                    return buffer_item.buffer_score +
                           config.haa_theta * (1.0 - degree_factor) * (1.0 / degree);
                } else {
                    // If neighboring node from which the update is triggered is an active ghost
                    // neighbor, adjust the score update accordingly
                    float ghost_neighbor_factor =
                        is_active_ghost_neighbor ? (config.inv_ghost_weight) : 1.0f;
                    float degree_factor = std::min(degree / b_score_dmax, 1.0f);
                    return buffer_item.buffer_score + ghost_neighbor_factor * config.haa_theta *
                                                          (1.0 - degree_factor) * (1.0 / degree);
                }
            }

            case BUFFER_SCORE_CMS:
                buffer_item.num_adj_partitioned = 0;
                return calc_buffer_score(node_id, adjacents, buffer_item.num_adj_partitioned);

            case BUFFER_SCORE_NSS:
                buffer_item.num_adj_partitioned = 0;
                return calc_buffer_score(node_id, adjacents, buffer_item.num_adj_partitioned);

            case BUFFER_SCORE_CBS:
            case BUFFER_SCORE_CBSQ:
            default:
                return buffer_item.buffer_score + config.cbs_theta / degree;
        }
    }

    float get_avg_degree() {
        return static_cast<float>(total_degree_sum) / node_counter;
    }

    // Adds a node to the buffer or partition directly if buffer score is higher than max score
    bool addNode(LongNodeID global_node_id, std::vector<LongNodeID>& adjacents) {
        TIMED_SCOPE(insert_time, "buffer_insert");

        unsigned num_adj_partitioned = 0;
        float buffer_score = calc_buffer_score(global_node_id, adjacents, num_adj_partitioned);
        unsigned bucket_idx = pq.discretize_score(buffer_score);

        // if (pq.size() >= config.max_buffer_size && bucket_idx > pq.maxValue()) { // test1
        // if (pq.size() >= config.max_buffer_size && bucket_idx >= pq.maxValue()) { // test2
        if (config.batch_size == 1 && pq.size() >= config.max_buffer_size &&
            bucket_idx >= pq.maxValue()) {  // test3
            return false;
        } else {
            pq.insert(global_node_id, buffer_score, adjacents, num_adj_partitioned, bucket_idx);
            return true;
        }
    }

    // Removes and partitions the node with the highest priority
    void partitionTopNode() {
        LongNodeID node_id;
        {
            TIMED_SCOPE(extract_time, "buffer_extract_top");
            node_id = pq.deleteMax();
        }

        // Partition the node
        auto& adjacents = pq.getBufferItem(node_id).get_adjacents();
        partition_single_node(config, pass_state, node_id, adjacents);

        // Update neighbors and clear buffer item
        update_neighbours_priority(adjacents);
        completely_remove_node(node_id);
    }

    bool max_value_above_1(float val = 0.0f) {
        return pq.maxValue() > val * config.bq_disc_factor;
    }

    // Loads the top nodes into a batch for MLP processing
    void loadTopNodesToBatch(
        std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>& batch_nodes,
        LongNodeID batch_size, size_t batch_id) {
        // Initialize the partition configuration

        config.nmbNodes = std::min<LongNodeID>(batch_size, static_cast<LongNodeID>(pq.size()));
        if (batch_nodes.capacity() < static_cast<size_t>(config.nmbNodes)) {
            batch_nodes.reserve(static_cast<size_t>(config.nmbNodes));
        }
        // batch_nodes = new std::vector<std::pair<LongNodeID,
        // std::vector<LongNodeID>>>(config.nmbNodes);

        PartitionID batch_marker = config.batch_manager->get_batch_marker(batch_id);

        std::vector<PartitionID>& node_to_batch_marker =
            config.sep_batch_marker ? (*pass_state.stream_nodes_batch_marker)
                                    : (*pass_state.stream_nodes_assign);

        // Extract the top batch_size number of nodes from the queue
        LongNodeID local_node_counter = 0;
        while (local_node_counter < config.nmbNodes && !pq.empty()) {
            LongNodeID node_id;
            {
                TIMED_SCOPE(extract_time, "buffer_extract_batch");
                node_id = pq.deleteMax();
            }
            auto& adjacents = pq.getBufferItem(node_id).get_adjacents();

            node_to_batch_marker[node_id - 1] = batch_marker;
            update_neighbours_priority(adjacents, false);

            batch_nodes.emplace_back(node_id, std::move(adjacents));
            completely_remove_node(node_id);

            // if (config.batch_extraction_strategy ==
            // BATCH_EXTRACTION_STRATEGY_COMPLETE_BATCH_WITH_ADJ) {
            //     for (LongNodeID &adj_id : adjacents) {
            //         if (pq.contains(adj_id) && local_node_counter < config.nmbNodes - 2) {
            //             pq.deleteNode(adj_id);
            //             std::vector<LongNodeID> adj_adjacents = std::move(get_adjacents(adj_id));
            //             (*pass_state.stream_nodes_assign)[adj_id - 1] = batch_marker;

            //             update_neighbours_priority(adj_adjacents, false);
            //             completely_remove_node(adj_id);
            //             // batch_nodes.emplace_back(adj_id, std::move(adj_adjacents));
            //             (*batch_nodes)[local_node_counter] = std::make_pair(adj_id,
            //             std::move(adj_adjacents));

            //             local_node_counter++;
            //         }
            //     }
            // }

            local_node_counter++;
        }
    }

    // Loads the top nodes into a batch for MLP processing
    void loadTopNodesAndNeighborsToBatch(
        std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>& batch_nodes,
        LongNodeID batch_size, size_t batch_id) {
        // Initialize the partition configuration
        (void)batch_size;
        if (batch_nodes.capacity() < static_cast<size_t>(config.nmbNodes)) {
            batch_nodes.reserve(static_cast<size_t>(config.nmbNodes));
        }

        PartitionID batch_marker = config.batch_manager->get_batch_marker(batch_id);
        std::vector<PartitionID>& node_to_batch_marker =
            config.sep_batch_marker ? (*pass_state.stream_nodes_batch_marker)
                                    : (*pass_state.stream_nodes_assign);
        LongNodeID local_node_counter = 0;
        while (local_node_counter < config.nmbNodes && !pq.empty()) {
            LongNodeID node_id;
            {
                TIMED_SCOPE(extract_time, "buffer_extract_batch_neighbors");
                node_id = pq.deleteMax();
            }
            std::vector<LongNodeID> adjacents = std::move(get_adjacents(node_id));

            node_to_batch_marker[node_id - 1] = batch_marker;

            update_neighbours_priority(adjacents, false);
            completely_remove_node(node_id);

            if (config.batch_extraction_strategy ==
                BATCH_EXTRACTION_STRATEGY_COMPLETE_BATCH_WITH_ADJ) {
                for (LongNodeID& adj_id : adjacents) {
                    if (pq.contains(adj_id) && local_node_counter < config.nmbNodes - 2) {
                        pq.deleteNode(adj_id);
                        std::vector<LongNodeID> adj_adjacents = std::move(get_adjacents(adj_id));
                        node_to_batch_marker[adj_id - 1] = batch_marker;

                        update_neighbours_priority(adj_adjacents, false);
                        completely_remove_node(adj_id);
                        batch_nodes.emplace_back(adj_id, std::move(adj_adjacents));

                        local_node_counter++;
                    }
                }
            }

            batch_nodes.emplace_back(node_id, std::move(adjacents));

            local_node_counter++;
        }

        if (local_node_counter < batch_nodes.size()) {
            batch_nodes.resize(local_node_counter);
        }
    }

    // Update the priority value of the neighbours of the node that was just partitioned in the
    // priority queue
    void update_neighbours_priority(std::vector<LongNodeID>& adjacents,
                                    bool part_adj_directly = false,
                                    bool emit_partition_tasks = false,
                                    bool is_active_ghost_neighbor = false) {
        const char* label =
            emit_partition_tasks ? "buffer_update_neighbors_parallel" : "buffer_update_neighbors";
        update_neighbours_priority_impl(adjacents, part_adj_directly, is_active_ghost_neighbor,
                                        emit_partition_tasks, label);
    }

    double get_update_adj_time() {
        return update_adj_time;
    }
    double get_insert_time() {
        return insert_time;
    }
    double get_extract_time() {
        return extract_time;
    }

    // Helper-Methoden
    bool isEmpty() const {
        return pq.empty();
    }
    size_t size() const {
        return pq.size();
    }

    std::vector<LongNodeID>& get_adjacents(LongNodeID node_id) {
        return pq.getBufferItem(node_id).get_adjacents();
    }

    void completely_remove_node(LongNodeID node_id) {
        pq.completely_remove_node(node_id);
    }

    unsigned get_max_value() {
        return pq.maxValue();
    }

    LongNodeID deleteMax() {
        LongNodeID node_id;
        {
            TIMED_SCOPE(extract_time, "buffer_extract_top");
            node_id = pq.deleteMax();
        }
        return node_id;
    }
};

#endif  // BUFFER_NKJSAF9
