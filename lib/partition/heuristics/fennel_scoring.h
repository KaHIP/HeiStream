/******************************************************************************
 * fennel_scoring.h
 *****************************************************************************/

#ifndef PARTITION_HEURISTICS_FENNEL_SCORING_H_
#define PARTITION_HEURISTICS_FENNEL_SCORING_H_

#include <limits>
#include <vector>

#include "definitions.h"
#include "random_functions.h"

namespace partition::heuristics::fennel {

// Computes the Fennel load penalty term for placing a node into a block.
inline double load_penalty(double fennel_alpha_gamma, NodeWeight block_load,
                           NodeWeight node_weight = 1) {
    return static_cast<double>(node_weight) * fennel_alpha_gamma *
           random_functions::approx_sqrt(block_load);
}

// Computes assignment score = edge_gain - fennel_weight * load_penalty.
inline double assignment_score(double edge_gain, double fennel_alpha_gamma, NodeWeight block_load,
                               NodeWeight node_weight = 1, double fennel_weight = 1.0) {
    return edge_gain -
           fennel_weight * load_penalty(fennel_alpha_gamma, block_load, node_weight);
}

// Returns the lightest feasible block. If none are feasible, returns 0.
inline PartitionID find_min_feasible_block(const std::vector<NodeWeight>& block_weights,
                                           NodeWeight max_block_weight) {
    PartitionID best = 0;
    NodeWeight min_weight = std::numeric_limits<NodeWeight>::max();
    for (PartitionID p = 0; p < static_cast<PartitionID>(block_weights.size()); ++p) {
        const NodeWeight w = block_weights[p];
        if (w < min_weight && w < max_block_weight) {
            min_weight = w;
            best = p;
        }
    }
    if (min_weight == std::numeric_limits<NodeWeight>::max()) {
        return 0;
    }
    return best;
}

// Chooses the best block for node-stream Fennel assignment.
// If use_queue_candidates is true, evaluate only touched neighbor blocks plus
// an optional queue-selected minimum-load block.
inline PartitionID choose_stream_block(const std::vector<float>& edge_gains,
                                       const std::vector<PartitionID>& touched_blocks,
                                       const std::vector<NodeWeight>& block_weights,
                                       NodeWeight max_block_weight, double fennel_alpha_gamma,
                                       bool use_queue_candidates,
                                       PartitionID queue_min_block = INVALID_PARTITION,
                                       bool require_queue_min_feasible = true) {
    const PartitionID k = static_cast<PartitionID>(block_weights.size());
    if (k == 0) {
        return 0;
    }

    PartitionID best_partition = 0;
    double best_score = std::numeric_limits<double>::lowest();

    if (use_queue_candidates) {
        for (const PartitionID p : touched_blocks) {
            if (p >= k) {
                continue;
            }
            if (block_weights[p] >= max_block_weight) {
                continue;
            }
            const double score = assignment_score(edge_gains[p], fennel_alpha_gamma, block_weights[p]);
            if (score > best_score || (score == best_score && p < best_partition)) {
                best_score = score;
                best_partition = p;
            }
        }

        if (queue_min_block < k &&
            (!require_queue_min_feasible || block_weights[queue_min_block] < max_block_weight)) {
            const double score = assignment_score(edge_gains[queue_min_block], fennel_alpha_gamma,
                                                  block_weights[queue_min_block]);
            if (score > best_score || (score == best_score && queue_min_block < best_partition)) {
                best_score = score;
                best_partition = queue_min_block;
            }
        }

        if (best_score == std::numeric_limits<double>::lowest()) {
            return 0;
        }
        return best_partition;
    }

    for (PartitionID p = 0; p < k; ++p) {
        if (block_weights[p] >= max_block_weight) {
            continue;
        }
        const double score = assignment_score(edge_gains[p], fennel_alpha_gamma, block_weights[p]);
        if (score > best_score || (score == best_score && p < best_partition)) {
            best_score = score;
            best_partition = p;
        }
    }
    return best_partition;
}

}  // namespace partition::heuristics::fennel

#endif  // PARTITION_HEURISTICS_FENNEL_SCORING_H_
