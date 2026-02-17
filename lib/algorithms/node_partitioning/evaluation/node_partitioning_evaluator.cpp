/****************************************************************************
 * node_partitioning_evaluator.cpp
 *****************************************************************************/

#include "algorithms/node_partitioning/evaluation/node_partitioning_evaluator.h"

#include "partition/metrics/node_partitioning_quality_eval.h"
#include "quality_metrics.h"
#include "tools/flat_buffer_writer.h"


namespace node_partitioning::evaluation {

void NodePartitioningEvaluator::evaluate(Config& config, const std::string& graph_filename,
                                         StreamContext& ctx, FlatBufferWriter& fb_writer) const {
    (void)ctx;
    const auto exec_cfg = config.stream_execution_config();
    if (!exec_cfg.evaluate) {
        return;
    }

    EdgeWeight edge_cut = 0;
    partition::metrics::node::evaluate_partition(config, *ctx.node_pass_state.stream_nodes_assign,
                                                 graph_filename, edge_cut);
    quality_metrics qm;
    const double imbalance = qm.balance_full_stream(*ctx.node_pass_state.stream_blocks_weight);
    fb_writer.updateVertexPartitionResults(edge_cut, imbalance);
}

} // namespace node_partitioning::evaluation

