/****************************************************************************
 * edge_partitioning_evaluator.cpp
 *****************************************************************************/

#include "algorithms/edge_partitioning/evaluation/edge_partitioning_evaluator.h"

#include "partition/metrics/edge_partitioning_quality_eval.h"
#include "tools/flat_buffer_writer.h"


namespace edge_partitioning::evaluation {

void EdgePartitioningEvaluator::evaluate(Config& config, const std::string& graph_filename,
                                         StreamContext& ctx, FlatBufferWriter& fb_writer) const {
    const auto exec_cfg = config.stream_execution_config();
    const auto edge_tuning_cfg = config.edge_tuning_config();
    if (exec_cfg.benchmark || !exec_cfg.evaluate) {
        return;
    }

    // In progress mode, assignments are streamed to disk during execution.
    // Close before reading for evaluation so metrics see the full partition.
    if (config.stream_output_progress && ctx.edge_pass_state.stream_out) {
        ctx.edge_pass_state.stream_out->flush();
        ctx.edge_pass_state.stream_out->close();
        ctx.edge_pass_state.stream_out.reset();
    }

    NodeID vertex_cut = 0;
    NodeID replicas = 0;
    double replication_factor = 0.0;
    double balance = 0.0;
    (void)edge_tuning_cfg;
    if (config.light_evaluator) {
        partition::metrics::edge::evaluate_partition(
            config, graph_filename, ctx.edge_pass_state.stream_nodes_assign.get(), vertex_cut,
            replicas, replication_factor, balance);
    } else {
        partition::metrics::edge::evaluate_partition_batch(
            config, graph_filename, ctx.edge_pass_state.stream_nodes_assign.get(), vertex_cut,
            replicas, replication_factor, balance);
    }
    fb_writer.updateEdgePartitionResults(vertex_cut, replicas, replication_factor, balance);
}

} // namespace edge_partitioning::evaluation

