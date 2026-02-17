/****************************************************************************
 * stream_node_mode_internal.h
 *****************************************************************************/

#ifndef STREAM_NODE_MODE_INTERNAL_H_
#define STREAM_NODE_MODE_INTERNAL_H_

#include <string>
#include <utility>
#include <vector>

#include "algorithms/node_partitioning/execution/batch_id_manager.h"
#include "core/context/stream_context.h"
#include "partition/partition_config.h"
#include "partition/state/node_partitioner_pass_state.h"

class NodeStreamReader;
class Buffer;

namespace stream_node_modes {

bool use_mlp_batching(const Config& config);
double begin_node_stream_pass(NodeStreamReader& reader,
                              partition::state::NodePartitionerPassState& pass_state,
                              const std::string& graph_filename);
bool read_next_node_from_stream(NodeStreamReader& reader,
                                partition::state::NodePartitionerPassState& pass_state,
                                std::vector<LongNodeID>& neighbors, LongNodeID& global_node_id,
                                double& io_delta_seconds);
void ensure_priority_buffer_supports_input(
    const partition::state::NodePartitionerPassState& pass_state);
void write_stream_progress_if_enabled(Config& config,
                                      const partition::state::NodePartitionerPassState& pass_state);
bool is_active_ghost_neighbor(const Config& config,
                              const partition::state::NodePartitionerPassState& pass_state,
                              LongNodeID node_id);
bool should_skip_restream_node(bool include_high_degree_nodes, LongNodeID d_max,
                               const std::vector<LongNodeID>& adjacents);
void mark_node_batch_assignment(const Config& config,
                                partition::state::NodePartitionerPassState& pass_state,
                                LongNodeID node_id, PartitionID batch_marker);
BatchNode extract_top_buffer_node_for_batch(
    Config& config, partition::state::NodePartitionerPassState& pass_state, Buffer& buffer,
    PartitionID batch_marker, bool use_parallel_update);
void advance_batch_window(Config& config, size_t& current_batch_id,
                          PartitionID& current_batch_marker);
void acquire_batch_window(Config& config, size_t& current_batch_id,
                          PartitionID& current_batch_marker);
bool build_batch_task_and_advance(Config& config, std::vector<BatchNode>& current_batch,
                                  size_t& current_batch_id, PartitionID& current_batch_marker,
                                  PartitionTask& task);
void collect_buffer_maintenance_totals(Buffer& buffer, double& add_node_time, double& extract_time,
                                       double& update_adj_time);
void assign_single_node_direct(Config& config,
                               partition::state::NodePartitionerPassState& pass_state,
                               LongNodeID global_node_id, std::vector<LongNodeID>& adjacents);
void extract_fennel_neighbors(int remaining_stream_ew, const std::vector<LongNodeID>& line_numbers,
                              LongNodeID current_global_node, NodeWeight& node_weight,
                              std::vector<LongNodeID>& adjacents, LongEdgeID& used_edges);
void perform_mlp_on_batch(Config& partition_config,
                          partition::state::NodePartitionerPassState& pass_state,
                          std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>& batch_nodes,
                          size_t batch_id, double& model_time, double& partition_time,
                          double& postprocess_time);

void run_node_stream_priority_buffer_parallel(
    Config& config, const std::string& graph_filename, StreamContext& ctx,
    partition::state::NodePartitionerPassState& pass_state);

void run_node_stream_priority_buffer_sequential(
    Config& config, const std::string& graph_filename, StreamContext& ctx,
    partition::state::NodePartitionerPassState& pass_state);

void run_node_stream_no_priority_buffer(Config& config, const std::string& graph_filename,
                                        StreamContext& ctx,
                                        partition::state::NodePartitionerPassState& pass_state);

void run_node_stream_direct_fennel(Config& config, const std::string& graph_filename,
                                   StreamContext& ctx,
                                   partition::state::NodePartitionerPassState& pass_state);

}  // namespace stream_node_modes

#endif  // STREAM_NODE_MODE_INTERNAL_H_
