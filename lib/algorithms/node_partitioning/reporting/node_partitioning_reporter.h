/****************************************************************************
 * node_partitioning_reporter.h
 *****************************************************************************/

#ifndef NODE_PARTITIONING_REPORTER_H_
#define NODE_PARTITIONING_REPORTER_H_

#include <cstddef>
#include <string>

#include "core/interfaces/algorithm_reporter.h"

class FlatBufferWriter;

namespace node_partitioning {
namespace reporting {

struct ParallelRuntimeTableData {
    double total_time = 0.0;
    double first_phase_time = 0.0;
    double second_phase_time = 0.0;
    double io_time = 0.0;
    double buffer_add_node_time = 0.0;
    double buffer_extract_time = 0.0;
    double updating_adj_time = 0.0;
    double part_single_node_time = 0.0;
    double mlp_time = 0.0;
    double io_thread_time = 0.0;
    double pq_thread_time = 0.0;
    double partition_thread_time = 0.0;
    double io_wait_time = 0.0;
    double pq_wait_time = 0.0;
    double partition_wait_time = 0.0;
    std::size_t nodes_processed_io = 0;
    std::size_t nodes_processed_pq = 0;
    std::size_t tasks_processed_partition = 0;
    std::size_t max_input_queue_size = 0;
    std::size_t max_partition_queue_size = 0;
};

void print_parallel_runtime_table(const ParallelRuntimeTableData& data);

class NodePartitioningReporter : public core::interfaces::IAlgorithmReporter {
   public:
    void print_run_configuration(const Config& config,
                                 const std::string& graph_filename) const override;
    void report(Config& config, const std::string& graph_filename, StreamContext& ctx,
                FlatBufferWriter& fb_writer) const override;
};

}  // namespace reporting
}  // namespace node_partitioning

#endif /* NODE_PARTITIONING_REPORTER_H_ */
