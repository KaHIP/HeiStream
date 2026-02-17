/****************************************************************************
 * edge_partitioning_reporter.h
 *****************************************************************************/

#ifndef EDGE_PARTITIONING_REPORTER_H_
#define EDGE_PARTITIONING_REPORTER_H_

#include <string>

#include "core/interfaces/algorithm_reporter.h"
#include "partition/partition_config.h"


namespace edge_partitioning::reporting {

class EdgePartitioningReporter : public core::interfaces::IAlgorithmReporter {
   public:
    void print_run_configuration(const Config& config,
                                 const std::string& graph_filename) const override;
    void report(Config& config, const std::string& graph_filename, StreamContext& ctx,
                FlatBufferWriter& fb_writer) const override;
};

} // namespace edge_partitioning::reporting


#endif /* EDGE_PARTITIONING_REPORTER_H_ */
