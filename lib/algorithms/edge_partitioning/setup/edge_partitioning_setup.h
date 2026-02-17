/****************************************************************************
 * edge_partitioning_setup.h
 *****************************************************************************/

#ifndef EDGE_PARTITIONING_SETUP_H_
#define EDGE_PARTITIONING_SETUP_H_

#include <string>

#include "core/interfaces/algorithm_setup.h"
#include "io/stream_reader.h"
#include "partition/partition_config.h"
#include "partition/state/edge_partitioner_pass_state.h"


namespace edge_partitioning::setup {

class EdgePartitioningSetup : public core::interfaces::IAlgorithmSetup {
   public:
    void prepare_run(Config& config, const std::string& graph_filename) const override;
    void initialize_pass(Config& partition_config,
                         partition::state::EdgePartitionerPassState& pass_state,
                         const io::stream::StreamHeader& header) const;
};

} // namespace edge_partitioning::setup


#endif /* EDGE_PARTITIONING_SETUP_H_ */
