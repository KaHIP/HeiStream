/****************************************************************************
 * node_partitioning_setup.h
 *****************************************************************************/

#ifndef NODE_PARTITIONING_SETUP_H_
#define NODE_PARTITIONING_SETUP_H_

#include <string>

#include "core/interfaces/algorithm_setup.h"
#include "io/stream_reader.h"
#include "partition/partition_config.h"
#include "partition/state/node_partitioner_pass_state.h"

namespace node_partitioning {
namespace setup {

class NodePartitioningSetup : public core::interfaces::IAlgorithmSetup {
   public:
    void prepare_run(Config& config, const std::string& graph_filename) const override;
    void initialize_pass(Config& partition_config,
                         partition::state::NodePartitionerPassState& pass_state,
                         const io::stream::StreamHeader& header) const;
};

}  // namespace setup
}  // namespace node_partitioning

#endif /* NODE_PARTITIONING_SETUP_H_ */
