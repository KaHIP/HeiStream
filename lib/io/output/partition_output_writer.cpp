/****************************************************************************
 * partition_output_writer.cpp
 *****************************************************************************/

#include "io/output/partition_output_writer.h"

#include <fstream>
#include <iostream>


namespace io::output {

void write_partition_output(Config& partition_config,
                            const std::vector<PartitionID>& stream_nodes_assign) {
    std::ofstream out(partition_config.filename_output.c_str());

    for (LongNodeID node = 0; node < static_cast<LongNodeID>(stream_nodes_assign.size()); ++node) {
        out << stream_nodes_assign[node] << "\n";
    }

    out.close();
}

} // namespace io::output

