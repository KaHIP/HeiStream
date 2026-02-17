/****************************************************************************
 * partition_output_writer.h
 *****************************************************************************/

#ifndef IO_OUTPUT_PARTITION_OUTPUT_WRITER_H_
#define IO_OUTPUT_PARTITION_OUTPUT_WRITER_H_

#include <vector>

#include "definitions.h"
#include "partition/partition_config.h"


namespace io::output {

// Write the final stream partition assignment file.
void write_partition_output(Config& partition_config,
                            const std::vector<PartitionID>& stream_nodes_assign);

} // namespace io::output


#endif /* IO_OUTPUT_PARTITION_OUTPUT_WRITER_H_ */
