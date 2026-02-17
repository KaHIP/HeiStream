/****************************************************************************
 * flatbuffer_log_writer.h
 *****************************************************************************/

#ifndef IO_OUTPUT_FLATBUFFER_LOG_WRITER_H_
#define IO_OUTPUT_FLATBUFFER_LOG_WRITER_H_

#include <string>

#include "partition/partition_config.h"

class FlatBufferWriter;


namespace io::output {

// Persist runtime/resource metrics to reporting sinks (console + optional bin log).
void write_flatbuffer_log(Config& partition_config, FlatBufferWriter& fb_writer,
                          const std::string& graph_filename, double read_graph_time,
                          double model_time, double postprocess_time,
                          double buffer_maintenance_time, double partition_time, double total_time,
                          long max_rss_kb);

} // namespace io::output


#endif /* IO_OUTPUT_FLATBUFFER_LOG_WRITER_H_ */
