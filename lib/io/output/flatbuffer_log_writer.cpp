/****************************************************************************
 * flatbuffer_log_writer.cpp
 *****************************************************************************/

#include "io/output/flatbuffer_log_writer.h"

#include "tools/flat_buffer_writer.h"


namespace io::output {

void write_flatbuffer_log(Config& partition_config, FlatBufferWriter& fb_writer,
                          const std::string& graph_filename, double read_graph_time,
                          double model_time, double postprocess_time,
                          double buffer_maintenance_time, double partition_time, double total_time,
                          long max_rss_kb) {
    fb_writer.updateResourceConsumption(read_graph_time, model_time, postprocess_time,
                                        buffer_maintenance_time, partition_time, total_time,
                                        max_rss_kb);

    if (partition_config.evaluate) {
        fb_writer.print_summary(graph_filename, partition_config);
    }

    if (partition_config.write_log) {
        fb_writer.write(graph_filename, partition_config);
    }
}

} // namespace io::output

