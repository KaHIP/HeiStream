/****************************************************************************
 * stream_node_reader.h
 *****************************************************************************/

#ifndef STREAM_NODE_READER_H_
#define STREAM_NODE_READER_H_

#include <memory>
#include <string>
#include <vector>

#include "io/stream_reader.h"
#include "partition/partition_config.h"

class NodeStreamReader {
   public:
    /**
     * Construct a node stream reader bound to a run configuration.
     *
     * @param config Mutable run configuration.
     */
    explicit NodeStreamReader(Config& config);

    /**
     * Initialize one stream pass.
     *
     * @param graph_filename Input graph path.
     */
    void begin_pass(const std::string& graph_filename);

    /**
     * Read one logical node line.
     *
     * @param line_numbers Output parsed line contents.
     * @return `true` on successful read, `false` on EOF.
     */
    bool next_node_line(std::vector<LongNodeID>& line_numbers);

    /**
     * Read up to `num_lines` logical node lines.
     *
     * @param num_lines Maximum number of lines to read.
     * @param batch_lines Output line buffer.
     * @return `true` if at least one line was read, otherwise `false`.
     */
    bool next_batch(LongNodeID num_lines, std::vector<std::vector<LongNodeID>>& batch_lines);

    /**
     * Release pass-local buffered state.
     */
    void end_pass();

    /**
     * @return Most recently parsed stream header.
     */
    const io::stream::StreamHeader& header() const {
        return header_;
    }
    /**
     * @return Mutable configuration reference.
     */
    Config& config_ref() {
        return config_;
    }
    /**
     * @return Accumulated read+parse time in seconds.
     */
    double disk_read_time() const {
        return stream_reader_.disk_read_time();
    }

   private:
    Config& config_;
    io::stream::StreamReader stream_reader_;
    io::stream::StreamHeader header_;
    std::unique_ptr<std::vector<std::vector<LongNodeID>>> ram_lines_;
    size_t ram_index_;
};

#endif  // STREAM_NODE_READER_H_
