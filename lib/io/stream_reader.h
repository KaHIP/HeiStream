/****************************************************************************
 * stream_reader.h
 *****************************************************************************/

#ifndef IO_STREAM_READER_H_
#define IO_STREAM_READER_H_

#include <fstream>
#include <string>
#include <vector>

#include "definitions.h"


namespace io::stream {

/**
 * Source options for opening a streamed graph input.
 */
struct StreamSourceConfig {
    std::string graph_filename;
    bool allow_binary = false;
};

/**
 * Parsed header information for a stream input source.
 */
struct StreamHeader {
    LongNodeID nodes = 0;
    LongEdgeID edges = 0;
    int edge_weight_flag = 0;
    bool binary_input = false;
};

struct StreamCursor {
    std::ifstream stream_in;
    bool binary_input = false;
    unsigned long long binary_start_pos = 3 * sizeof(unsigned long long);
};

/**
 * Stateful reader for streamed graph input (text or binary).
 */
class StreamReader {
   public:
    StreamReader() = default;

    /**
     * Open stream input and parse the first header line.
     *
     * @param source_config Stream source options and path.
     * @return Parsed stream header.
     */
    StreamHeader read_header(const StreamSourceConfig& source_config);

    /**
     * Read next non-comment line and parse integer tokens.
     *
     * @param numbers Output buffer filled with parsed integer fields.
     * @return `true` on successful read, `false` on EOF.
     */
    bool read_next_parsed_line(std::vector<LongNodeID>& numbers);

    /**
     * @return Current binary offset cursor for block reads.
     */
    unsigned long long binary_start_pos() const {
        return cursor_.binary_start_pos;
    }
    /**
     * Set binary offset cursor for subsequent block reads.
     *
     * @param start_pos Byte offset into the open binary file.
     */
    void set_binary_start_pos(unsigned long long start_pos) {
        cursor_.binary_start_pos = start_pos;
    }
    /**
     * @return Underlying file stream.
     */
    std::ifstream& stream() {
        return cursor_.stream_in;
    }
    /**
     * @return `true` when current source is binary.
     */
    bool is_binary_input() const {
        return cursor_.binary_input;
    }
    /**
     * @return Accumulated read+parse time in seconds.
     */
    double disk_read_time() const {
        return disk_read_time_;
    }
    /**
     * Add measured read+parse time.
     *
     * @param seconds Time to accumulate.
     */
    void add_disk_read_time(double seconds) {
        disk_read_time_ += seconds;
    }

   private:
    StreamCursor cursor_;
    double disk_read_time_ = 0.0;
    // Reused buffer to avoid per-line string allocations in hot streaming loops.
    std::string parse_line_buffer_;
};

} // namespace io::stream


#endif /* IO_STREAM_READER_H_ */
