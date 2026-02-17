/****************************************************************************
 * stream_reader.cpp
 *****************************************************************************/

#include "io/stream_reader.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "core/timing/scoped_stage_timer.h"

namespace {

bool has_ending(const std::string& value, const std::string& ending) {
    if (value.length() < ending.length()) {
        return false;
    }
    return value.compare(value.length() - ending.length(), ending.length(), ending) == 0;
}

void open_text_stream(io::stream::StreamCursor& cursor, const std::string& graph_filename) {
    if (cursor.stream_in.is_open()) {
        cursor.stream_in.close();
    }
    cursor.stream_in = std::ifstream(graph_filename.c_str());
    if (!cursor.stream_in) {
        std::cerr << "Error opening " << graph_filename << '\n';
        exit(1);
    }
    cursor.binary_input = false;
}

void open_binary_stream(io::stream::StreamCursor& cursor, const std::string& graph_filename) {
    if (cursor.stream_in.is_open()) {
        cursor.stream_in.close();
    }
    cursor.stream_in = std::ifstream(graph_filename.c_str(), std::ios::binary | std::ios::in);
    if (!cursor.stream_in) {
        std::cerr << "Error opening " << graph_filename << '\n';
        exit(1);
    }
    cursor.binary_input = true;
    cursor.binary_start_pos = 3 * sizeof(unsigned long long);
}

bool read_next_non_comment_line(std::ifstream& stream_in, std::string& line) {
    while (std::getline(stream_in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '%') {
            continue;
        }
        return true;
    }
    return false;
}

io::stream::StreamHeader parse_text_header(const std::string& line) {
    io::stream::StreamHeader header;
    std::stringstream ss(line);
    ss >> header.nodes;
    ss >> header.edges;
    if (!(ss >> header.edge_weight_flag)) {
        header.edge_weight_flag = 0;
    }
    header.binary_input = false;
    return header;
}

io::stream::StreamHeader parse_binary_header(std::ifstream& stream_in) {
    std::vector<unsigned long long> buffer(3, 0);
    stream_in.read(reinterpret_cast<char*>(&buffer[0]), 3 * sizeof(unsigned long long));
    unsigned long long version = buffer[0];
    (void)version;

    io::stream::StreamHeader header;
    header.nodes = static_cast<LongNodeID>(buffer[1]);
    header.edges = static_cast<LongEdgeID>(buffer[2]) / 2;
    header.edge_weight_flag = 0;
    header.binary_input = true;
    return header;
}

void parse_numeric_tokens_fast(const std::string& line, std::vector<LongNodeID>& numbers) {
    numbers.clear();

    const char* ptr = line.data();
    const char* end = ptr + line.size();
    while (ptr < end) {
        while (ptr < end && (*ptr < '0' || *ptr > '9')) {
            ++ptr;
        }
        if (ptr >= end) {
            break;
        }
        LongNodeID value = 0;
        while (ptr < end && *ptr >= '0' && *ptr <= '9') {
            value = value * 10 + static_cast<LongNodeID>(*ptr - '0');
            ++ptr;
        }
        numbers.push_back(value);
    }
}

}  // namespace


namespace io::stream {

StreamHeader StreamReader::read_header(const StreamSourceConfig& source_config) {
    const bool is_binary =
        source_config.allow_binary && (has_ending(source_config.graph_filename, ".bin") ||
                                       has_ending(source_config.graph_filename, ".parhip"));

    if (is_binary) {
        open_binary_stream(cursor_, source_config.graph_filename);
        TIMED_SCOPE(disk_read_time_, "stream_read_header_binary");
        StreamHeader header = parse_binary_header(cursor_.stream_in);
        return header;
    }

    open_text_stream(cursor_, source_config.graph_filename);
    std::string line;
    if (!read_next_non_comment_line(cursor_.stream_in, line)) {
        std::cerr << "Failed to read stream header from " << source_config.graph_filename << '\n';
        exit(1);
    }
    return parse_text_header(line);
}

bool StreamReader::read_next_parsed_line(std::vector<LongNodeID>& numbers) {
    // I/O stage includes both line fetch and tokenization.
    TIMED_SCOPE(disk_read_time_, "stream_read_next_parsed_line");
    if (!read_next_non_comment_line(cursor_.stream_in, parse_line_buffer_)) {
        return false;
    }
    // Fast path parser: avoid constructing temporary scanner state for each line.
    parse_numeric_tokens_fast(parse_line_buffer_, numbers);
    return true;
}

} // namespace io::stream

