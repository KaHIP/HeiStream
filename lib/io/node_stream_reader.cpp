/****************************************************************************
 * stream_node_reader.cpp
 *****************************************************************************/

#include "io/node_stream_reader.h"

#include <utility>

#include "io/stream_reader.h"

namespace {

std::unique_ptr<std::vector<std::vector<LongNodeID>>> load_lines_from_stream_to_binary(
    io::stream::StreamReader& stream_reader, LongNodeID num_lines) {
    std::unique_ptr<std::vector<std::vector<LongNodeID>>> input(
        new std::vector<std::vector<LongNodeID>>(num_lines));
    LongNodeID node_counter = 0;
    while (node_counter < num_lines) {
        if (!stream_reader.read_next_parsed_line((*input)[node_counter])) {
            break;
        }
        ++node_counter;
    }
    input->resize(static_cast<size_t>(node_counter));
    return input;
}

}  // namespace

NodeStreamReader::NodeStreamReader(Config& config)
    : config_(config), stream_reader_(), ram_index_(0) {}

void NodeStreamReader::begin_pass(const std::string& graph_filename) {
    io::stream::StreamSourceConfig source_config;
    source_config.graph_filename = graph_filename;
    source_config.allow_binary = false;
    header_ = stream_reader_.read_header(source_config);
    ram_index_ = 0;
    ram_lines_.reset();
}

bool NodeStreamReader::next_node_line(std::vector<LongNodeID>& line_numbers) {
    line_numbers.clear();

    if (config_.ram_stream && !ram_lines_) {
        ram_lines_ = load_lines_from_stream_to_binary(stream_reader_, header_.nodes);
    }

    if (ram_lines_) {
        if (ram_index_ >= ram_lines_->size()) {
            return false;
        }
        line_numbers = std::move((*ram_lines_)[ram_index_++]);
        return true;
    }

    return stream_reader_.read_next_parsed_line(line_numbers);
}

bool NodeStreamReader::next_batch(LongNodeID num_lines,
                                  std::vector<std::vector<LongNodeID>>& batch_lines) {
    batch_lines.clear();
    batch_lines.reserve(static_cast<size_t>(num_lines));

    std::vector<LongNodeID> line_numbers;
    line_numbers.reserve(64);
    for (LongNodeID i = 0; i < num_lines; ++i) {
        if (!next_node_line(line_numbers)) {
            break;
        }
        batch_lines.push_back(std::move(line_numbers));
        line_numbers.clear();
    }
    return !batch_lines.empty();
}

void NodeStreamReader::end_pass() {
    ram_lines_.reset();
    ram_index_ = 0;
}
