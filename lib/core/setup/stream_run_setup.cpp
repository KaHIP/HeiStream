/****************************************************************************
 * stream_run_setup.cpp
 *****************************************************************************/

#include "core/setup/stream_run_setup.h"

#include <fstream>
#include <iostream>

namespace core::setup {

bool prepare_framework_run(Config& config, const std::string& graph_filename) {
    if (graph_filename.empty()) {
        std::cerr << "Input graph file path is empty.\n";
        return false;
    }

    std::ifstream in(graph_filename.c_str(), std::ios::binary | std::ios::in);
    if (!in) {
        std::cerr << "Unable to open input graph file: " << graph_filename << '\n';
        return false;
    }

    // Framework-level run metadata shared by all algorithms.
    config.graph_filename = graph_filename;
    config.stream_input = true;
    return true;
}

} // namespace core::setup
