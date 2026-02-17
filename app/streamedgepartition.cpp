/******************************************************************************
 * streamedgepartition.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <adilchhabra7@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <regex.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "algorithms/bootstrap/register_builtin_algorithms.h"
#include "cli/parse_parameters.h"
#include "core/pipeline/stream_engine.h"
#include "partition/partition_config.h"
#include "random_functions.h"

int main(int argn, char** argv) {
    std::cout << R"(
██   ██ ███████ ██ ███████ ████████ ██████  ███████  █████  ███    ███
██   ██ ██      ██ ██         ██    ██   ██ ██      ██   ██ ████  ████
███████ █████   ██ ███████    ██    ██████  █████   ███████ ██ ████ ██
██   ██ ██      ██      ██    ██    ██   ██ ██      ██   ██ ██  ██  ██
██   ██ ███████ ██ ███████    ██    ██   ██ ███████ ██   ██ ██      ██


███████ ██████   ██████  ███████
██      ██   ██ ██       ██
█████   ██   ██ ██   ███ █████
██      ██   ██ ██    ██ ██
███████ ██████   ██████  ███████
    )" << '\n';
    Config partition_config;
    std::string graph_filename;
    bool is_graph_weighted = false;
    bool suppress_output = false;
    bool recursive = false;

    int ret_code = parse_parameters_edge(argn, argv, partition_config, graph_filename,
                                         is_graph_weighted, suppress_output, recursive);

    if (ret_code) {
        return 0;
    }

    std::ofstream ofs;
    ofs.open("/dev/null");
    if (suppress_output) {
        std::cout.rdbuf(ofs.rdbuf());
    }
    srand(partition_config.seed);
    random_functions::setSeed(partition_config.seed);
    algorithms::bootstrap::register_builtin_algorithms();

    partition_config.LogDump(stdout);
    core::pipeline::run_stream_engine(StreamMode::Edge, partition_config, graph_filename);

    return 0;
}
