/******************************************************************************
 * streamedgepartition.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <adilchhabra7@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <memory>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "graph_io_stream.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "tools/flat_buffer_writer.h"
#include <sys/resource.h>
#include "cpi/run_length_compression.hpp"

#include "data_structure/compression_vectors/CompressionDataStructure.h"
#include "data_structure/compression_vectors/RunLengthCompressionVector.h"
#include "data_structure/compression_vectors/BatchRunLengthCompression.h"

#define MIN(A, B) ((A) < (B)) ? (A) : (B)
#define MAX(A, B) ((A) > (B)) ? (A) : (B)

void config_multibfs_initial_partitioning(PartitionConfig &partition_config);

long getMaxRSS();

std::string extractBaseFilename(const std::string &fullPath);

int main(int argn, char **argv) {
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
    )" << std::endl;
    PartitionConfig partition_config;
    std::string graph_filename;
    timer processing_t, model_t, partition_t, mapping_t;
    double model_construction_time = 0;
    double partition_time = 0;
    double mapping_time = 0;
    EdgeWeight total_edge_cut = 0;
    quality_metrics qm;
    balance_configuration bc;

    bool is_graph_weighted = false;
    bool suppress_output = false;
    bool recursive = false;

    int ret_code =
            parse_parameters(argn, argv, partition_config, graph_filename,
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

    partition_config.LogDump(stdout);
    partition_config.graph_filename = graph_filename;
    partition_config.stream_input = true;
    auto *G = new graph_access();
    partition_config.global_cycle_iterations = partition_config.reps;
    partition_config.edge_partition = true;
    partition_config.use_queue = true;
    partition_config.initial_partitioning_repetitions=4;

    if(partition_config.edge_partition && partition_config.rle_length != -1) {
        std::cout << "HeiStreamEdge does not currently support run length compression." << std::endl;
        exit(0);
    }

    int &passes = partition_config.num_streams_passes;
    for (partition_config.restream_number = 0;
         partition_config.restream_number < passes;
         partition_config.restream_number++) {

        // load number of nodes and edges from the text input
        // and initialize global objects
        graph_io_stream::readFirstLineStreamEdge(partition_config, graph_filename,
                                                 total_edge_cut);

        // while we have not streamed all nodes of the input graph
        while (partition_config.remaining_stream_nodes) {
            // nmbNodes = buffer size (delta) or remaining nodes in final batch
            partition_config.nmbNodes =
                            MIN(partition_config.stream_buffer_len,
                                partition_config.remaining_stream_graph_nodes);

            // ****** Graph IO + Model Construction ******
            G->set_partition_count(partition_config.k);
            model_t.restart();
            graph_io_stream::constructBatchModel(partition_config, *G);
            model_construction_time += model_t.elapsed();

            // ****** partition batch ******
            partition_t.restart();
            graph_io_stream::countAssignedNodes(partition_config);
            graph_io_stream::prescribeBufferInbalance(partition_config);
            bool already_fully_partitioned = (partition_config.restream_vcycle &&
                                              partition_config.restream_number);

            bc.configurate_balance(partition_config, *G,
                                   already_fully_partitioned ||
                                   !partition_config.stream_initial_bisections);
            config_multibfs_initial_partitioning(partition_config);
            graph_partitioner partitioner;
            partitioner.perform_partitioning(partition_config, *G);
            partition_time += partition_t.elapsed();
            ofs.close();

            // ****** permanent assignment ******

            mapping_t.restart();
            graph_io_stream::generalizeStreamPartitionEdge(partition_config, *G);
            mapping_time += mapping_t.elapsed();
            if (partition_config.stream_output_progress) {
                graph_io_stream::writeStreamOutput(partition_config, *G);
            }

        }
        partition_config.lower_global_node = 0;
        partition_config.upper_global_node = 0;
    }

    double total_time = processing_t.elapsed();
    delete G;
    if (partition_config.stream_output_progress) {
        (*partition_config.stream_out).close();
    }
    long maxRSS = getMaxRSS();
    FlatBufferWriter fb_writer;

    if (!partition_config.benchmark) {
        NodeID vertexCut = 0;
        NodeID replicas = 0;
        double replication_factor = 0;
        double balance = 0;

        if (partition_config.stream_output_progress) {
            partition_config.evaluate_mode = true;
        } else {
            partition_config.evaluate_mode = false;
        }

        if (partition_config.light_evaluator) {
            graph_io_stream::streamEvaluateEdgePartition(partition_config, graph_filename,
                                                         vertexCut, replicas, replication_factor, balance);
        } else {
            graph_io_stream::streamEvaluatePartitionEdgeBatch(partition_config, graph_filename, vertexCut,
                                                              replicas, replication_factor, balance);
        }
        fb_writer.updateEdgePartitionResults(vertexCut, replicas, replication_factor, balance);
    }

    // write the partition to the disc
    std::stringstream filename;
    std::string baseFilename = extractBaseFilename(graph_filename);
    if (!partition_config.benchmark && !partition_config.stream_output_progress) {
        graph_io_stream::writePartitionStream(partition_config);
    }

    fb_writer.updateResourceConsumption(partition_config.read_graph_time,model_construction_time, mapping_time, partition_time, total_time, maxRSS);
    fb_writer.write(graph_filename, partition_config);

    return 0;
}

void config_multibfs_initial_partitioning(PartitionConfig &partition_config) {
    if (partition_config.initial_part_multi_bfs &&
        partition_config.curr_batch >= 2) {
        partition_config.initial_partitioning_type = INITIAL_PARTITIONING_MULTIBFS;
    }
}

long getMaxRSS() {
    struct rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // The maximum resident set size is in kilobytes
        return usage.ru_maxrss;
    } else {
        std::cerr << "Error getting resource usage information." << std::endl;
        // Return a sentinel value or handle the error in an appropriate way
        return -1;
    }
}

// Function to extract the base filename without path and extension
std::string extractBaseFilename(const std::string &fullPath) {
    size_t lastSlash = fullPath.find_last_of('/');
    size_t lastDot = fullPath.find_last_of('.');

    if (lastSlash != std::string::npos) {
        // Found a slash, extract the substring after the last slash
        return fullPath.substr(lastSlash + 1, lastDot - lastSlash - 1);
    } else {
        // No slash found, just extract the substring before the last dot
        return fullPath.substr(0, lastDot);
    }
}
