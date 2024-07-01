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

#include "graph_io_stream.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "flatbuffers/flatbuffers.h"
#include "PartitionInfo_generated.h"

#define MIN(A, B) ((A) < (B)) ? (A) : (B)
#define MAX(A, B) ((A) > (B)) ? (A) : (B)

void config_multibfs_initial_partitioning(PartitionConfig &partition_config);
std::string extractBaseFilename(const std::string& fullPath);

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
    std::cout << "Running evaluator..." << std::endl;

    PartitionConfig partition_config;
    std::string graph_filename;
    EdgeWeight total_edge_cut = 0;
    quality_metrics qm;


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
    partition_config.stream_input = true;
    partition_config.evaluate_mode = true;
    partition_config.graph_filename = graph_filename;

    std::string baseFilename = extractBaseFilename(graph_filename);
    std::string outputFileNameStream;
    outputFileNameStream = baseFilename + "_" + std::to_string(partition_config.k) + "_" + std::to_string(partition_config.stream_buffer_len) + ".bin";

    std::fstream file(outputFileNameStream, std::ios::binary | std::ios::in | std::ios::out);

    if (!file.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    // Read the contents of the file into a vector<char>
    file.seekg(0, std::ios::end);
    size_t length = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(length);
    file.read(buffer.data(), length);

    // Deserialize the FlatBuffer from the buffer
    flatbuffers::Verifier verifier(reinterpret_cast<const uint8_t*>(buffer.data()), length);
    if (!PartitionInfo::VerifyEdgePartitionBuffer(verifier)) {
        std::cerr << "Error verifying FlatBuffer." << std::endl;
        return 1;
    }

    // Access the root table
    auto edgePartition = PartitionInfo::GetEdgePartition(buffer.data());

    // Access individual fields of the EdgePartition table
    auto metadata = edgePartition->graph_metadata();
    const std::string& filename = metadata->filename()->str();
    uint64_t numNodes = metadata->num_nodes();
    uint64_t numEdges = metadata->num_edges();
    std::cout << "Graph: " << filename << std::endl;
    std::cout << "Nodes (n): " << numNodes << std::endl;
    std::cout << "Edges (m): " << numEdges << std::endl;

    // Access PartitionConfiguration
    const PartitionInfo::PartitionConfiguration* config = edgePartition->partition_configuration();
    uint32_t k = config->k();
    int seed = config->seed();
    uint64_t streamBuffer = config->stream_buffer();
    int modelMode = config->model_mode();
    int alpha = config->alpha();
    std::cout << "Blocks = k: " << k << std::endl;
    std::cout << "Seed: " << seed << std::endl;
    std::cout << "Stream Buffer: " << streamBuffer << std::endl;
    if(modelMode == -1) std::cout << "Mode: Minimal" << std::endl;
    if(alpha == 0) std::cout << "Alpha: Batch Alpha" << std::endl;

    // Access RunTime
    const PartitionInfo::RunTime* runtime = edgePartition->runtime();
    double ioTime = runtime->io_time();
    double partitionTime = runtime->partition_time();
    double modelConstructionTime = runtime->model_construction_time();
    double mappingTime = runtime->mapping_time();
    double totalTime = runtime->total_time();
    std::cout << "IO time: " << ioTime << std::endl;
    std::cout << "Model Construction time: " << modelConstructionTime << std::endl;
    std::cout << "Mapping time: " << mappingTime << std::endl;
    std::cout << "Partition time: " << partitionTime << std::endl;
    std::cout << "Total time: " << totalTime << std::endl;

    // read partition from input - store in stream_nodes_assign
    if (partition_config.evaluate_mode) {
        if (!partition_config.filename_output.compare("")) {
            partition_config.filename_output = "tmp_output.txt";
        }
    }

    if (partition_config.num_streams_passes >
        1 + partition_config.restream_number) {
        partition_config.stream_total_upperbound = ceil(
                ((100 + 1.5 * partition_config.imbalance) / 100.) *
                (numEdges / (double)partition_config.k));
    } else {
        partition_config.stream_total_upperbound = ceil(
                ((100 + partition_config.imbalance) / 100.) *
                (numEdges / (double)partition_config.k));
    }

    // Access MemoryConsumption
    const PartitionInfo::MemoryConsumption* memory_consumption = edgePartition->memory_consumption();
    long maxRSS = memory_consumption->max_rss();
    std::cout << "Maximum Resident Set Size (KB): " << maxRSS << std::endl;

    // Access PartitionMetrics
    const PartitionInfo::PartitionMetrics* metrics = edgePartition->metrics();
    EdgeWeight edgeCut = metrics->edge_cut();
    NodeID vertexCut = metrics->vertex_cut();
    NodeID replicas = metrics->replicas();
    double balance = metrics->balance();
    double replicationFactor = metrics->replication_factor();

    if(replicationFactor == 0) {
        vertexCut = 0;
        replicas = 0;
        if(partition_config.light_evaluator) {
            graph_io_stream::streamEvaluateEdgePartition(partition_config, graph_filename,
                                                         vertexCut, replicas, replicationFactor, balance);
        } else {
            graph_io_stream::streamEvaluatePartitionEdgeBatch(partition_config, graph_filename, vertexCut,
                                                              replicas, replicationFactor, balance);
        }

        // update metrics
        flatbuffers::FlatBufferBuilder builder(1024);
        auto updatedMetrics = PartitionInfo::CreatePartitionMetrics(builder,edgeCut,
                                                                vertexCut,
                                                                replicas,
                                                                replicationFactor,
                                                                balance);

        // Create a new EdgePartition with the updated metrics
        auto filenameOffset = builder.CreateString(baseFilename);
        auto updatedEdgePartition = PartitionInfo::CreateEdgePartition(
                builder,
                PartitionInfo::CreateGraphMetadata(builder, filenameOffset, edgePartition->graph_metadata()->num_nodes(),
                                               edgePartition->graph_metadata()->num_edges()),
                GraphPartitionInfoInfo::CreatePartitionConfiguration(builder, edgePartition->partition_configuration()->k(),
                                                        edgePartition->partition_configuration()->seed(),
                                                        edgePartition->partition_configuration()->stream_buffer(),
                                                        edgePartition->partition_configuration()->model_mode(),
                                                        edgePartition->partition_configuration()->alpha()),
                PartitionInfo::CreateRunTime(builder, edgePartition->runtime()->io_time(),
                                         edgePartition->runtime()->partition_time(),
                                         edgePartition->runtime()->model_construction_time(),
                                         edgePartition->runtime()->mapping_time(),
                                         edgePartition->runtime()->total_time()),
                PartitionInfo::CreateMemoryConsumption(builder, edgePartition->memory_consumption()->max_rss()),
                updatedMetrics
        );

        // Finish building the FlatBuffer
        builder.Finish(updatedEdgePartition);
        //Step 4: Write to File
        const uint8_t *bufferPointer = builder.GetBufferPointer();
        int bufferSize = builder.GetSize();
        file.close();

        std::string outputFileNameStream_w = baseFilename + "_" + std::to_string(partition_config.k) + "_" + std::to_string(partition_config.stream_buffer_len) + "_stats.bin";
        const char* outputFileName = outputFileNameStream_w.c_str();
        FILE* file_w = fopen(outputFileName, "wb");
        fwrite(bufferPointer, 1, bufferSize, file_w);
        fclose(file_w);
    }

    std::cout << "Vertex Cut: " << vertexCut << std::endl;
    std::cout << "Replicas: " << replicas << std::endl;
    std::cout << "Replication Factor: " << replicationFactor << std::endl;
    std::cout << "Balance: " << balance << std::endl;

    return 0;
}

void config_multibfs_initial_partitioning(PartitionConfig &partition_config) {
    if (partition_config.initial_part_multi_bfs &&
        partition_config.curr_batch >= 2) {
        partition_config.initial_partitioning_type = INITIAL_PARTITIONING_MULTIBFS;
    }
}

// Function to extract the base filename without path and extension
std::string extractBaseFilename(const std::string& fullPath) {
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