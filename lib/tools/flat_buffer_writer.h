/******************************************************************************
 * flat_buffer_writer.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <adilchhabra7@gmail.com>
 *****************************************************************************/

#ifndef KAHIP_FLATBUFFERWRITER_H
#define KAHIP_FLATBUFFERWRITER_H

#include <fstream>
#include <iostream>
#include <vector>
#include "PartitionInfo_generated.h"
#include "partition/partition_config.h"
#include <iomanip>

class FlatBufferWriter {
private:
    double io_time_;
    double model_construction_time_;
    double mapping_time_;
    double partition_time_;
    double total_time_;
    long maxRSS_;
    EdgeWeight edge_cut_;
    NodeID vertex_cut_;
    NodeID replicas_;
    double replication_factor_ ;
    double balance_;

public:
    FlatBufferWriter()
            : io_time_(0.0), model_construction_time_(0.0), mapping_time_(0.0),
              partition_time_(0.0), total_time_(0.0), maxRSS_(0), edge_cut_(-1),
              vertex_cut_(-1), replicas_(-1), replication_factor_(0.0), balance_ (0.0){}

    void updateResourceConsumption(double &io_time, double &model_construction_time,
                                    double &mapping_time, double &partition_time, double &total_time, long &maxRSS) {
        io_time_ = io_time;
        model_construction_time_ = model_construction_time;
        mapping_time_ = mapping_time;
        partition_time_ = partition_time;
        total_time_ = total_time;
        maxRSS_ = maxRSS;
    }

    void updateVertexPartitionResults(EdgeWeight & edge_cut, double balance) {
        edge_cut_ = edge_cut;
        balance_ = balance;
    }

    void updateEdgePartitionResults(NodeID & vertex_cut, NodeID & replicas, double replication_factor, double balance) {
        vertex_cut_ = vertex_cut;
        replicas_ = replicas;
        replication_factor_ = replication_factor;
        balance_ = balance;
    }

    // Function to extract the base filename without path and extension
    static std::string extractBaseFilename(const std::string& fullPath) {
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

    void write(const std::string& graph_filename, const PartitionConfig &partition_config ) const {
        flatbuffers::FlatBufferBuilder builder(1024);
        std::stringstream filename;
        std::string baseFilename = extractBaseFilename(graph_filename);
        PartitionInfo::GraphMetadataBuilder metadata_builder(builder);
        auto filenameOffset = builder.CreateString(baseFilename);
        metadata_builder.add_filename(filenameOffset);
        metadata_builder.add_num_nodes(partition_config.total_nodes);
        metadata_builder.add_num_edges(partition_config.total_edges);
        auto metadata = metadata_builder.Finish();
        builder.Finish(metadata);
        std::cout << "Graph: " << baseFilename << std::endl;
        std::cout << "Nodes (n): " << partition_config.total_nodes << std::endl;
        std::cout << "Edges (m): " << partition_config.total_edges << std::endl;

        PartitionInfo::PartitionConfigurationBuilder config_builder(builder);
        config_builder.add_k(partition_config.k);
        config_builder.add_seed(partition_config.seed);
        config_builder.add_stream_buffer(partition_config.stream_buffer_len);
        if(partition_config.minimal_mode) config_builder.add_model_mode(-1);
        if(partition_config.batch_alpha) config_builder.add_alpha(0);
        auto configdata = config_builder.Finish();
        builder.Finish(configdata);
        std::cout << "Blocks = k: " << partition_config.k << std::endl;
        std::cout << "Seed: " << partition_config.seed << std::endl;
        std::cout << "Stream Buffer: " << partition_config.stream_buffer_len << std::endl;
        if(partition_config.edge_partition) {
            if (partition_config.minimal_mode) std::cout << "Mode: Minimal" << std::endl;
            if (partition_config.batch_alpha) std::cout << "Alpha: Batch Alpha" << std::endl;
        }

        PartitionInfo::RunTimeBuilder runtime_builder(builder);
        runtime_builder.add_io_time(io_time_);
        std::cout << "IO time: " << io_time_ << std::endl;
        runtime_builder.add_model_construction_time(model_construction_time_);
        std::cout << "Model Construction time: " << model_construction_time_ << std::endl;
        runtime_builder.add_mapping_time(mapping_time_);
        std::cout << "Mapping time: " << mapping_time_ << std::endl;
        runtime_builder.add_partition_time(partition_time_);
        std::cout << "Partition time: " << partition_time_ << std::endl;
        runtime_builder.add_total_time(total_time_);
        std::cout << "Total time: " << total_time_ << std::endl;
        auto runtimedata = runtime_builder.Finish();
        builder.Finish(runtimedata);

        // Create PartitionMetrics
        PartitionInfo::PartitionMetricsBuilder partition_metrics_builder(builder);
        if(!partition_config.edge_partition) {
            partition_metrics_builder.add_edge_cut(edge_cut_);
            std::cout << "Edge Cut: " << edge_cut_ << std::endl;
        } else {
            partition_metrics_builder.add_vertex_cut(vertex_cut_);
            partition_metrics_builder.add_replicas(replicas_);
            partition_metrics_builder.add_replication_factor(replication_factor_);
            if(replication_factor_ != 0){
                std::cout << "Vertex Cut: " << vertex_cut_ << std::endl;
                std::cout << "Replicas: " << replicas_ << std::endl;
                std::cout << "Replication Factor: " << replication_factor_ << std::endl;
            }
        }
        partition_metrics_builder.add_balance(balance_);
        std::cout << "Balance: " << balance_ << std::endl;
        auto partition_metrics = partition_metrics_builder.Finish();
        builder.Finish(partition_metrics);

        // Create MemoryConsumption
        if (maxRSS_ != -1) {
            std::cout << "Maximum Resident Set Size (KB): " << maxRSS_ << std::endl;
        }
        auto memory_consumption = PartitionInfo::CreateMemoryConsumption(builder, maxRSS_);

        // Create PartitionLog
        PartitionInfo::PartitionLogBuilder partition_log_builder(builder);
        partition_log_builder.add_graph_metadata(metadata);
        partition_log_builder.add_partition_configuration(configdata);
        partition_log_builder.add_runtime(runtimedata);
        partition_log_builder.add_memory_consumption(memory_consumption);
        partition_log_builder.add_metrics(partition_metrics);
        auto partition_log = partition_log_builder.Finish();
        builder.Finish(partition_log);

        //Step 4: Write to File
        const uint8_t* bufferPointer = builder.GetBufferPointer();
        int bufferSize = builder.GetSize();

        std::string outputFileNameStream;
        outputFileNameStream = partition_config.output_path + baseFilename + "_" + std::to_string(partition_config.k) + "_" + std::to_string(partition_config.stream_buffer_len) + ".bin";
        const char* outputFileName = outputFileNameStream.c_str();
        if(partition_config.write_results) {
            FILE *file = fopen(outputFileName, "wb");
            fwrite(bufferPointer, 1, bufferSize, file);
            fclose(file);
        }
    }
};


#endif //KAHIP_FLATBUFFERWRITER_H