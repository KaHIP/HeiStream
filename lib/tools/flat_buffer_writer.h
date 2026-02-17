/******************************************************************************
 * flat_buffer_writer.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <adilchhabra7@gmail.com>
 *****************************************************************************/

#ifndef KAHIP_FLATBUFFERWRITER_H
#define KAHIP_FLATBUFFERWRITER_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "PartitionInfo_generated.h"
#include "partition/partition_config.h"
#include "tools/stream_util.h"

class FlatBufferWriter {
   private:
    double io_time_;
    double model_construction_time_;
    double postprocess_time_;
    double buffer_maintenance_time_;
    double partition_time_;
    double total_time_;
    long maxRSS_;
    EdgeWeight edge_cut_;
    NodeID vertex_cut_;
    NodeID replicas_;
    double replication_factor_;
    double balance_;

   public:
    FlatBufferWriter()
        : io_time_(0.0),
          model_construction_time_(0.0),
          postprocess_time_(0.0),
          buffer_maintenance_time_(0.0),
          partition_time_(0.0),
          total_time_(0.0),
          maxRSS_(0),
          edge_cut_(-1),
          vertex_cut_(-1),
          replicas_(-1),
          replication_factor_(0.0),
          balance_(0.0) {}

    void updateResourceConsumption(double& io_time, double& model_construction_time,
                                   double& postprocess_time, double& buffer_maintenance_time,
                                   double& partition_time, double& total_time, long& maxRSS) {
        io_time_ = io_time;
        model_construction_time_ = model_construction_time;
        postprocess_time_ = postprocess_time;
        buffer_maintenance_time_ = buffer_maintenance_time;
        partition_time_ = partition_time;
        total_time_ = total_time;
        maxRSS_ = maxRSS;
    }

    void updateVertexPartitionResults(EdgeWeight& edge_cut, double balance) {
        edge_cut_ = edge_cut;
        balance_ = balance;
    }

    void updateEdgePartitionResults(NodeID& vertex_cut, NodeID& replicas, double replication_factor,
                                    double balance) {
        vertex_cut_ = vertex_cut;
        replicas_ = replicas;
        replication_factor_ = replication_factor;
        balance_ = balance;
    }

    void print_summary(const std::string& graph_filename, const Config& partition_config) const {
        (void)graph_filename;
        const int key_width = 28;
        const int val_width = 22;

        auto print_line = [&](const std::string& key, const std::string& value) {
            std::cout << "| " << std::left << std::setw(key_width) << key << " | " << std::right
                      << std::setw(val_width) << value << " |\n";
        };

        auto print_sep = [&]() {
            std::cout << "+" << std::string(key_width + 2, '-') << "+"
                      << std::string(val_width + 2, '-') << "+\n";
        };

        std::ostringstream ratio_ss;
        if (partition_config.total_edges > 0 && edge_cut_ >= 0) {
            ratio_ss << std::fixed << std::setprecision(6)
                     << (edge_cut_ / static_cast<double>(partition_config.total_edges));
        } else {
            ratio_ss << "n/a";
        }

        std::cout << "\n================ Stream Run Summary ================\n";
        print_sep();
        print_line("Output Filename", partition_config.filename_output);
        print_sep();
#ifdef ENABLE_TIME_MEASUREMENTS
        print_line("IO time (s)", to_fixed(io_time_));
        print_line("Model time (s)", to_fixed(model_construction_time_));
        print_line("Postprocess time (s)", to_fixed(postprocess_time_));
        print_line("Buffer maintenance (s)", to_fixed(buffer_maintenance_time_));
        print_line("Partition time (s)", to_fixed(partition_time_));
#else
        print_line("Scope Timing", "Profiling Disabled");
#endif
        print_line("Total time (s)", to_fixed(total_time_));
        print_line("Max RSS (KB)", std::to_string(maxRSS_));
        print_sep();
        if (!partition_config.edge_partition) {
            print_line("Edge cut", std::to_string(edge_cut_));
            print_line("Cut ratio", ratio_ss.str());
        } else {
            print_line("Vertex cut", std::to_string(vertex_cut_));
            print_line("Replicas", std::to_string(replicas_));
            print_line("Replication factor", to_fixed(replication_factor_));
        }
        print_line("Balance", to_fixed(balance_));
        print_sep();
    }

    void write(const std::string& graph_filename, const Config& partition_config) const {
        flatbuffers::FlatBufferBuilder builder(1024);
        std::stringstream filename;
        std::string baseFilename = ::extractBaseFilename(graph_filename);
        PartitionInfo::GraphMetadataBuilder metadata_builder(builder);
        auto filenameOffset = builder.CreateString(baseFilename);
        metadata_builder.add_filename(filenameOffset);
        metadata_builder.add_num_nodes(partition_config.total_nodes);
        metadata_builder.add_num_edges(partition_config.total_edges);
        auto metadata = metadata_builder.Finish();
        builder.Finish(metadata);

        PartitionInfo::PartitionConfigurationBuilder config_builder(builder);
        uint64_t stream_buffer_value = partition_config.batch_size;
        config_builder.add_k(partition_config.k);
        config_builder.add_seed(partition_config.seed);
        config_builder.add_stream_buffer(stream_buffer_value);
        if (partition_config.minimal_mode)
            config_builder.add_model_mode(-1);
        if (partition_config.batch_alpha)
            config_builder.add_alpha(0);
        config_builder.add_batch_size(partition_config.batch_size);
        config_builder.add_max_buffer_size(partition_config.max_buffer_size);
        config_builder.add_buffer_score_type(static_cast<int>(partition_config.buffer_score_type));
        config_builder.add_bq_disc_factor(partition_config.bq_disc_factor);
        config_builder.add_haa_theta(partition_config.haa_theta);
        config_builder.add_cbs_theta(partition_config.cbs_theta);
        config_builder.add_bb_ratio(partition_config.bb_ratio);
        auto configdata = config_builder.Finish();
        builder.Finish(configdata);

        PartitionInfo::RunTimeBuilder runtime_builder(builder);
        runtime_builder.add_io_time(io_time_);
        runtime_builder.add_model_construction_time(model_construction_time_);
        runtime_builder.add_postprocess_time(postprocess_time_);
        runtime_builder.add_buffer_maintenance_time(buffer_maintenance_time_);
        runtime_builder.add_partition_time(partition_time_);
        runtime_builder.add_total_time(total_time_);
        auto runtimedata = runtime_builder.Finish();
        builder.Finish(runtimedata);

        // Create PartitionMetrics
        PartitionInfo::PartitionMetricsBuilder partition_metrics_builder(builder);
        if (!partition_config.edge_partition) {
            partition_metrics_builder.add_edge_cut(edge_cut_);
        } else {
            partition_metrics_builder.add_vertex_cut(vertex_cut_);
            partition_metrics_builder.add_replicas(replicas_);
            partition_metrics_builder.add_replication_factor(replication_factor_);
        }
        partition_metrics_builder.add_balance(balance_);
        auto partition_metrics = partition_metrics_builder.Finish();
        builder.Finish(partition_metrics);

        // Create MemoryConsumption
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

        // Step 4: Write to File
        const uint8_t* bufferPointer = builder.GetBufferPointer();
        int bufferSize = builder.GetSize();

        std::string outputFileNameStream;
        if (!partition_config.edge_partition && partition_config.max_buffer_size > 1) {
            outputFileNameStream = baseFilename + "_" + std::to_string(partition_config.k) + "_" +
                                   std::to_string(partition_config.batch_size) + "_" +
                                   std::to_string(partition_config.max_buffer_size) + ".bin";
        } else {
            outputFileNameStream = baseFilename + "_" + std::to_string(partition_config.k) + "_" +
                                   std::to_string(partition_config.batch_size) + ".bin";
        }
        const char* outputFileName = outputFileNameStream.c_str();
        FILE* file = fopen(outputFileName, "wb");
        fwrite(bufferPointer, 1, bufferSize, file);
        fclose(file);
    }

   private:
    static std::string to_fixed(double value) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(6) << value;
        return os.str();
    }
};

#endif  // KAHIP_FLATBUFFERWRITER_H
