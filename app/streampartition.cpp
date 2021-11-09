/******************************************************************************
 * streampartition.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj <marcelofaraj@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 
#include <vector>
#include <fstream>
#include <sstream>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "graph_io_stream.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "partition/uncoarsening/refinement/mixed_refinement.h"
#include "partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.h"

#define MIN(A,B) ((A)<(B))?(A):(B)
#define MAX(A,B) ((A)>(B))?(A):(B)


void config_multibfs_initial_partitioning(PartitionConfig & partition_config);


int main(int argn, char **argv) {
        PartitionConfig partition_config;
        std::string graph_filename;
	LINE_BUFFER lines = NULL;
        timer t, processing_t, io_t;
	EdgeWeight total_edge_cut = 0;
        double global_mapping_time = 0;
	double buffer_mapping_time = 0;
	double buffer_io_time = 0;
        quality_metrics qm;
	balance_configuration bc;

        bool is_graph_weighted = false;
        bool suppress_output   = false;
        bool recursive         = false;

        int ret_code = parse_parameters(argn, argv, 
                                        partition_config, 
                                        graph_filename, 
                                        is_graph_weighted, 
                                        suppress_output, recursive); 

        if(ret_code) {
                return 0;
        }

        std::ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
                std::cout.rdbuf(ofs.rdbuf()); 
        }
	srand(partition_config.seed);
	random_functions::setSeed(partition_config.seed);

        partition_config.LogDump(stdout);
	partition_config.stream_input = true;
	graph_access *G = new graph_access(); 


	int &passes = partition_config.num_streams_passes;
	for (partition_config.restream_number=0; partition_config.restream_number<passes; partition_config.restream_number++) {

		io_t.restart();
		graph_io_stream::readFirstLineStream(partition_config, graph_filename, total_edge_cut);
		graph_io_stream::loadRemainingLines(partition_config, lines);
		buffer_io_time += io_t.elapsed();

		while (partition_config.remaining_stream_nodes) {
			io_t.restart();
			int streamed_nodes = graph_io_stream::readStreamBuffer(partition_config, *G, lines);
			G->set_partition_count(partition_config.k); 
			buffer_io_time += io_t.elapsed();

			if (streamed_nodes == 0) {
				break;
			}

			buffer_mapping_time = 0;
			graph_io_stream::countAssignedNodes(partition_config);
			graph_io_stream::prescribeBufferInbalance(partition_config);

			bool already_fully_partitioned = (partition_config.restream_vcycle && partition_config.restream_number) ;
			bc.configurate_balance( partition_config, *G, already_fully_partitioned || !partition_config.stream_initial_bisections);
			config_multibfs_initial_partitioning(partition_config);


			// ***************************** perform partitioning ***************************************       
			t.restart();
			graph_partitioner partitioner;



			if(partition_config.time_limit == 0) {
				partitioner.perform_partitioning(partition_config, *G);
			} else {
				PartitionID* map = new PartitionID[G->number_of_nodes()];
				EdgeWeight best_cut = std::numeric_limits<EdgeWeight>::max();
				while(t.elapsed() < partition_config.time_limit) {
					partition_config.graph_allready_partitioned = false;
					partitioner.perform_partitioning(partition_config, *G);
					EdgeWeight cut = qm.edge_cut(*G);
					if(cut < best_cut) {
						best_cut = cut;
						forall_nodes((*G), node) {
							map[node] = G->getPartitionIndex(node);
						} endfor
					}
				}

				forall_nodes((*G), node) {
					G->setPartitionIndex(node, map[node]);
				} endfor
				delete[] map;
			}


			ofs.close();
			
			buffer_mapping_time += t.elapsed();

			io_t.restart();
			graph_io_stream::generalizeStreamPartition(partition_config, *G);
			buffer_io_time += io_t.elapsed();
			global_mapping_time += buffer_mapping_time;

			 
			// write batch partition to the disc 
			if (partition_config.stream_output_progress) {
				std::stringstream filename;
				if(!partition_config.filename_output.compare("")) {
					filename << "tmppartition" << partition_config.k << "-" << partition_config.curr_batch;
				} else {
					filename << partition_config.filename_output << "-" << partition_config.curr_batch;
				}
				if(partition_config.restream_number) {
					filename << ".R" << partition_config.restream_number;
				}
				graph_io_stream::writePartitionStream(partition_config, filename.str());
			}
	 
		}

		if (partition_config.ram_stream) {
			delete lines;
		}
        }

	double total_time = processing_t.elapsed();

        std::cout << "Total processing time: " << total_time  << std::endl;
        std::cout << "io time: " << buffer_io_time  << std::endl;
        std::cout << "partitioning/mapping time in total: " << global_mapping_time  << std::endl;

	graph_io_stream::streamEvaluatePartition(partition_config, graph_filename, total_edge_cut);

        std::cout << "cut \t\t"         << total_edge_cut                        << std::endl;
        std::cout << "finalobjective  " << total_edge_cut			 << std::endl;
        std::cout << "balance \t"  << qm.balance_full_stream(*partition_config.stream_blocks_weight) << std::endl;

        // write the partition to the disc 
        std::stringstream filename;
        if(!partition_config.filename_output.compare("")) {
                filename << "tmppartition" << partition_config.k;
        } else {
                filename << partition_config.filename_output;
        }

        if (!partition_config.suppress_output) {
                graph_io_stream::writePartitionStream(partition_config, filename.str());
        } else {
                std::cout << "No partition will be written as output." << std::endl;
        }

	delete G;
	
	if (partition_config.ghostkey_to_edges != NULL) {
		delete partition_config.ghostkey_to_edges;
	}

	return 0;
}



void config_multibfs_initial_partitioning(PartitionConfig & partition_config) {
	if (partition_config.initial_part_multi_bfs && partition_config.curr_batch >= 2) {
		partition_config.initial_partitioning_type = INITIAL_PARTITIONING_MULTIBFS;
	}
}


