/******************************************************************************
 * initial_partitioner.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef INITIAL_PARTITIONER_TJKC6RWY
#define INITIAL_PARTITIONER_TJKC6RWY

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"

class initial_partitioner {
   public:
    initial_partitioner();
    virtual ~initial_partitioner();

    virtual void initial_partition(Config& config, const unsigned int seed, graph_access& G,
                                   int* xadj, int* adjncy, int* vwgt, int* adjwgt,
                                   int* partition_map, int ismultisec = 0,
                                   int stream_pass_index = 0) = 0;

    virtual void initial_partition(Config& config, const unsigned int seed, graph_access& G,
                                   int* partition_map, int ismultisec = 0,
                                   int stream_pass_index = 0) = 0;
};

#endif /* end of include guard: INITIAL_PARTITIONER_TJKC6RWY */
