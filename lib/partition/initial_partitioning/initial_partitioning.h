/******************************************************************************
 * initial_partitioning.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef INITIAL_PARTITIONING_D7VA0XO9
#define INITIAL_PARTITIONING_D7VA0XO9

#include "data_structure/graph_hierarchy.h"
#include "partition/partition_config.h"

class initial_partitioning {
   public:
    initial_partitioning();
    virtual ~initial_partitioning();
    void perform_initial_partitioning(Config& config, graph_hierarchy& hierarchy,
                                      int stream_pass_index = 0);
    void perform_initial_partitioning(Config& config, graph_access& G, int stream_pass_index = 0);
};

#endif /* end of include guard: INITIAL_PARTITIONING_D7VA0XO9 */
