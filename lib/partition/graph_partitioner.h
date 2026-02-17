/******************************************************************************
 * graph_partitioner.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARTITION_OL9XTLU4
#define PARTITION_OL9XTLU4

#include "coarsening/coarsening.h"
#include "coarsening/stop_rules/stop_rules.h"
#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "uncoarsening/refinement/refinement.h"

class graph_partitioner {
   public:
    graph_partitioner();
    virtual ~graph_partitioner();

    void perform_partitioning(Config& graph_partitioner_config, graph_access& G,
                              int stream_pass_index = 0);
    void perform_recursive_partitioning(Config& graph_partitioner_config, graph_access& G,
                                        int ismultisec = 0);

   private:
    void perform_recursive_partitioning_internal(Config& graph_partitioner_config, graph_access& G,
                                                 PartitionID lb, PartitionID ub,
                                                 int ismultisec = 0);

    void single_run(Config& config, graph_access& G, int stream_pass_index = 0);

    unsigned m_global_k;
    int m_global_upper_bound;
    int m_rnd_bal;
};

#endif /* end of include guard: PARTITION_OL9XTLU4 */
