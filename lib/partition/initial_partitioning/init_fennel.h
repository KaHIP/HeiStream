/******************************************************************************
 * init_fennel.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef INITFENNEL_7I4IR31Y
#define INITFENNEL_7I4IR31Y

#include <algorithm>

#include "definitions.h"
#include "initial_partitioner.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "uncoarsening/refinement/refinement.h"

class init_fennel : public initial_partitioner {
   public:
    init_fennel();
    virtual ~init_fennel();

    void initial_partition(Config& config, const unsigned int seed, graph_access& G,
                           int* partition_map, int ismultisec = 0,
                           int stream_pass_index = 0) override;

    void initial_partition(Config& config, const unsigned int seed, graph_access& G, int* xadj,
                           int* adjncy, int* vwgt, int* adjwgt, int* partition_map,
                           int ismultisec = 0, int stream_pass_index = 0) override;

   private:
    EdgeWeight fennel(Config& partition_config, graph_access& G, bool is_restream_pass);
};

#endif /* end of include guard: BIPARTITION_7I4IR31Y */
