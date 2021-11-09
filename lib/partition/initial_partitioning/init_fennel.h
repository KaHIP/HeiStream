/******************************************************************************
 * init_fennel.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef INITFENNEL_7I4IR31Y
#define INITFENNEL_7I4IR31Y

#include "initial_partitioner.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "uncoarsening/refinement/refinement.h"
#include "definitions.h"
#include <algorithm>

class init_fennel : public initial_partitioner {
        public:
                init_fennel();
                virtual ~init_fennel();
                
                void initial_partition( PartitionConfig & config, const unsigned int seed,  
                                graph_access & G, int* partition_map, int ismultisec=0); 

                void initial_partition( PartitionConfig & config, const unsigned int seed,  
                                graph_access & G, 
                                int* xadj,
                                int* adjncy, 
                                int* vwgt, 
                                int* adjwgt,
                                int* partition_map, 
                                int ismultisec=0); 

        private:
		EdgeWeight fennel(PartitionConfig & partition_config, graph_access & G);

};


#endif /* end of include guard: BIPARTITION_7I4IR31Y */
