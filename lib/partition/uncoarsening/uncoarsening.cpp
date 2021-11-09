/******************************************************************************
 * uncoarsening.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "graph_partition_assertions.h"
#include "misc.h"
#include "quality_metrics.h"
#include "refinement/mixed_refinement.h"
#include "refinement/label_propagation_refinement/label_propagation_refinement.h"
#include "refinement/refinement.h"
#include "tools/random_functions.h"
#include "uncoarsening.h"


uncoarsening::uncoarsening() {

}

uncoarsening::~uncoarsening() {

}

int uncoarsening::perform_uncoarsening(PartitionConfig & config, graph_hierarchy & hierarchy) {

	return perform_uncoarsening_cut(config, hierarchy);
}

int uncoarsening::perform_uncoarsening_cut(PartitionConfig & config, graph_hierarchy & hierarchy) {
        int improvement = 0;

        PartitionConfig cfg     = config;
        refinement* refine      = NULL;

	refine      = new mixed_refinement();

        graph_access * coarsest = hierarchy.get_coarsest();
        PRINT(std::cout << "log>" << "unrolling graph with " << coarsest->number_of_nodes() << std::endl;)

        complete_boundary* finer_boundary   = NULL;
        complete_boundary* coarser_boundary = NULL;
        double factor = config.balance_factor;
        cfg.upper_bound_partition = ((!hierarchy.isEmpty()) * factor +1.0)*config.upper_bound_partition;

	improvement += (int)refine->perform_refinement(cfg, *coarsest, *coarser_boundary);

        graph_access* to_delete = NULL;
        unsigned int hierarchy_deepth = hierarchy.size();

        while(!hierarchy.isEmpty()) {
                bool use_delta_gains = false;
                if (!cfg.initial_partitioning) {
                        use_delta_gains = cfg.use_delta_gains;
                }
                graph_access* G = hierarchy.pop_finer_and_project(config, use_delta_gains);

                PRINT(std::cout << "log>" << "unrolling graph with " << G->number_of_nodes()<<  std::endl;)
                
                //call refinement
                double cur_factor = factor/(hierarchy_deepth-hierarchy.size());
                cfg.upper_bound_partition = ((!hierarchy.isEmpty()) * cur_factor+1.0)*config.upper_bound_partition;
                PRINT(std::cout <<  "cfg upperbound " <<  cfg.upper_bound_partition  << std::endl;)

		improvement += (int)refine->perform_refinement(cfg, *G, *finer_boundary);

                ASSERT_TRUE(graph_partition_assertions::assert_graph_has_kway_partition(config, *G));

                if(config.use_balance_singletons && !config.label_propagation_refinement) {
                        if (config.adapt_bal && config.initial_partitioning) {
                                finer_boundary->balance_singletons_het_bal( config, *G );
                        } else {
                                finer_boundary->balance_singletons( config, *G );
                        }
                }

                // update boundary pointers
                coarser_boundary = finer_boundary;

		//clean up 
		if(to_delete != NULL) {
			delete to_delete;
		}
		if(!hierarchy.isEmpty()) {
			to_delete = G;
		}

        }

	 
        delete refine;
        if(finer_boundary != NULL) delete finer_boundary;
	delete coarsest;

        return improvement;
}


