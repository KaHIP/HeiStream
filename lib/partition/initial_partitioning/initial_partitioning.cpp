/******************************************************************************
 * initial_partitioning.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "initial_partitioning.h"

#include <math.h>

#include "core/timing/scoped_stage_timer.h"
#include "graph_partition_assertions.h"
#include "graph_partitioner.h"
#include "init_fennel.h"
#include "initial_refinement/initial_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"

initial_partitioning::initial_partitioning() {}

initial_partitioning::~initial_partitioning() {}

void initial_partitioning::perform_initial_partitioning(Config& config, graph_hierarchy& hierarchy,
                                                        int stream_pass_index) {
    graph_access& G = *hierarchy.get_coarsest();
    perform_initial_partitioning(config, G, stream_pass_index);
}

void initial_partitioning::perform_initial_partitioning(Config& config, graph_access& G,
                                                        int stream_pass_index) {
    graph_partitioner partitioner;
    initial_partitioner* partition = nullptr;
    switch (config.initial_partitioning_type) {
        case INITIAL_PARTITIONING_FENNEL:
            partition = new init_fennel();
            break;
    }

    quality_metrics qm;
    EdgeWeight best_cut;
    int* best_map = new int[G.number_of_nodes()];
    if (config.graph_allready_partitioned && !config.omit_given_partitioning) {
        best_cut = qm.edge_cut(G);
        forall_nodes(G, n) {
            best_map[n] = G.getPartitionIndex(n);
        }
        endfor
    } else {
        best_cut = std::numeric_limits<EdgeWeight>::max();
    }

    double t = 0.0;
    TIMED_SCOPE(t, "initial_partitioning_total");
    int* partition_map = new int[G.number_of_nodes()];
    unsigned reps_to_do = (unsigned)std::max(
        (int)ceil(config.initial_partitioning_repetitions / (double)log2(config.k)), 2);

    if (config.initial_partitioning_repetitions == 0) {
        reps_to_do = 1;
    }
    if (config.eco) {
        // bound the number of initial partitioning repetions
        reps_to_do = std::min((int)config.minipreps, (int)reps_to_do);
    }

    PRINT(std::cout << "no of initial partitioning repetitions = " << reps_to_do << '\n';);
    PRINT(std::cout << "no of nodes for partition = " << G.number_of_nodes() << '\n';);
    if (!((config.graph_allready_partitioned && config.no_new_initial_partitioning) ||
          config.omit_given_partitioning)) {
        for (unsigned int rep = 0; rep < reps_to_do; rep++) {
            unsigned seed = random_functions::nextInt(0, std::numeric_limits<int>::max());
            Config working_config = config;
            working_config.combine = false;

            partition->initial_partition(working_config, seed, G, partition_map, 0,
                                         stream_pass_index);

            EdgeWeight cur_cut;

            if (config.use_fennel_objective) {
                cur_cut = qm.fennel_objective(config, G, partition_map, config.fennel_gamma,
                                              config.fennel_alpha);
            } else {
                cur_cut = qm.edge_cut(G, partition_map);
            }

            if (cur_cut < best_cut) {
                PRINT(std::cout << "log>" << "improved the current initial partitiong from "
                                << best_cut << " to " << cur_cut << '\n';)

                forall_nodes(G, n) {
                    best_map[n] = partition_map[n];
                }
                endfor

                    best_cut = cur_cut;
                if (best_cut == 0)
                    break;
            }
        }

        forall_nodes(G, n) {
            G.setPartitionIndex(n, best_map[n]);
        }
        endfor
    }

    G.set_partition_count(config.k);

    PRINT(std::cout << "initial partitioning took " << t << '\n';)
    PRINT(std::cout << "log>" << "current initial balance " << qm.balance(G) << '\n';)

    if (config.initial_partition_optimize || config.combine) {
        initial_refinement iniref;
        iniref.optimize(config, G, best_cut);
    }

    PRINT(std::cout << "log>" << "final current initial partitiong from " << best_cut << " to "
                    << best_cut << '\n';)

    if (!(config.graph_allready_partitioned && config.no_new_initial_partitioning)) {
        PRINT(std::cout << "finalinitialcut " << best_cut << '\n';)
        PRINT(std::cout << "log>" << "final current initial balance " << qm.balance(G) << '\n';)
    }

    ASSERT_TRUE(graph_partition_assertions::assert_graph_has_kway_partition(config, G));

    delete[] partition_map;
    delete[] best_map;
    delete partition;
}
