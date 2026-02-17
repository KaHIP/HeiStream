/******************************************************************************
 * mixed_refinement.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "mixed_refinement.h"

#include "label_propagation_refinement/label_propagation_refinement.h"

mixed_refinement::mixed_refinement() {}

mixed_refinement::~mixed_refinement() {}

EdgeWeight mixed_refinement::perform_refinement(Config& config, graph_access& G,
                                                complete_boundary& boundary) {
    label_propagation_refinement* label_refine = new label_propagation_refinement();

    EdgeWeight overall_improvement = 0;
    // call refinement

    if (!config.initial_partitioning && !config.skip_outer_ls) {
        if (config.stream_input) {
            overall_improvement +=
                (int)label_refine->perform_refinement_fennel(config, G, boundary);
        }
    }

    delete label_refine;

    return overall_improvement;
}
