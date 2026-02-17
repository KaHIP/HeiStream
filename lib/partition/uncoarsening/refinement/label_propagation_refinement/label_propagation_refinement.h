/******************************************************************************
 * label_propagation_refinement.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef LABEL_PROPAGATION_REFINEMENT_R4XW141Y
#define LABEL_PROPAGATION_REFINEMENT_R4XW141Y

#include "../refinement.h"
#include "definitions.h"

class label_propagation_refinement : public refinement {
   public:
    label_propagation_refinement();
    virtual ~label_propagation_refinement();

    virtual EdgeWeight perform_refinement(Config& partition_config, graph_access& G,
                                          complete_boundary& boundary);

    virtual EdgeWeight perform_refinement_fennel(Config& config, graph_access& G,
                                                 complete_boundary& boundary);

    bool is_boundary(NodeID node, graph_access& G);
};

#endif /* end of include guard: LABEL_PROPAGATION_REFINEMENT_R4XW141Y */
