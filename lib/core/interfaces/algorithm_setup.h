/****************************************************************************
 * algorithm_setup.h
 *****************************************************************************/

#ifndef CORE_ALGORITHM_SETUP_H_
#define CORE_ALGORITHM_SETUP_H_

#include <string>

#include "partition/partition_config.h"


namespace core::interfaces {

/**
 * Algorithm-level setup hook before execution begins.
 *
 * @details Implementations should perform run-scoped initialization that is
 *          independent of individual stream passes.
 */
class IAlgorithmSetup {
   public:
    virtual ~IAlgorithmSetup() {}

    /**
     * Prepare algorithm state and configuration before scheduling.
     *
     * @param config Mutable run configuration.
     * @param graph_filename Input graph path.
     */
    virtual void prepare_run(Config& config, const std::string& graph_filename) const = 0;
};

} // namespace core::interfaces


#endif /* CORE_ALGORITHM_SETUP_H_ */
