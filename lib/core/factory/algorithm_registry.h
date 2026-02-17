/****************************************************************************
 * algorithm_registry.h
 *****************************************************************************/

#ifndef CORE_ALGORITHM_REGISTRY_H_
#define CORE_ALGORITHM_REGISTRY_H_

#include <functional>
#include <memory>

#include "definitions.h"

class IStreamAlgorithm;


namespace core::factory {

using AlgorithmCreator = std::function<std::unique_ptr<IStreamAlgorithm>()>;

// Register one concrete algorithm implementation for a stream mode.
void register_algorithm(StreamMode mode, AlgorithmCreator creator);
// Create an algorithm instance for the requested mode; returns null if unregistered.
std::unique_ptr<IStreamAlgorithm> create_registered_algorithm(StreamMode mode);

} // namespace core::factory


#endif /* CORE_ALGORITHM_REGISTRY_H_ */
