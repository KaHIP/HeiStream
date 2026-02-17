/****************************************************************************
 * register_builtin_algorithms.cpp
 *****************************************************************************/

#include "algorithms/bootstrap/register_builtin_algorithms.h"

#include <mutex>

#include "algorithms/edge_partitioning/factory/edge_partitioning_algorithm_factory.h"
#include "algorithms/node_partitioning/factory/node_partitioning_algorithm_factory.h"
#include "core/factory/algorithm_registry.h"
#include "core/interfaces/stream_algorithm.h"


namespace algorithms::bootstrap {

void register_builtin_algorithms() {
    static std::once_flag once;
    std::call_once(once, []() {
        // Built-in stream algorithm registrations.
        // Add new algorithms here after introducing a new StreamMode value.
        core::factory::register_algorithm(
            StreamMode::Node, []() { return node_partitioning::factory::create_algorithm(); });
        core::factory::register_algorithm(
            StreamMode::Edge, []() { return edge_partitioning::factory::create_algorithm(); });
    });
}

} // namespace algorithms::bootstrap
