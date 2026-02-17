/****************************************************************************
 * algorithm_registry.cpp
 *****************************************************************************/

#include "core/factory/algorithm_registry.h"

#include <unordered_map>
#include <utility>

#include "core/interfaces/stream_algorithm.h"

namespace {

struct StreamModeHash {
    std::size_t operator()(StreamMode mode) const {
        return static_cast<std::size_t>(mode);
    }
};

std::unordered_map<StreamMode, core::factory::AlgorithmCreator, StreamModeHash>& registry() {
    static std::unordered_map<StreamMode, core::factory::AlgorithmCreator, StreamModeHash> map;
    return map;
}

}  // namespace


namespace core::factory {

void register_algorithm(StreamMode mode, AlgorithmCreator creator) {
    registry()[mode] = std::move(creator);
}

std::unique_ptr<IStreamAlgorithm> create_registered_algorithm(StreamMode mode) {
    auto it = registry().find(mode);
    if (it == registry().end()) {
        return {};
    }
    return it->second();
}

} // namespace core::factory

