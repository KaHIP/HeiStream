/****************************************************************************
 * stream_engine.cpp
 *****************************************************************************/

#include "core/pipeline/stream_engine.h"

#include "core/context/stream_context.h"
#include "core/factory/algorithm_registry.h"
#include "core/interfaces/stream_algorithm.h"
#include "core/setup/stream_run_setup.h"


namespace core::pipeline {

void run_stream_engine(StreamMode mode, Config& config, const std::string& graph_filename) {
    if (!setup::prepare_framework_run(config, graph_filename)) {
        return;
    }

    StreamContext ctx(config, mode);
    std::unique_ptr<IStreamAlgorithm> algorithm = factory::create_registered_algorithm(mode);
    if (!algorithm) {
        return;
    }
    algorithm->prepare_run(config, graph_filename, ctx);
    algorithm->run(config, graph_filename, ctx);
    algorithm->finalize_run(config, graph_filename, ctx);
}

} // namespace core::pipeline
