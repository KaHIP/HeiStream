/****************************************************************************
 * algorithm_pipeline.cpp
 *****************************************************************************/

#include "core/pipeline/algorithm_pipeline.h"

#include "core/interfaces/algorithm_evaluator.h"
#include "core/interfaces/algorithm_reporter.h"
#include "core/interfaces/algorithm_scheduler.h"
#include "core/timing/scoped_stage_timer.h"
#include "tools/flat_buffer_writer.h"
#include "tools/stream_util.h"


namespace core::pipeline {

void execute_algorithm_pipeline(Config& config, const std::string& graph_filename,
                                StreamContext& ctx, const interfaces::IAlgorithmScheduler& scheduler,
                                const interfaces::IAlgorithmEvaluator& evaluator,
                                const interfaces::IAlgorithmReporter& reporter) {
    FlatBufferWriter fb_writer;
    {
        // Total runtime is owned centrally by the core pipeline.
        TIMED_SCOPE(ctx.timing, core::timing::StreamTimingStage::Total);
        scheduler.execute(config, graph_filename, ctx);
    }
    // RSS snapshot is taken before evaluation.
    ctx.max_rss_kb = getMaxRSS();
    evaluator.evaluate(config, graph_filename, ctx, fb_writer);
    reporter.report(config, graph_filename, ctx, fb_writer);
}

} // namespace core::pipeline
