/****************************************************************************
 * scoped_stage_timer.h
 *****************************************************************************/

#ifndef CORE_SCOPED_STAGE_TIMER_H_
#define CORE_SCOPED_STAGE_TIMER_H_

#include <string>
#include <utility>

#include "core/timing/stream_timing_recorder.h"
#include "timer.h"


namespace core::timing {

// Unified scoped timing helper used by TIMED_SCOPE.
class TimedScope {
   public:
#ifdef ENABLE_TIME_MEASUREMENTS
    TimedScope(StreamTimingRecorder& recorder, StreamTimingStage stage)
        : recorder_(&recorder), stage_(stage) {
        t_.restart();
    }

    TimedScope(StreamTimingRecorder& recorder, const char* name)
        : recorder_(&recorder), named_stage_(name) {
        t_.restart();
    }

    TimedScope(StreamTimingRecorder& recorder, const std::string& name)
        : recorder_(&recorder), named_stage_(name) {
        t_.restart();
    }

    TimedScope(double& accumulator, const char* name) : accumulator_(&accumulator) {
        (void)name;
        t_.restart();
    }

    TimedScope(double& accumulator, const std::string& name) : accumulator_(&accumulator) {
        (void)name;
        t_.restart();
    }

    ~TimedScope() {
        const double elapsed = t_.elapsed();
        if (recorder_ != nullptr) {
            if (!named_stage_.empty()) {
                recorder_->add_named(named_stage_, elapsed);
            } else if (stage_ == StreamTimingStage::Total) {
                recorder_->set(stage_, elapsed);
            } else {
                recorder_->add(stage_, elapsed);
            }
        } else if (accumulator_ != nullptr) {
            *accumulator_ += elapsed;
        }
    }

    double elapsed() {
        return t_.elapsed();
    }
#else
    TimedScope(StreamTimingRecorder& recorder, StreamTimingStage stage)
        : recorder_(&recorder), measure_total_(stage == StreamTimingStage::Total) {
        if (measure_total_) {
            t_.restart();
        }
    }

    TimedScope(StreamTimingRecorder&, const char*) {}

    TimedScope(StreamTimingRecorder&, const std::string&) {}

    TimedScope(double&, const char*) {}

    TimedScope(double&, const std::string&) {}

    ~TimedScope() {
        if (measure_total_ && recorder_ != nullptr) {
            recorder_->set(StreamTimingStage::Total, t_.elapsed());
        }
    }

    double elapsed() {
        return measure_total_ ? t_.elapsed() : 0.0;
    }
#endif

   private:
#ifdef ENABLE_TIME_MEASUREMENTS
    StreamTimingRecorder* recorder_ = nullptr;
    double* accumulator_ = nullptr;
    StreamTimingStage stage_ = StreamTimingStage::ReadInput;
    std::string named_stage_;
    timer t_;
#else
    StreamTimingRecorder* recorder_ = nullptr;
    bool measure_total_ = false;
    timer t_;
#endif
};

} // namespace core::timing


#define TIMED_SCOPE_JOIN_IMPL(a, b) a##b
#define TIMED_SCOPE_JOIN(a, b) TIMED_SCOPE_JOIN_IMPL(a, b)
#define TIMED_SCOPE(target, name_or_stage)                                            \
    ::core::timing::TimedScope TIMED_SCOPE_JOIN(_timed_scope_, __COUNTER__)((target), \
                                                                            (name_or_stage))

#endif /* CORE_SCOPED_STAGE_TIMER_H_ */
