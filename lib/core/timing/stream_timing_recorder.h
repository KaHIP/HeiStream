/****************************************************************************
 * stream_timing_recorder.h
 *****************************************************************************/

#ifndef CORE_STREAM_TIMING_RECORDER_H_
#define CORE_STREAM_TIMING_RECORDER_H_

#include <array>
#include <cassert>
#include <cstddef>
#include <string>
#include <unordered_map>


namespace core::timing {

enum class StreamTimingStage : std::size_t {
    ReadInput = 0,
    ModelBuild,
    Partition,
    Postprocess,
    Total,
    BufferMaintenance,
    Count
};

class StreamTimingRecorder {
   public:
    struct Snapshot {
        double read_input = 0.0;
        double model_build = 0.0;
        double partition = 0.0;
        double postprocess = 0.0;
        double total = 0.0;
        double buffer_maintenance = 0.0;
    };

    void reset() {
        values_.fill(0.0);
        written_.fill(false);
        named_values_.clear();
    }

    void set(StreamTimingStage stage, double seconds) {
        values_[index(stage)] = seconds;
        written_[index(stage)] = true;
    }

    void add(StreamTimingStage stage, double seconds) {
        values_[index(stage)] += seconds;
        written_[index(stage)] = true;
    }

    void add_delta(StreamTimingStage stage, double before_seconds, double after_seconds) {
        add(stage, after_seconds - before_seconds);
    }

    double get(StreamTimingStage stage) const {
        return values_[index(stage)];
    }

    bool has(StreamTimingStage stage) const {
        return written_[index(stage)];
    }

    void add_named(const std::string& name, double seconds) {
        named_values_[name] += seconds;
    }

    double get_named(const std::string& name) const {
        auto it = named_values_.find(name);
        if (it == named_values_.end()) {
            return 0.0;
        }
        return it->second;
    }

    Snapshot snapshot() const {
        Snapshot s;
        s.read_input = get(StreamTimingStage::ReadInput);
        s.model_build = get(StreamTimingStage::ModelBuild);
        s.partition = get(StreamTimingStage::Partition);
        s.postprocess = get(StreamTimingStage::Postprocess);
        s.total = get(StreamTimingStage::Total);
        s.buffer_maintenance = get(StreamTimingStage::BufferMaintenance);
        return s;
    }

    bool validate_totals() const {
        const double non_total =
            get(StreamTimingStage::ReadInput) + get(StreamTimingStage::ModelBuild) +
            get(StreamTimingStage::Partition) + get(StreamTimingStage::Postprocess);
        return get(StreamTimingStage::Total) + 1e-9 >= non_total;
    }

   private:
    static constexpr std::size_t index(StreamTimingStage stage) {
        return static_cast<std::size_t>(stage);
    }

    std::array<double, static_cast<std::size_t>(StreamTimingStage::Count)> values_{};
    std::array<bool, static_cast<std::size_t>(StreamTimingStage::Count)> written_{};
    std::unordered_map<std::string, double> named_values_{};
};

} // namespace core::timing


#endif /* CORE_STREAM_TIMING_RECORDER_H_ */
