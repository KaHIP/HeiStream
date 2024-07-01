//
// Created by adilchhabra on 16.05.24.
//

#pragma once

#include <vector>

#include "data_structure/compression_vectors/CompressionDataStructure.h"
#include "data_structure/compression_vectors/RunLengthCompressionVector.h"
#include "cpi/run_length_compression.hpp"


template <typename T, bool enable_buffer = false>
class BatchRunLengthCompression : public CompressionDataStructure<T> {
private:
    std::vector<RunLengthCompressionVector<T, enable_buffer>> batch_compression_;

public:
    BatchRunLengthCompression(int x) : batch_compression_(x) {}

    void Append(T value) override {
      // not applicable since batch_index not specified
      throw std::logic_error("Append without batch index is undefined for BatchRunLengthCompression");
    }

    void BatchAppend(int batch_index, T value) override {
        batch_compression_[batch_index].Append(value);
    }

    T GetValueByIndex(int index) const override {
      // not applicable since batch_index not specified
      throw std::logic_error("GetValueByIndex without batch index is undefined for BatchRunLengthCompression");
    }

    T GetValueByBatchIndex(int batch_index, int data_index) const {
        if (batch_index < batch_compression_.size()) {
            return batch_compression_[batch_index].GetValueByIndex(data_index);
        }
        return T();
    }
};