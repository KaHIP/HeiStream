//
// Created by adilchhabra on 16.05.24.
//

#pragma once

#include <vector>

#include "data_structure/compression_vectors/CompressionDataStructure.h"
#include "cpi/run_length_compression.hpp"

// Template class for Run-Length Compression Vector
template <typename T, bool enable_buffer = false>
class RunLengthCompressionVector : public CompressionDataStructure<T> {
private:
  // Instance of the cpi RunLengthCompression class with template parameters
  // for the type T, and fixed byte sizes 16 and 64 uncompressed runs in memory
  cpi::RunLengthCompression<T, 16, 64, enable_buffer> compressed_data_;

public:
  // Override the Append method to add a value to the compressed data structure
  void Append(T value) override {
    compressed_data_.push_back(value);
  }

  void BatchAppend(int batch_index, T value) override {
    compressed_data_.push_back(value);
  }

  // Override the GetValueByIndex method to retrieve a value by index from the
  // compressed data structure
  T GetValueByIndex(int index) const override {
    return compressed_data_[index];
  }

  T GetValueByBatchIndex(int batch_index, int data_index) const override {
    return compressed_data_[data_index];
  }

};
