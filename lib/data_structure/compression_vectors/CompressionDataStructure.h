//
// Created by adilchhabra on 16.05.24.
//

#pragma once

template <typename T, bool enable_buffer = false>
class CompressionDataStructure {
public:
  virtual void Append(T value) = 0;
  virtual void BatchAppend(int batch_index, T value) = 0;
  virtual T GetValueByIndex(int index) const = 0;
  virtual T GetValueByBatchIndex(int batch_index, int data_index) const = 0;
};
