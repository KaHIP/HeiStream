/******************************************************************************
 * buffered_map.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BUFFERED_MAP_EFRXO4X2
#define BUFFERED_MAP_EFRXO4X2

#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "definitions.h"

//////// Vector Version
class buffered_map {
   public:
    buffered_map(std::vector<PartitionID>* external_vector, bool vec_preassigned) {
        this->key_to_value = external_vector;
        this->vec_preassigned = vec_preassigned;
    }
    ~buffered_map() {}

    bool has_key(LongNodeID key) {
        return (this->vec_preassigned) || (*this->key_to_value)[key] != INVALID_PARTITION;
    }

    // We make no test here. Programmer has to call has_key first
    LongNodeID const operator[](LongNodeID key) const {
        return (*this->key_to_value)[key];
    }

    // We make no test here. has_key has to return false before push_back and value ==
    // this->value_to_key.size()
    void push_back(LongNodeID key, LongNodeID value) {
        (*this->key_to_value)[key] = value;
        this->value_to_key.push_back(key);
    }

    void clear() {
        if (!this->vec_preassigned) {
            for (auto& key : this->value_to_key) {
                (*this->key_to_value)[key] = INVALID_PARTITION;
            }
        }
        this->value_to_key.clear();
    }

    bool vec_preassigned;
    std::vector<PartitionID>* key_to_value;
    std::vector<LongNodeID> value_to_key;
};

// macros
#define forall_pairs_bufferedmap(obj, key, value)                                \
    {                                                                            \
        for (LongNodeID value = 0; value < (obj).value_to_key.size(); value++) { \
            LongNodeID key = (obj).value_to_key[value];
#define endfor_bufferedmap \
    }                      \
    }

class buffered_input {
   public:
    buffered_input(std::vector<std::string>* input_lines) {
        this->lines = input_lines;
        this->row = 0;
        this->n_buffer_lines = 1;
        this->column = 0;
    }
    buffered_input(std::vector<std::string>* input_lines, LongNodeID input_cursor,
                   LongNodeID n_buffer_lines) {
        this->lines = input_lines;
        this->row = input_cursor;
        this->n_buffer_lines = n_buffer_lines;
        this->column = 0;
    }
    ~buffered_input() {}
    void simple_scan_line(std::vector<LongNodeID>& vec) {
        vec.clear();

        const std::string& line = (*this->lines)[this->row];
        const char* ptr = line.c_str();
        const char* end = ptr + line.length();

        while (ptr < end) {
            // Skip non-digits
            while (ptr < end && (*ptr < '0' || *ptr > '9')) {
                ptr++;
            }

            if (ptr >= end)
                break;

            LongNodeID num = 0;
            while (ptr < end && *ptr >= '0' && *ptr <= '9') {
                num = num * 10 + (*ptr - '0');
                ptr++;
            }

            vec.push_back(num);
        }
        this->row++;
    }

   private:
    std::vector<std::string>* lines;
    LongNodeID row, column, n_buffer_lines;
};

#endif /* end of include guard: BUFFERED_MAP_EFRXO4X2 */
// macros
