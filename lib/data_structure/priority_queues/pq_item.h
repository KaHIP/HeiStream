#ifndef PQ_ITEM_H
#define PQ_ITEM_H

#include <limits>
#include <vector>

#include "definitions.h"

const float INVALID_SCORE = std::numeric_limits<float>::max();
const unsigned INVALID_POS = std::numeric_limits<unsigned>::max();

class PQItem {
   public:
    float buffer_score;
    unsigned bucket_idx;
    unsigned num_adj_partitioned;
    std::vector<LongNodeID>* adjacents;
    unsigned pos_in_bucket;

    // Default constructor
    PQItem()
        : buffer_score(INVALID_SCORE),
          bucket_idx(0),
          num_adj_partitioned(0),
          adjacents(nullptr),
          pos_in_bucket(INVALID_POS) {}

    // Constructor for normal references (lvalue)
    PQItem(float score, std::vector<LongNodeID>& l, unsigned adj_part, unsigned pos_in_bucket,
           unsigned bucket_idx)
        : buffer_score(score),
          bucket_idx(bucket_idx),
          num_adj_partitioned(adj_part),
          pos_in_bucket(pos_in_bucket) {
        adjacents = new std::vector<LongNodeID>(l);
    }

    // Move constructor for rvalue references
    PQItem(float score, std::vector<LongNodeID>&& l, unsigned adj_part, unsigned pos_in_bucket,
           unsigned bucket_idx)
        : buffer_score(score),
          bucket_idx(bucket_idx),
          num_adj_partitioned(adj_part),
          pos_in_bucket(pos_in_bucket) {
        // Take ownership of the vector without copying
        adjacents = new std::vector<LongNodeID>(std::move(l));
    }

    // Copy constructor
    PQItem(const PQItem& other)
        : buffer_score(other.buffer_score),
          bucket_idx(other.bucket_idx),
          num_adj_partitioned(other.num_adj_partitioned),
          pos_in_bucket(other.pos_in_bucket) {
        // Copy neighbors if they exist
        adjacents = other.adjacents ? new std::vector<LongNodeID>(*other.adjacents) : nullptr;
    }

    // Move constructor
    PQItem(PQItem&& other) noexcept
        : buffer_score(other.buffer_score),
          bucket_idx(other.bucket_idx),
          num_adj_partitioned(other.num_adj_partitioned),
          adjacents(other.adjacents),
          pos_in_bucket(other.pos_in_bucket) {
        // Avoid double deletion
        other.adjacents = nullptr;
    }

    // Copy assignment operator
    PQItem& operator=(const PQItem& other) {
        if (this != &other) {
            buffer_score = other.buffer_score;
            bucket_idx = other.bucket_idx;
            num_adj_partitioned = other.num_adj_partitioned;
            pos_in_bucket = other.pos_in_bucket;

            // Delete old vector
            delete adjacents;

            // Copy new vector if it exists
            adjacents = other.adjacents ? new std::vector<LongNodeID>(*other.adjacents) : nullptr;
        }
        return *this;
    }

    // Move assignment operator
    PQItem& operator=(PQItem&& other) noexcept {
        if (this != &other) {
            buffer_score = other.buffer_score;
            bucket_idx = other.bucket_idx;
            num_adj_partitioned = other.num_adj_partitioned;
            pos_in_bucket = other.pos_in_bucket;

            // Delete old vector
            delete adjacents;

            // Take ownership of pointer
            adjacents = other.adjacents;
            other.adjacents = nullptr;
        }
        return *this;
    }

    // Destructor with nullptr check
    ~PQItem() {
        if (adjacents != nullptr) {
            delete adjacents;  // Safe for nullptr
        }
    }

    std::vector<LongNodeID>& get_adjacents() {
        return *adjacents;
    }

    unsigned get_address() const {
        return bucket_idx;
    }

    unsigned get_pos_in_bucket() const {
        return pos_in_bucket;
    }

    void set_pos_in_bucket(unsigned pos) {
        pos_in_bucket = pos;
    }

    void set_buffer_score(float new_buffer_score, unsigned new_bucket_idx) {
        buffer_score = new_buffer_score;
        bucket_idx = new_bucket_idx;
    }
};

#endif
