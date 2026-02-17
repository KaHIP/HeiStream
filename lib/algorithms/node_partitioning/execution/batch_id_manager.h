#ifndef BATCH_ID_MANAGER_H
#define BATCH_ID_MANAGER_H

#include <condition_variable>
#include <mutex>
#include <queue>

#include "definitions.h"

struct ParsedLine {
    LongNodeID node_id;
    std::vector<LongNodeID> neighbors;
};

using BatchNode = std::pair<LongNodeID, std::vector<LongNodeID>>;

struct PartitionTask {
    int batch_id;

    std::vector<BatchNode> nodes;

    PartitionTask() : batch_id(-1) {}

    PartitionTask(int bid, std::vector<BatchNode> single_nodes)
        : batch_id(bid), nodes(std::move(single_nodes)) {}
};

class BatchIDManager {
   public:
    BatchIDManager(size_t max_batches) {
        max_active_batches = max_batches;
        for (size_t i = 0; i < max_active_batches; ++i) {
            free_ids.push(i);
        }
    }

    PartitionID get_batch_marker(size_t batch_id) const {
        return TO_BE_PARTITIONED - batch_id - 1;
    }

    // Last valid PartitionID, used to determine the maximum valid partition ID that is not a batch
    // marker.
    PartitionID get_max_valid_partition_id() const {
        return TO_BE_PARTITIONED - max_active_batches * 2;
    }

    // Hole eine freie batch_id (blockierend, wenn keine verfügbar)
    size_t acquire_id() {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]() { return !free_ids.empty(); });
        size_t id = free_ids.front();
        free_ids.pop();
        return id;
    }

    // Gib eine batch_id wieder frei
    void release_id(size_t id) {
        {
            std::lock_guard<std::mutex> lock(mtx);
            free_ids.push(id);
        }
        cv.notify_one();
    }

   private:
    std::queue<size_t> free_ids;
    std::mutex mtx;
    std::condition_variable cv;
    size_t max_active_batches;
};

#endif  // BATCH_ID_MANAGER_H
