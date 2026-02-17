/******************************************************************************
 * maxNodeHeap.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef MAX_NODE_HEAP_39CK1B8I
#define MAX_NODE_HEAP_39CK1B8I

#include <execinfo.h>

#include <limits>
#include <unordered_map>
#include <vector>

#include "data_structure/priority_queues/priority_queue_interface.h"

typedef int Key;

template <typename Data>
class QElement {
   public:
    QElement(Data data, Key key, int index) : m_data(data), m_key(key), m_index(index) {};
    virtual ~QElement() {};

    Data& get_data() {
        return m_data;
    }
    const Data& get_data() const {
        return m_data;
    }

    void set_data(Data& data) {
        m_data = data;
    }

    Key get_key() {
        return m_key;
    }

    void set_key(Key key) {
        m_key = key;
    }

    int get_index() {
        return m_index;
    }

    void set_index(int index) {
        m_index = index;
    }

   private:
    Data m_data;
    Key m_key;
    int m_index;  // the index of the element in the heap
};

class maxNodeHeap : public priority_queue_interface {
   public:
    struct Data {
        NodeID node;
        Data(NodeID node) : node(node) {};
    };

    typedef QElement<Data> PQElement;

    maxNodeHeap() {};
    ~maxNodeHeap() override = default;

    NodeID size() override;
    bool empty() override;

    bool contains(NodeID node) override;
    void insert(NodeID id, Gain gain) override;

    NodeID deleteMax() override;
    void deleteNode(NodeID node) override;
    NodeID maxElement() override;
    Gain maxValue() override;

    void decreaseKey(NodeID node, Gain gain) override;
    void increaseKey(NodeID node, Gain gain) override;
    void changeKey(NodeID node, Gain gain) override;
    void subtractKey(NodeID node, Gain gain);
    Gain getKey(NodeID node) override;

   private:
    std::vector<PQElement> m_elements;  // elements that contain the data
    std::unordered_map<NodeID, int>
        m_element_index;                      // stores index of the node in the m_elements array
    std::vector<std::pair<Key, int>> m_heap;  // key and index in elements (pointer)

    bool higher_priority(int lhs, int rhs) const;
    void siftUp(int pos);
    void siftDown(int pos);
};

inline bool maxNodeHeap::higher_priority(int lhs, int rhs) const {
    const Key lhs_key = m_heap[lhs].first;
    const Key rhs_key = m_heap[rhs].first;
    if (lhs_key != rhs_key) {
        return lhs_key > rhs_key;
    }
    const NodeID lhs_node = m_elements[m_heap[lhs].second].get_data().node;
    const NodeID rhs_node = m_elements[m_heap[rhs].second].get_data().node;
    return lhs_node < rhs_node;
}

inline Gain maxNodeHeap::maxValue() {
    return m_heap[0].first;
};

inline NodeID maxNodeHeap::maxElement() {
    return m_elements[m_heap[0].second].get_data().node;
};

inline void maxNodeHeap::siftDown(int pos) {
    int lhsChild = 2 * pos + 1;
    int rhsChild = 2 * pos + 2;
    if (rhsChild < (int)m_heap.size()) {
        if (!higher_priority(lhsChild, pos) && !higher_priority(rhsChild, pos)) {
            return;  // we are done
        } else {
            // exchange with the larger one (maxHeap)
            int swap_pos = higher_priority(lhsChild, rhsChild) ? lhsChild : rhsChild;
            std::swap(m_heap[pos], m_heap[swap_pos]);

            int element_pos = m_heap[pos].second;
            m_elements[element_pos].set_index(pos);

            element_pos = m_heap[swap_pos].second;
            m_elements[element_pos].set_index(swap_pos);

            siftDown(swap_pos);
            return;
        }

    } else if (lhsChild < (int)m_heap.size()) {
        if (higher_priority(lhsChild, pos)) {
            std::swap(m_heap[pos], m_heap[lhsChild]);

            int element_pos = m_heap[pos].second;
            m_elements[element_pos].set_index(pos);

            element_pos = m_heap[lhsChild].second;
            m_elements[element_pos].set_index(lhsChild);

            siftDown(lhsChild);
            return;
        } else {
            return;  // we are done
        }
    }
}

inline void maxNodeHeap::siftUp(int pos) {
    if (pos > 0) {
        int parentPos = (int)(pos - 1) / 2;
        if (higher_priority(pos, parentPos)) {
            // heap condition not fulfulled
            std::swap(m_heap[parentPos], m_heap[pos]);

            int element_pos = m_heap[pos].second;
            m_elements[element_pos].set_index(pos);

            // update the heap index in the element
            element_pos = m_heap[parentPos].second;
            m_elements[element_pos].set_index(parentPos);

            siftUp(parentPos);
        }
    }
}

inline NodeID maxNodeHeap::size() {
    return m_heap.size();
}

inline bool maxNodeHeap::empty() {
    return m_heap.empty();
}

inline void maxNodeHeap::insert(NodeID node, Gain gain) {
    if (m_element_index.find(node) == m_element_index.end()) {
        int element_index = m_elements.size();
        int heap_size = m_heap.size();

        m_elements.push_back(PQElement(Data(node), gain, heap_size));
        m_heap.push_back(std::pair<Key, int>(gain, element_index));
        m_element_index[node] = element_index;
        siftUp(heap_size);
    }
}

inline void maxNodeHeap::deleteNode(NodeID node) {
    int element_index = m_element_index[node];
    int heap_index = m_elements[element_index].get_index();

    m_element_index.erase(node);

    std::swap(m_heap[heap_index], m_heap[m_heap.size() - 1]);
    // update the position of its element in the element array
    m_elements[m_heap[heap_index].second].set_index(heap_index);

    // we dont want holes in the elements array -- delete the deleted element from the array
    if (element_index != (int)(m_elements.size() - 1)) {
        std::swap(m_elements[element_index], m_elements[m_elements.size() - 1]);
        m_heap[m_elements[element_index].get_index()].second = element_index;
        int cnode = m_elements[element_index].get_data().node;
        m_element_index[cnode] = element_index;
    }

    m_elements.pop_back();
    m_heap.pop_back();

    if (m_heap.size() > 1 && heap_index < (int)m_heap.size()) {
        // fix the max heap property
        siftDown(heap_index);
        siftUp(heap_index);
    }
};

inline NodeID maxNodeHeap::deleteMax() {
    if (m_heap.size() > 0) {
        int element_index = m_heap[0].second;
        int node = m_elements[element_index].get_data().node;
        m_element_index.erase(node);

        m_heap[0] = m_heap[m_heap.size() - 1];
        // update the position of its element in the element array
        m_elements[m_heap[0].second].set_index(0);

        // we dont want holes in the elements array -- delete the deleted element from the array
        if (element_index != (int)(m_elements.size() - 1)) {
            m_elements[element_index] = m_elements[m_elements.size() - 1];
            m_heap[m_elements[element_index].get_index()].second = element_index;
            int cnode = m_elements[element_index].get_data().node;
            m_element_index[cnode] = element_index;
        }

        m_elements.pop_back();
        m_heap.pop_back();

        if (m_heap.size() > 1) {
            // fix the heap property
            siftDown(0);
        }

        return node;
    }

    return -1;
}

inline void maxNodeHeap::changeKey(NodeID node, Gain gain) {
    Gain old_gain = m_heap[m_elements[m_element_index[node]].get_index()].first;
    if (old_gain > gain) {
        decreaseKey(node, gain);
    } else if (old_gain < gain) {
        increaseKey(node, gain);
    }
};

inline void maxNodeHeap::decreaseKey(NodeID node, Gain gain) {
    ASSERT_TRUE(m_element_index.find(node) != m_element_index.end());
    int queue_idx = m_element_index[node];
    int heap_idx = m_elements[queue_idx].get_index();
    m_elements[queue_idx].set_key(gain);
    m_heap[heap_idx].first = gain;
    siftDown(heap_idx);
}

inline void maxNodeHeap::increaseKey(NodeID node, Gain gain) {
    ASSERT_TRUE(m_element_index.find(node) != m_element_index.end());
    int queue_idx = m_element_index[node];
    int heap_idx = m_elements[queue_idx].get_index();
    m_elements[queue_idx].set_key(gain);
    m_heap[heap_idx].first = gain;
    siftUp(heap_idx);
}

inline void maxNodeHeap::subtractKey(NodeID node, Gain gain) {
    ASSERT_TRUE(m_element_index.find(node) != m_element_index.end());
    int queue_idx = m_element_index[node];
    int heap_idx = m_elements[queue_idx].get_index();
    // std::cout << "Set " << m_heap[heap_idx].first << " to "
    //           << m_heap[heap_idx].first + gain << '\n';
    m_elements[queue_idx].set_key(m_heap[heap_idx].first + gain);
    m_heap[heap_idx].first = m_heap[heap_idx].first + gain;
    siftDown(heap_idx);
}

inline Gain maxNodeHeap::getKey(NodeID node) {
    return m_heap[m_elements[m_element_index[node]].get_index()].first;
};

inline bool maxNodeHeap::contains(NodeID node) {
    return m_element_index.find(node) != m_element_index.end();
}

#endif
