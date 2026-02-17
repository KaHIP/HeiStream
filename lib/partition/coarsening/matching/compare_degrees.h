/******************************************************************************
 * compare_degrees.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef COMPARE_DEGREES_750FUZ7Z
#define COMPARE_DEGREES_750FUZ7Z

#include <vector>

#include "definitions.h"

class compare_degrees {
   public:
    explicit compare_degrees(std::vector<EdgeWeight>* degrees) : m_node_degrees(degrees) {}
    virtual ~compare_degrees() = default;

    bool operator()(const EdgeWeight left, const EdgeWeight right) {
        return (*m_node_degrees)[left] < (*m_node_degrees)[right];
    }

   private:
    std::vector<EdgeWeight>* m_node_degrees;
};

#endif /* end of include guard: COMPARE_DEGREES_750FUZ7Z */
