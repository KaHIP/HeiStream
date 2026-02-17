/****************************************************************************
 * edge_partitioning_algorithm_factory.h
 *****************************************************************************/

#ifndef EDGE_PARTITIONING_ALGORITHM_FACTORY_H_
#define EDGE_PARTITIONING_ALGORITHM_FACTORY_H_

#include <memory>

class IStreamAlgorithm;


namespace edge_partitioning::factory {

std::unique_ptr<IStreamAlgorithm> create_algorithm();

} // namespace edge_partitioning::factory


#endif /* EDGE_PARTITIONING_ALGORITHM_FACTORY_H_ */
