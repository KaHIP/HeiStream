/****************************************************************************
 * node_partitioning_algorithm_factory.h
 *****************************************************************************/

#ifndef NODE_PARTITIONING_ALGORITHM_FACTORY_H_
#define NODE_PARTITIONING_ALGORITHM_FACTORY_H_

#include <memory>

class IStreamAlgorithm;


namespace node_partitioning::factory {

std::unique_ptr<IStreamAlgorithm> create_algorithm();

} // namespace node_partitioning::factory


#endif /* NODE_PARTITIONING_ALGORITHM_FACTORY_H_ */
