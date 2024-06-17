/******************************************************************************
 * graph_io_stream.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <math.h>
#include <sstream>
#include "graph_io_stream.h"
#include "timer.h"
#include "definitions.h"
#include <chrono>
#include <cmath>
#include <iomanip>

#define MIN(A, B) ((A)<(B))?(A):(B)
#define MAX(A, B) ((A)>(B))?(A):(B)


graph_io_stream::graph_io_stream() {

}

graph_io_stream::~graph_io_stream() {

}


NodeID
graph_io_stream::createModel(PartitionConfig &config, graph_access &G, std::vector <std::vector<LongNodeID>> *&input) {
    NodeWeight total_nodeweight = 0;
    NodeID node_counter = 0;
    EdgeID edge_counter = 0;
    LongEdgeID used_edges = 0;
    bool read_ew = false;
    bool read_nw = false;
    LongEdgeID nmbEdges = 2 * config.remaining_stream_edges;
    LongNodeID target;
    NodeWeight weight;
    std::vector < std::vector < std::pair < NodeID, EdgeWeight>>> all_edges;
    std::vector <NodeWeight> all_nodes;

    std::vector <NodeWeight> all_assigned_ghost_nodes(config.nmbNodes + config.quotient_nodes, 0);
    all_edges.resize(config.nmbNodes + config.quotient_nodes);
    all_nodes.resize(config.nmbNodes + config.quotient_nodes);
    config.lower_global_node = config.total_stream_nodecounter + 1; // Bounds below start from 1 instead of 0
    config.upper_global_node = config.total_stream_nodecounter + config.nmbNodes;
    LongNodeID cursor = 0;
    NodeID node = 0;

    config.curr_batch++;

    if (config.ram_stream) {
        cursor = input->size() - config.remaining_stream_nodes;
    }


    if (nmbEdges > std::numeric_limits<EdgeWeight>::max() || config.nmbNodes > std::numeric_limits<LongNodeID>::max()) {
#ifdef MODE64BITEDGES
        std::cerr <<  "The graph is too large. Currently only 64bits supported!"  << std::endl;
#else
        std::cerr << "The graph is too large. Currently only 32bits supported!" << std::endl;
#endif
        exit(0);
    }

    switch (config.remaining_stream_ew) {
        case 1:
            read_ew = true;
            break;
        case 10:
            read_nw = true;
            break;
        case 11:
            read_ew = true;
            read_nw = true;
            break;
    }

    /* config.degree_nodeBlock = new std::vector<std::vector<EdgeWeight>> (config.nmbNodes, std::vector<EdgeWeight>(config.k,0)); */
    /* config.edge_block_nodes = new std::vector<std::vector<NodeID>> (config.k, std::vector<NodeID>()); */
    config.edge_block_nodes = new std::vector <std::vector<std::pair < NodeID, NodeWeight>> >
                              (config.k, std::vector < std::pair < NodeID, NodeWeight >> ());

    setupForGhostNeighbors(config);

    for (node_counter = 0; node_counter < config.nmbNodes; node_counter++) {
        std::vector <LongNodeID> &line_numbers = (*input)[cursor];
        LongNodeID col_counter = 0;
        node = (NodeID)
        node_counter;
        weight = 1;
        if (read_nw) {
            weight = line_numbers[col_counter++];
            if (total_nodeweight > std::numeric_limits<NodeWeight>::max()) {
                std::cerr << "The sum of the node weights is too large (it exceeds the node weight type)." << std::endl;
                std::cerr << "Currently not supported. Please scale your node weights." << std::endl;
                exit(0);
            }
        }
        total_nodeweight += weight;
        processNodeWeight(config, all_nodes, node, weight);

        while (col_counter < line_numbers.size()) {
            target = line_numbers[col_counter++];
            EdgeWeight edge_weight = 1;
            if (read_ew) {
                edge_weight = line_numbers[col_counter++];
            }

            if (target > config.upper_global_node) { // edge to future batch
                processGhostNeighborInBatch(config, node, target, edge_weight);
            } else if (target < config.lower_global_node) { // edge to previous batch
                used_edges++;
                processQuotientEdgeInBatch(config, node, target, edge_weight);
            } else { // edge to current batch
                used_edges += ((NodeID)(target - config.lower_global_node) <
                               node); // used_edges only counts arcs to previus nodes
                edge_counter += insertRegularEdgeInBatch(config, all_edges, node, target, edge_weight);
            }
        }

        cursor++;
    }
    if (!config.ram_stream) {
        delete input;
    }

    NodeID uncontracted_ghost_nodes = mapGhostKeysToNodesInBatch(config, all_edges, all_nodes, all_assigned_ghost_nodes,
                                                                 node_counter);
    insertQuotientNodesInBatch(config, all_nodes, uncontracted_ghost_nodes, node_counter);
    edge_counter += insertGhostEdgesInBatch(config, all_edges);
    edge_counter += insertQuotientEdgesInBatch(config, all_edges, uncontracted_ghost_nodes);

    createGraphForBatch(config, G, node_counter, edge_counter, all_edges, all_nodes, all_assigned_ghost_nodes);

    /* delete config.degree_nodeBlock; */
    delete config.edge_block_nodes;


    config.total_stream_nodecounter += config.nmbNodes;
    config.total_stream_nodeweight += total_nodeweight;
    config.remaining_stream_nodes -= config.nmbNodes;
    config.remaining_stream_edges -= used_edges;

    if (node_counter != (NodeID) config.nmbNodes + uncontracted_ghost_nodes + config.quotient_nodes) {
        std::cerr << "number of specified nodes mismatch" << std::endl;
        std::cerr << (config.nmbNodes + uncontracted_ghost_nodes + config.quotient_nodes) << " " << node_counter
                  << std::endl;
        exit(0);
    }

    return node_counter;
}


void graph_io_stream::insertQuotientNodesInBatch(PartitionConfig &config, std::vector <NodeWeight> &all_nodes,
                                                 NodeID uncontracted_ghost_nodes, NodeID &node_counter) {
    while (node_counter < config.nmbNodes + uncontracted_ghost_nodes +
                          config.quotient_nodes) { // Add artificial nodes representing current blocks
        NodeID node = node_counter++;
        NodeID targetPar = node - (NodeID)
        config.nmbNodes - uncontracted_ghost_nodes;
        NodeWeight weight = (*config.stream_blocks_weight)[targetPar];
        all_nodes[node] = weight;
    }
}


/* EdgeID graph_io_stream::insertQuotientEdgesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, NodeID uncontracted_ghost_nodes) { */
/* 	EdgeID inserted_edges = 0; */
/* 	for (NodeID node=0; node < config.nmbNodes; node++) { */
/* 		for (PartitionID block=0; block < config.k; block++) { */
/* 			NodeID target = config.nmbNodes + uncontracted_ghost_nodes + block; */
/* 			EdgeWeight edge_weight = (*config.degree_nodeBlock)[node][block]; */
/* 			if (edge_weight > 0) { */
/* 				inserted_edges += includeEdgeInBatch(all_edges, node, target, (1+config.double_non_ghost_edges)*edge_weight); */
/* 			} */
/* 		} */
/* 	} */
/* 	return inserted_edges; */
/* } */
/* EdgeID graph_io_stream::insertQuotientEdgesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, NodeID uncontracted_ghost_nodes) { */
/* 	EdgeID inserted_edges = 0; */
/* 	NodeID prev_node; */
/* 	EdgeWeight edge_weight; */
/* 	for (PartitionID block=0; block < config.k; block++) { */
/* 		NodeID target = config.nmbNodes + uncontracted_ghost_nodes + block; */
/* 		if ((*config.edge_block_nodes)[block].size() >= 1) { */
/* 			prev_node = (*config.edge_block_nodes)[block][0]; */
/* 			edge_weight = 0; */
/* 		} else continue; // if block has no neighbors in batch, continue outer for loop */
/* 		for (auto& node : (*config.edge_block_nodes)[block]) { */
/* 			if (node == prev_node) { */
/* 				edge_weight++; */
/* 			} else { */
/* 				inserted_edges += includeEdgeInBatch(all_edges, prev_node, target, (1+config.double_non_ghost_edges)*edge_weight); */
/* 				prev_node = node; */
/* 				edge_weight = 1; */
/* 			} */
/* 		} */
/* 		// necessary because last neighbor of each block is not directly included inside the inner for loop */
/* 		inserted_edges += includeEdgeInBatch(all_edges, prev_node, target, (1+config.double_non_ghost_edges)*edge_weight); */
/* 	} */
/* 	return inserted_edges; */
/* } */
EdgeID graph_io_stream::insertQuotientEdgesInBatch(PartitionConfig &config,
                                                   std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

>& all_edges,
NodeID uncontracted_ghost_nodes
) {
EdgeID inserted_edges = 0;
EdgeWeight edge_weight;
for (
PartitionID block = 0;
block<config.
k;
block++) {
NodeID target = config.nmbNodes + uncontracted_ghost_nodes + block;
if ((*config.edge_block_nodes)[block].

size()

< 1) continue; // if block has no neighbors in batch, continue outer for loop
for (
auto &element
: (*config.edge_block_nodes)[block]) {
inserted_edges +=
includeEdgeInBatch(all_edges, element
.first,
target,
(1+config.double_non_ghost_edges)*element.second);
}}
return
inserted_edges;
}


EdgeID graph_io_stream::insertGhostEdgesInBatch(PartitionConfig &config,
                                                std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

>& all_edges) {
EdgeID inserted_edges = 0;
if (config.stream_allow_ghostnodes || config.restream_number) {
for (
LongNodeID ghost_key = 0;
ghost_key<config.ghostkey_to_edges->

size();

ghost_key++) {
NodeID ghost_node = (*config.ghostkey_to_node)[ghost_key];
auto &neighbors_list = (*config.ghostkey_to_edges)[ghost_key];
for (
auto &edge
: neighbors_list) {
inserted_edges +=
includeEdgeInBatch(all_edges, ghost_node, edge
.first, edge.second);
}}}
return
inserted_edges;
}

NodeID graph_io_stream::mapGhostKeysToNodesInBatch(PartitionConfig &config,
                                                   std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

>& all_edges,
std::vector <NodeWeight> &all_nodes, std::vector<NodeWeight>
& all_assigned_ghost_nodes,
NodeID &node_counter
) {
NodeID inserted_nodes = 0;
if (config.restream_number) {
inserted_nodes = restreamMapGhostKeysToNodes(config);
config.ghostglobal_to_ghostkey->

clear();

} else if (config.stream_allow_ghostnodes) {
inserted_nodes = greedyMapGhostKeysToNodes(config, all_edges, all_nodes, all_assigned_ghost_nodes, node_counter);
config.ghostglobal_to_ghostkey->

clear();

}
return
inserted_nodes;
}

NodeID graph_io_stream::restreamMapGhostKeysToNodes(PartitionConfig &config) {
    if (config.ghostkey_to_node != NULL) {
        delete config.ghostkey_to_node;
    }
    config.ghostkey_to_node = new std::vector<NodeID>(config.ghost_nodes, 0);
    for (PartitionID targetPar = 0; targetPar < config.ghost_nodes; targetPar++) {
        NodeID node = targetPar + config.nmbNodes;
        (*config.ghostkey_to_node)[targetPar] = node; // PS: weight of quotient nodes is already up-to-date
    }
    return 0;
}

NodeID graph_io_stream::greedyMapGhostKeysToNodes(PartitionConfig &config,
                                                  std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

>& all_edges,
std::vector <NodeWeight> &all_nodes, std::vector<NodeWeight>
& all_assigned_ghost_nodes,
NodeID &node_counter
) {
if (config.ghostkey_to_node != NULL)
delete config.
ghostkey_to_node;
config.
ghostkey_to_node = new std::vector<NodeID>(config.ghost_nodes, 0);
NodeID inserted_nodes = MIN(config.ghost_nodes, config.ghost_nodes_threshold);
all_nodes.
resize(config
.nmbNodes + inserted_nodes + config.quotient_nodes);
all_edges.
resize(config
.nmbNodes + inserted_nodes + config.quotient_nodes);
LongNodeID ghost_key;
for (
ghost_key = 0;
ghost_key<config.
ghost_nodes;
ghost_key++) {
if (ghost_key < (LongNodeID) inserted_nodes) { // uncontracted ghost node
NodeID node = node_counter++;
(*config.ghostkey_to_node)[ghost_key] =
node;
all_nodes[node] = 1;
} else { // contracted ghost node
auto &list_edges = (*config.ghostkey_to_edges)[ghost_key];
//			int contr_pos = random_functions::nextInt(0, list_edges.size()-1);
LongNodeID contr_pos = ((LongNodeID) ghost_key * config.curr_batch) % list_edges.size();
NodeID node = return_and_delete_element(list_edges, contr_pos).first; // remove self edge
(*config.ghostkey_to_node)[ghost_key] =
node;
all_nodes[node]++;
all_assigned_ghost_nodes[node]++;
}}
return
inserted_nodes;
}

void graph_io_stream::processGhostNeighborInBatch(PartitionConfig &config, NodeID node, LongNodeID ghost_target,
                                                  EdgeWeight edge_weight) {
    LongNodeID ghost_key;
    PartitionID targetGlobalPar = (*config.stream_nodes_assign)[ghost_target - 1];
    if (config.stream_allow_ghostnodes || (config.restream_number && targetGlobalPar != INVALID_PARTITION)) {
        if (!config.ghostglobal_to_ghostkey->has_key(ghost_target - 1)) {
            ghost_key = config.ghost_nodes++;
            config.ghostglobal_to_ghostkey->push_back(ghost_target - 1, ghost_key);
            config.ghostkey_to_edges->push_back(std::vector < std::pair < NodeID, ShortEdgeWeight >> ());
        } else {
            ghost_key = (*config.ghostglobal_to_ghostkey)[ghost_target - 1];
        }
        (*config.ghostkey_to_edges)[ghost_key].push_back(std::make_pair(node, (ShortEdgeWeight) edge_weight));
    }
}


/* void graph_io_stream::processQuotientEdgeInBatch(PartitionConfig & config, NodeID node, LongNodeID global_target, EdgeWeight edge_weight) { */
/* 	PartitionID targetGlobalPar = (*config.stream_nodes_assign)[global_target-1]; */
/* 	if ( targetGlobalPar == config.k ) { // SIGNAL: neighbor not yet assigned */
/* 		processGhostNeighborInBatch(config, node, global_target, edge_weight); */
/* 		return; */
/* 	} */
/* 	if ( targetGlobalPar > config.k ) { */
/* 		std::cerr << "ERROR regarding number of parts.\n"; */
/* 		exit(0); */
/* 	} */
/* 	(*config.degree_nodeBlock)[node][targetGlobalPar] += edge_weight; */
/* } */
/* void graph_io_stream::processQuotientEdgeInBatch(PartitionConfig & config, NodeID node, LongNodeID global_target, EdgeWeight edge_weight) { */
/* 	PartitionID targetGlobalPar = (*config.stream_nodes_assign)[global_target-1]; */
/* 	if ( targetGlobalPar == config.k ) { // SIGNAL: neighbor not yet assigned */
/* 		processGhostNeighborInBatch(config, node, global_target, edge_weight); */
/* 		return; */
/* 	} */
/* 	if ( targetGlobalPar > config.k ) { */
/* 		std::cerr << "ERROR regarding number of parts.\n"; */
/* 		exit(0); */
/* 	} */
/* 	for (EdgeWeight i=0; i<edge_weight; i++) (*config.edge_block_nodes)[targetGlobalPar].push_back(node); */
/* } */
void graph_io_stream::processQuotientEdgeInBatch(PartitionConfig &config, NodeID node, LongNodeID global_target,
                                                 EdgeWeight edge_weight) {
    PartitionID targetGlobalPar = (*config.stream_nodes_assign)[global_target - 1];
    if (targetGlobalPar == config.k) { // SIGNAL: neighbor not yet assigned
        processGhostNeighborInBatch(config, node, global_target, edge_weight);
        return;
    }
    if (targetGlobalPar > config.k) {
        std::cerr << "ERROR regarding number of parts.\n";
        exit(0);
    }
    if ((*config.edge_block_nodes)[targetGlobalPar].size() >= 1) {
        auto &curr_element = (*config.edge_block_nodes)[targetGlobalPar][0];
        if (curr_element.first == node) {
            curr_element.second += edge_weight;
            return;
        }
    }
    std::pair <NodeID, NodeWeight> element = {node, edge_weight};
    (*config.edge_block_nodes)[targetGlobalPar].push_back(element);
}


void graph_io_stream::processNodeWeight(PartitionConfig &config, std::vector <NodeWeight> &all_nodes, NodeID node,
                                        NodeWeight weight) {
    if (config.edge_partition) {
        all_nodes[node] = weight;
        return;
    }
    LongNodeID global_node = config.lower_global_node + (LongNodeID)
    node - 1;
    PartitionID nodeGlobalPar = (*config.stream_nodes_assign)[global_node];
    if (config.restream_number) {
        (*config.stream_blocks_weight)[nodeGlobalPar] -= weight;
    }
    all_nodes[node] = weight;
}


EdgeID graph_io_stream::insertRegularEdgeInBatch(PartitionConfig &config,
                                                 std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

>& all_edges,
NodeID node, LongNodeID
global_target,
EdgeWeight edge_weight
) {
NodeID target = (NodeID) (global_target - config.lower_global_node);
if(target == node) {
std::cerr <<  "The graph file contains self-loops, which are not supported. Please remove them from the file."  <<
std::endl;
exit(0);
}
if (target > node) {
return 0;
}

return
includeEdgeInBatch(all_edges, node, target,
(1+config.double_non_ghost_edges)*edge_weight);
}


EdgeID graph_io_stream::includeEdgeInBatch(std::vector < std::vector < std::pair < NodeID, EdgeWeight
>>>& all_edges,
NodeID node, NodeID
target,
EdgeWeight edge_weight
) {
all_edges[node].
push_back(std::make_pair(target, edge_weight)
);
all_edges[target].
push_back(std::make_pair(node, edge_weight)
);

return 2;
}


void graph_io_stream::setupForGhostNeighbors(PartitionConfig &config) {
    if (config.stream_allow_ghostnodes || config.restream_number) {
        if (config.ghostglobal_to_ghostkey == NULL) {
            config.ghostglobal_to_ghostkey = new buffered_map(config.stream_nodes_assign, config.restream_number);
        }
        if (config.ghostkey_to_edges != NULL) {
            config.ghostkey_to_edges->clear();
        } else {
            config.ghostkey_to_edges = new std::vector <std::vector<std::pair < NodeID, ShortEdgeWeight>> > ();
        }
        config.ghost_nodes = 0 + (config.restream_number > 0) * (config.k);
        if (config.restream_number > 0) {
            config.ghost_nodes = config.k;
            config.ghostkey_to_edges->resize(config.ghost_nodes);
        }
    }
}


void
graph_io_stream::createGraphForBatch(PartitionConfig &config, graph_access &G, NodeID node_counter, EdgeID edge_counter,
                                     std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

>& all_edges,
std::vector <NodeWeight> &all_nodes, std::vector<NodeWeight>
& all_assigned_ghost_nodes) {
G.
stream_repeat_construction(node_counter, edge_counter
);
G.
resizeImplicitGhostNodesList(node_counter);
std::vector <EdgeWeight> neighbor_edgeweight(node_counter, 0);
std::vector <EdgeID> neighbor_edgeid(node_counter, 0);
std::vector <NodeID> neighbors;
EdgeID e;
for(
NodeID i = 0;
i<node_counter;
i++) { // actual insertion of nodes and edges
NodeID node = G.new_node();
NodeWeight weight = all_nodes[node];
NodeWeight ghost_weight = all_assigned_ghost_nodes[node];
G.
setNodeWeight(node, weight
);
G.
setImplicitGhostNodes(node, ghost_weight
);
G.
setPartitionIndex(node,
0);
recoverBlockAssignedToNode(config, G, node, node_counter
);
for (auto &[target, edge_weight] : all_edges[node]) {
if (neighbor_edgeweight[target] == 0) {
e = G.new_edge(node, target);
neighbor_edgeid[target] =
e;
neighbors.
push_back(target);
} else {
e = neighbor_edgeid[target];
}
neighbor_edgeweight[target] +=
edge_weight;
G.
setEdgeWeight(e, neighbor_edgeweight[target]
);
}
for (
auto &target
: neighbors) {
neighbor_edgeweight[target] = 0;
}
neighbors.

clear();

}
G.

stream_finish_construction();

}


void graph_io_stream::recoverBlockAssignedToNode(PartitionConfig &config, graph_access &G, NodeID node,
                                                 NodeID node_counter) {
    if (node >= node_counter - config.quotient_nodes) { // quotient nodes
        NodeID targetPar = node - (node_counter - config.quotient_nodes);
        G.setPartitionIndex(node, targetPar);
    } else if (config.restream_number &&
               (config.restream_vcycle || config.initial_partitioning_type == INITIAL_PARTITIONING_FENNEL)) {
        if (node < config.nmbNodes) { // regular nodes
            LongNodeID global_node = config.lower_global_node + (LongNodeID)
            node - 1;
            PartitionID targetPar = (*config.stream_nodes_assign)[global_node];
            G.setPartitionIndex(node, targetPar);
        } else { // ghost node. PS: there are no isolated ghost nodes in restream
            std::cerr << "Unexpected branch.\n";
            exit(0);
        }
    }
}


void graph_io_stream::generalizeStreamPartition(PartitionConfig &config, graph_access &G_local) {
    for (NodeID node = 0, end = config.nmbNodes; node < end; node++) {
        PartitionID block = G_local.getPartitionIndex(node);
        LongNodeID global_node = (LongNodeID)
        node + config.lower_global_node - 1;
        (*config.stream_nodes_assign)[global_node] = block;
        (*config.stream_blocks_weight)[block] += G_local.getNodeWeight(node) - G_local.getImplicitGhostNodes(node);
    }
}

void graph_io_stream::countAssignedNodes(PartitionConfig &config) {
    config.stream_assigned_nodes = 0;
    for (int i = 0; i < (int) config.stream_blocks_weight->size(); i++) {
        config.stream_assigned_nodes += (*config.stream_blocks_weight)[i];
    }
}

// This method modifies the matrix edges_virtualReal 
void
graph_io_stream::onePassPartition(PartitionConfig &config, std::vector <std::vector<EdgeWeight>> &edges_virtualReal,
                                  std::vector <PartitionID> &blockVirtualToReal,
                                  std::vector <NodeWeight> &weight_VirtualBlocks) {
    int len = edges_virtualReal.size();
    for (int i = 0; i < len; i++) {
        PartitionID decision = onePassDecide(config, i, edges_virtualReal[i]);
        blockVirtualToReal[i] = decision;
        (*config.stream_blocks_weight)[decision] += weight_VirtualBlocks[i];
    }
}

int graph_io_stream::onePassDecide(PartitionConfig &config, NodeID node, std::vector <EdgeWeight> &edges_i_real) {
    PartitionID decision;
    double best = std::numeric_limits<double>::lowest();
    LongNodeID global_node = (LongNodeID)
    node + config.lower_global_node - 1;
    double score = 0;
    double fennel_weight = getFennelWeight(config);
    int blocks = config.k;
    EdgeWeight block_weight;
    switch (config.one_pass_algorithm) {
        case ONEPASS_HASHING:
            decision = fnv1a(global_node) % config.k;
            break;
        case ONEPASS_GREEDY:
            for (int j = 0; j < blocks; j++) {
                if ((double) edges_i_real[j] > best) {
                    decision = j;
                    best = (double) edges_i_real[j];
                }
            }
            break;
        case ONEPASS_LDG:
            for (int j = 0; j < blocks; j++) {
                block_weight = (*config.stream_blocks_weight)[j] + (*config.add_blocks_weight)[j];
                if (block_weight >= (EdgeWeight) config.stream_total_upperbound) {
                    continue;
                }
                score = (0.1 + edges_i_real[j]) * (config.stream_total_upperbound - block_weight) /
                        (double) config.stream_total_upperbound;
                if (score > best) {
                    decision = j;
                    best = score;
                }
            }
            break;
        case ONEPASS_FENNEL:
            for (int j = 0; j < blocks; j++) {
                block_weight = (*config.stream_blocks_weight)[j] + (*config.add_blocks_weight)[j];
                if (block_weight >= (EdgeWeight) config.stream_total_upperbound) {
                    continue;
                }
                score = (0.1 + edges_i_real[j]) -
                        fennel_weight * (config.fennel_alpha_gamma * std::pow(block_weight, config.fennel_gamma - 1));
                if (score > best) {
                    decision = j;
                    best = score;
                }
            }
            break;
        case ONEPASS_CHUNKING:
            decision = global_node % config.k;
            break;
        case ONEPASS_FRACTIONAL_GREEDY:
            for (int j = 0; j < blocks; j++) {
                score = (0.1 + edges_i_real[j]) - (double) config.stream_total_upperbound /
                                                  (config.stream_total_upperbound - (*config.stream_blocks_weight)[j] -
                                                   (*config.add_blocks_weight)[j]);
                if (score > best) {
                    decision = j;
                    best = score;
                }
            }
            break;
    }
    return decision;
}


double graph_io_stream::getFennelWeight(PartitionConfig &partition_config) {
    double fennel_weight = 2;
    double fennel_tmp = 0;
    switch (partition_config.fennel_dynamics) {
        case FENNELADP_ORIGINAL:
            fennel_weight = 1;
            break;
        case FENNELADP_DOUBLE:
            fennel_weight = 2;
            break;
        case FENNELADP_LINEAR:
            fennel_weight =
                    2 * partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            break;
        case FENNELADP_MID_LINEAR:
            fennel_tmp =
                    2 * partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            if (fennel_tmp <= 1) {
                fennel_weight = 2 * (fennel_tmp);
            }
            break;
        case FENNELADP_QUADRATIC:
            fennel_tmp =
                    partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            fennel_weight = 2 * fennel_tmp * fennel_tmp;
            break;
        case FENNELADP_MID_QUADRATIC:
            fennel_tmp =
                    2 * partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            if (fennel_tmp <= 1) {
                fennel_weight = 2 * fennel_tmp * fennel_tmp;
            }
            break;
        case FENNELADP_MID_CONSTANT:
            fennel_tmp =
                    partition_config.remaining_stream_nodes / (double) partition_config.total_stream_edges;
            if (fennel_tmp <= 1.5) {
                fennel_weight = 0.5;
            }
            break;
        case FENNELADP_EDGE_CUT:
            fennel_weight = 0;
            break;
    }
    return fennel_weight;
}


void graph_io_stream::writePartitionStream(PartitionConfig &config) {
    std::ofstream f(config.filename_output.c_str());
    std::cout << "writing partition to " << config.filename_output << " ... " << std::endl;

    for (LongNodeID node = 0; node < config.stream_nodes_assign->size(); node++) {
        f << (*config.stream_nodes_assign)[node] << "\n";
    }

    f.close();
}


void graph_io_stream::readFirstLineStream(PartitionConfig &partition_config, std::string graph_filename,
                                          EdgeWeight &total_edge_cut) {
    if (partition_config.stream_in != NULL) {
        delete partition_config.stream_in;
    }
    partition_config.stream_in = new std::ifstream(graph_filename.c_str());
    if (!(*(partition_config.stream_in))) {
        std::cerr << "Error opening " << graph_filename << std::endl;
        exit(1);
    }
    std::vector <std::string> *lines;

    lines = new std::vector<std::string>(1);
    std::getline(*(partition_config.stream_in), (*lines)[0]);

    //skip comments
    while ((*lines)[0][0] == '%') {
        std::getline(*(partition_config.stream_in), (*lines)[0]);
    }

    std::stringstream ss((*lines)[0]);
    ss >> partition_config.remaining_stream_nodes;
    ss >> partition_config.remaining_stream_edges;
    ss >> partition_config.remaining_stream_ew;

    if (partition_config.stream_nodes_assign == NULL) {
        partition_config.stream_nodes_assign = new std::vector<PartitionID>(partition_config.remaining_stream_nodes,
                                                                            INVALID_PARTITION);
    }
    if (partition_config.stream_blocks_weight == NULL) {
        partition_config.stream_blocks_weight = new std::vector<NodeWeight>(partition_config.k, 0);
    }
    if (partition_config.add_blocks_weight == NULL) {
        partition_config.add_blocks_weight = new std::vector<NodeWeight>(partition_config.k, 0);
    }
    partition_config.total_stream_nodeweight = 0;
    partition_config.total_stream_nodecounter = 0;
    partition_config.stream_n_nodes = partition_config.remaining_stream_nodes;
    partition_config.total_stream_edges = partition_config.remaining_stream_nodes;

    auto total_weight = (partition_config.balance_edges) ? (partition_config.remaining_stream_nodes +
                                                            2 * partition_config.remaining_stream_edges)
                                                         : partition_config.remaining_stream_nodes;

    if (partition_config.num_streams_passes > 1 + partition_config.restream_number) {
        partition_config.stream_total_upperbound = ceil(
                ((100 + 1.5 * partition_config.imbalance) / 100.) * (total_weight / (double) partition_config.k));
    } else {
        partition_config.stream_total_upperbound = ceil(
                ((100 + partition_config.imbalance) / 100.) * (total_weight / (double) partition_config.k));
    }

    partition_config.fennel_alpha =
            partition_config.remaining_stream_edges * std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
            (std::pow(partition_config.remaining_stream_nodes, partition_config.fennel_gamma));

    partition_config.fennel_alpha_gamma = partition_config.fennel_alpha * partition_config.fennel_gamma;

    partition_config.quotient_nodes = partition_config.k;

    total_edge_cut = 0;
    if (partition_config.stream_buffer_len == 0) { // signal of partial restream standard buffer size
        partition_config.stream_buffer_len = (LongNodeID)
        ceil(partition_config.remaining_stream_nodes / (double) partition_config.k);
    }
    partition_config.nmbNodes = MIN(partition_config.stream_buffer_len, partition_config.remaining_stream_nodes);
    partition_config.n_batches = ceil(partition_config.remaining_stream_nodes / (double) partition_config.nmbNodes);
    partition_config.curr_batch = 0;
//	partition_config.stream_global_epsilon = (partition_config.imbalance)/100.;

    delete lines;
}


void graph_io_stream::loadRemainingLines(PartitionConfig &partition_config, LINE_BUFFER &lines) {
    if (partition_config.ram_stream) {
        lines = graph_io_stream::loadLinesFromStream(partition_config, partition_config.remaining_stream_nodes);
    }
}


void graph_io_stream::loadBufferLines(PartitionConfig &partition_config, LINE_BUFFER &lines, LongNodeID num_lines) {
    if (!partition_config.ram_stream) {
        lines = graph_io_stream::loadLinesFromStream(partition_config, num_lines);
    }
}


std::vector <std::string> *
graph_io_stream::loadLinesFromStream(PartitionConfig &partition_config, LongNodeID num_lines) {
    std::vector <std::string> *lines;
    lines = new std::vector<std::string>(num_lines);
    LongNodeID node_counter = 0;
    while (node_counter < num_lines) {
        std::getline(*(partition_config.stream_in), (*lines)[node_counter]);
        if ((*lines)[node_counter][0] == '%') { // a comment in the file
            continue;
        }
        node_counter++;
    }
    return lines;
}

void graph_io_stream::prescribeBufferInbalance(PartitionConfig &partition_config) {
    double &global_epsilon = partition_config.stream_global_epsilon;
    int &passes = partition_config.num_streams_passes;
    if (partition_config.restream_number &&
        partition_config.restream_number >= partition_config.num_streams_passes - 1) {
        partition_config.imbalance = 100 * global_epsilon;
    } else {
        double current_nodes = (double) partition_config.stream_assigned_nodes + partition_config.nmbNodes +
                               partition_config.ghost_nodes;
        partition_config.imbalance = 100 * partition_config.stream_n_nodes * (1 + global_epsilon) / current_nodes - 100;
        if (passes > 1) {
            partition_config.imbalance = MIN(partition_config.imbalance, partition_config.batch_inbalance);
        } else {
            partition_config.imbalance = MIN(MAX(100 * global_epsilon, 0.75 * partition_config.imbalance),
                                             partition_config.batch_inbalance);
        }
    }
}

void
graph_io_stream::streamEvaluatePartition(PartitionConfig &config, const std::string &filename, EdgeWeight &edgeCut) {
    std::vector <std::vector<LongNodeID>> *input;
    std::vector <std::string> *lines;
    lines = new std::vector<std::string>(1);
    LongNodeID node_counter = 0;
    buffered_input *ss2 = NULL;
    std::string line;
    std::ifstream in(filename.c_str());
    if (!in) {
        std::cerr << "Error opening " << filename << std::endl;
        return 1;
    }
    long nmbNodes;
    long nmbEdges;
    int ew = 0;
    std::getline(in, (*lines)[0]);
    while ((*lines)[0][0] == '%') std::getline(in, (*lines)[0]); // a comment in the file

    std::stringstream ss((*lines)[0]);
    ss >> nmbNodes;
    ss >> nmbEdges;
    ss >> ew;
    bool read_ew = false;
    bool read_nw = false;
    if (ew == 1) {
        read_ew = true;
    } else if (ew == 11) {
        read_ew = true;
        read_nw = true;
    } else if (ew == 10) {
        read_nw = true;
    }
    NodeID target;
    NodeWeight total_nodeweight = 0;
    EdgeWeight total_edgeweight = 0;
    edgeCut = 0;

    while (std::getline(in, (*lines)[0])) {
        if ((*lines)[0][0] == '%') continue; // a comment in the file
        NodeID node = node_counter++;
        PartitionID partitionIDSource = (*config.stream_nodes_assign)[node];

        input = new std::vector <std::vector<LongNodeID>>(1);
        ss2 = new buffered_input(lines);
        ss2->simple_scan_line((*input)[0]);
        std::vector <LongNodeID> &line_numbers = (*input)[0];
        LongNodeID col_counter = 0;

        NodeWeight weight = 1;
        if (read_nw) {
            weight = line_numbers[col_counter++];
            total_nodeweight += weight;
        }
        while (col_counter < line_numbers.size()) {
            target = line_numbers[col_counter++];
            target = target - 1;
            EdgeWeight edge_weight = 1;
            if (read_ew) {
                edge_weight = line_numbers[col_counter++];
            }
            total_edgeweight += edge_weight;
            PartitionID partitionIDTarget = (*config.stream_nodes_assign)[target];
            if (partitionIDSource != partitionIDTarget) {
                edgeCut += edge_weight;
            }
        }
        (*lines)[0].clear();
        delete ss2;
        delete input;
        if (in.eof()) {
            break;
        }
    }
    edgeCut = edgeCut / 2; // Since every edge is counted twice
    delete lines;
}


template<typename T>
T graph_io_stream::return_and_delete_element(std::vector <T> &vec, LongNodeID pos) {
    T val = vec[pos];
    vec[pos] = vec[vec.size() - 1];
    vec.erase(vec.begin() + vec.size() - 1);
    return val;
}

// for edge partitioning
void graph_io_stream::readFirstLineStreamEdge(PartitionConfig &partition_config, std::string graph_filename,
                                              EdgeWeight &total_edge_cut) {

    std::string bin_ending(".bin");
    std::string parhip_ending(".parhip");
    if (hasEnding(graph_filename, bin_ending) || hasEnding(graph_filename, parhip_ending)) {
        std::vector<unsigned long long> buffer(3, 0);
        partition_config.stream_in = new std::ifstream(graph_filename.c_str(), std::ios::binary | std::ios::in);;
        if ((*partition_config.stream_in)) {
            (*partition_config.stream_in).read((char *) (&buffer[0]), 3 * sizeof(unsigned long long));
        }

        unsigned long long version = buffer[0];
        partition_config.remaining_stream_nodes = static_cast<NodeID>(buffer[1]);
        partition_config.remaining_stream_edges = static_cast<NodeID>(buffer[2]) / 2;
    } else {
        if (partition_config.stream_in != NULL) {
            delete partition_config.stream_in;
        }
        partition_config.stream_in = new std::ifstream(graph_filename.c_str());
        if (!(*(partition_config.stream_in))) {
            std::cerr << "Error opening " << graph_filename << std::endl;
            exit(1);
        }
        std::vector <std::string> *lines;

        lines = new std::vector<std::string>(1);
        std::getline(*(partition_config.stream_in), (*lines)[0]);

        // skip comments
        while ((*lines)[0][0] == '%') {
            std::getline(*(partition_config.stream_in), (*lines)[0]);
        }

        std::stringstream ss((*lines)[0]);
        ss >> partition_config.remaining_stream_nodes;
        ss >> partition_config.remaining_stream_edges;
        ss >> partition_config.remaining_stream_ew;

        delete lines;
    }

    if (!partition_config.filename_output.compare("")) {
        partition_config.filename_output = "tmp_output.txt";
    }

    if (partition_config.stream_output_progress) {
        partition_config.stream_out = new std::ofstream(partition_config.filename_output.c_str());
    }

    // storing number of nodes of input graph
    partition_config.remaining_stream_nodes_OG = partition_config.remaining_stream_nodes;
    partition_config.remaining_stream_graph_nodes = partition_config.remaining_stream_nodes;
    // setting num. nodes of input graph = num. edges as necessitated by model construction
    partition_config.remaining_stream_nodes = partition_config.remaining_stream_edges;
    partition_config.total_stream_edges = partition_config.remaining_stream_nodes;
    std::cout << "Input has " << partition_config.remaining_stream_nodes_OG << " nodes " << " and "  << partition_config.remaining_stream_nodes << " edges" << std::endl;

    if (partition_config.stream_nodes_assign == NULL && !partition_config.benchmark &&
        !partition_config.stream_output_progress) {
        partition_config.stream_nodes_assign = new std::vector<PartitionID>(partition_config.remaining_stream_nodes,
                                                                            INVALID_PARTITION);
    }

    // when not using minimal mode, store blocks incident on nodes of input graph in a google dense hash set
    if (partition_config.blocks_on_node == NULL && !partition_config.minimal_mode) {
        partition_config.blocks_on_node = new std::vector <google::dense_hash_set<PartitionID>>(
                partition_config.remaining_stream_nodes_OG, google::dense_hash_set<PartitionID>());
    }

    // when using minimal mode, store blocks incident on nodes of input graph in 1D vector
    if (partition_config.blocks_on_node_minimal == NULL && partition_config.minimal_mode) {
        std::cout << "Minimal mode" << std::endl;
        partition_config.blocks_on_node_minimal = new std::vector<NodeID>(partition_config.remaining_stream_nodes_OG,
                                                                          -1);
    }

    if (partition_config.stream_blocks_weight == NULL) {
        partition_config.stream_blocks_weight = new std::vector<NodeWeight>(partition_config.k, 0);
    }
    if (partition_config.add_blocks_weight == NULL) {
        partition_config.add_blocks_weight = new std::vector<NodeWeight>(partition_config.k, 0);
    }
    partition_config.total_stream_nodeweight = 0;
    partition_config.total_stream_nodecounter = 0;
    partition_config.stream_n_nodes = partition_config.remaining_stream_nodes;

    if (partition_config.num_streams_passes > 1 + partition_config.restream_number) {
        partition_config.stream_total_upperbound = ceil(((100 + 1.5 * partition_config.imbalance) / 100.) *
                                                        (partition_config.remaining_stream_nodes /
                                                         (double) partition_config.k));
    } else {
        partition_config.stream_total_upperbound = ceil(((100 + partition_config.imbalance) / 100.) *
                                                        (partition_config.remaining_stream_nodes /
                                                         (double) partition_config.k));
    }

    // configurations for dynamic / static / batch Fennel Alpha
    if (partition_config.dynamic_alpha) {
        partition_config.fennel_edges = 6 * partition_config.remaining_stream_edges;
    } else {
        partition_config.fennel_edges = 5 * partition_config.remaining_stream_edges;
    }
    if (partition_config.num_split_edges == -1) {
        partition_config.fennel_alpha =
                partition_config.fennel_edges * std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
                (std::pow(partition_config.remaining_stream_nodes, partition_config.fennel_gamma));
    } else {
        partition_config.fennel_alpha = partition_config.remaining_stream_edges *
                                        std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
                                        (std::pow(partition_config.remaining_stream_nodes_OG,
                                                  partition_config.fennel_gamma));
    }

    partition_config.fennel_alpha_gamma = partition_config.fennel_alpha * partition_config.fennel_gamma;

    partition_config.quotient_nodes = partition_config.k;

    total_edge_cut = 0;
    if (partition_config.stream_buffer_len == 0) { // signal of partial restream standard buffer size
        partition_config.stream_buffer_len = (LongNodeID)
        ceil(partition_config.remaining_stream_nodes / (double) partition_config.k);
    }
    partition_config.nmbNodes = MIN(partition_config.stream_buffer_len, partition_config.remaining_stream_nodes);
    partition_config.n_batches = ceil(partition_config.remaining_stream_nodes / (double) partition_config.nmbNodes);
    partition_config.curr_batch = 0;

}

void graph_io_stream::generalizeStreamPartitionEdge(PartitionConfig &config, graph_access &G_local) {
    for (NodeID node = 0, end = config.nmbNodes; node < end; node++) {
        PartitionID block = G_local.getPartitionIndex(node);
        // Here, global_node controls the mapping of the node in the graph model
        // G_local to the node in the original graph G
        LongNodeID global_node = (LongNodeID)
        node + config.lower_global_node - 1;

        // if maintaining vector of partition IDs for all edges
        if (!config.benchmark && !config.stream_output_progress) {
            (*config.stream_nodes_assign)[global_node] = block;
        }
        (*config.stream_blocks_weight)[block] += G_local.getNodeWeight(node) - G_local.getImplicitGhostNodes(node);
        // assign block to incident nodes of the edge that "node" corresponds to
        timer ppt;
        ppt.restart();
        //  incident_node_list contains the original nodes u_OG and v_OG from the
        //  original graph that are incident on the edge (node) we are currently
        //  looking at
        std::vector <NodeID> incident_node_list;
        incident_node_list = (*config.nodes_on_edge_conv)[node];

        if (config.past_subset_size != 0) {
            for (auto &incident_node: incident_node_list) {

                if (config.minimal_mode == false) {
                    if ((*config.blocks_on_node)[incident_node].empty() == true) {
                        (*config.blocks_on_node)[incident_node].set_empty_key(-1);
                        (*config.blocks_on_node)[incident_node].insert(block);
                        continue;
                    }
                    // if blocks are already there
                    if ((*config.blocks_on_node)[incident_node].empty() == false) {
                        // add block only if not in unordered set
                        if ((*config.blocks_on_node)[incident_node].find(block) ==
                            (*config.blocks_on_node)[incident_node].end()) {
                            (*config.blocks_on_node)[incident_node].insert(block);
                        }
                    }
                } else {
                    (*config.blocks_on_node_minimal)[incident_node] = block;
                }
            }
        }
        config.finding_past_assignments_time += ppt.elapsed();
    }
    //f.close();
}

bool graph_io_stream::hasEnding(std::string const &string, std::string const &ending) {
    if (string.length() >= ending.length()) {
        return (0 == string.compare(string.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void graph_io_stream::find_reverse_edges(graph_access &G, std::vector <EdgeID> &m_reverse_edge) {
    // reverse edges were already calculated, thus nothing to do
    if (!m_reverse_edge.empty()) {
        return;
    }

    const EdgeID guard = std::numeric_limits<EdgeID>::max();
    m_reverse_edge.resize(G.number_of_edges(), guard);

    for (NodeID u = 0; u < G.number_of_nodes(); ++u) {
        for (EdgeID e_uv = G.get_first_edge(u); e_uv < G.get_first_invalid_edge(u); ++e_uv) {
            NodeID v = G.getEdgeTarget(e_uv);
            if (u < v) {
                assert(m_reverse_edge[e_uv] == guard);
            } else {
                assert(m_reverse_edge[e_uv] != guard);
                continue;
            }

            for (EdgeID e_vu = G.get_first_edge(v); e_vu < G.get_first_invalid_edge(v); ++e_vu) {
                if (G.getEdgeTarget(e_vu) == u) {
                    m_reverse_edge[e_uv] = e_vu;
                    m_reverse_edge[e_vu] = e_uv;

                    break;
                }
            }

            assert(m_reverse_edge[e_uv] != guard);
        }
    }
}

EdgeID graph_io_stream::insertRegularEdgeInBatchEdge(PartitionConfig &config,
                                                     std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

> &all_edges,
NodeID node, LongNodeID
global_target,
EdgeWeight edge_weight
) {
NodeID target = (NodeID) (global_target - config.lower_global_node);
if (target == node) {
std::cerr << "The graph file contains self-loops, which are not supported. "
"Please remove them from the file."
<<
std::endl;
exit(0);
}
if (target > node) {
return 0;
}

return
includeEdgeInBatchEdge(all_edges, node, target,
(1 + config.double_non_ghost_edges) *
edge_weight);
}

EdgeID graph_io_stream::includeEdgeInBatchEdge(std::vector < std::vector < std::pair < NodeID, EdgeWeight
>>> &all_edges,
NodeID node, NodeID
target,
EdgeWeight edge_weight
) {
// std::cout << node << " -> " << target << ": " << edge_weight << std::endl;
all_edges[node].
push_back(std::make_pair(target, edge_weight)
);
all_edges[target].
push_back(std::make_pair(node, edge_weight)
);

return 2;
}

void
graph_io_stream::processQuotientEdgeInBatchEdge(PartitionConfig &config, NodeID node, LongNodeID global_target,
                                                EdgeWeight edge_weight) {
    PartitionID targetGlobalPar = global_target;
    if ((*config.edge_block_nodes)[targetGlobalPar].size() >= 1) {
        auto &curr_element = (*config.edge_block_nodes)[targetGlobalPar][0];
        // this scenario did not occur in my testing
        if (curr_element.first == node) {
            curr_element.second += edge_weight;
            return;
        }
    }
    std::pair <NodeID, NodeWeight> element = {node, edge_weight};
    (*config.edge_block_nodes)[targetGlobalPar].push_back(element);
}

void
graph_io_stream::insertQuotientNodesInBatchEdge(PartitionConfig &config, std::vector <NodeWeight> &all_nodes,
                                                NodeID &node_counter) {
    while (node_counter < config.nmbNodes + config.quotient_nodes) { // Add artificial nodes representing
        // current blocks
        NodeID node = node_counter++;
        NodeID targetPar = node - (NodeID)
        config.nmbNodes;
        NodeWeight weight = (*config.stream_blocks_weight)[targetPar];
        all_nodes[node] = weight;
    }
}

EdgeID graph_io_stream::insertQuotientEdgesInBatchEdge(PartitionConfig &config,
                                                       std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

> &all_edges) {
// std::cout << "Quotient edges: " << std::endl;
EdgeID inserted_edges = 0;
EdgeWeight edge_weight;
for (
PartitionID block = 0;
block<config.
k;
block++) {
NodeID target = config.nmbNodes + block;
if ((*config.edge_block_nodes)[block].

size()

< 1)
continue; // if block has no neighbors in batch, continue outer for loop
for (
auto &element
: (*config.edge_block_nodes)[block]) {
// if((1+config.double_non_ghost_edges)*element.second > 2)std::cout <<
// "Weight of quotient edge: " <<
// (1+config.double_non_ghost_edges)*element.second << std::endl;
inserted_edges +=
includeEdgeInBatchEdge(
        all_edges, element
.first, target,
(1 + config.double_non_ghost_edges) * element.second);
}}
return
inserted_edges;
}

void graph_io_stream::createGraphForBatchEdge(PartitionConfig &config, graph_access &G, NodeID node_counter,
                                              EdgeID edge_counter,
                                              std::vector <std::vector<std::pair < NodeID, EdgeWeight>>

> &all_edges,
std::vector <NodeWeight> &all_nodes
) {
G.
stream_repeat_construction(node_counter, edge_counter
);
G.
resizeImplicitGhostNodesList(node_counter);
std::vector <EdgeWeight> neighbor_edgeweight(node_counter, 0);
std::vector <EdgeID> neighbor_edgeid(node_counter, 0);
std::vector <NodeID> neighbors;
EdgeID e;
for (
NodeID i = 0;
i<node_counter;
i++) { // actual insertion of nodes and edges
NodeID node = G.new_node();
NodeWeight weight = all_nodes[node];
G.
setNodeWeight(node, weight
);
G.
setImplicitGhostNodes(node,
0);
G.
setPartitionIndex(node,
0);
recoverBlockAssignedToNode(config, G, node, node_counter
);
for (auto &[target, edge_weight]: all_edges[node]) {
if (neighbor_edgeweight[target] == 0) {
e = G.new_edge(node, target);
neighbor_edgeid[target] =
e;
neighbors.
push_back(target);
} else {
e = neighbor_edgeid[target];
}
neighbor_edgeweight[target] +=
edge_weight;
G.
setEdgeWeight(e, neighbor_edgeweight[target]
);
}
for (
auto &target
: neighbors) {
neighbor_edgeweight[target] = 0;
}
neighbors.

clear();

}
G.

stream_finish_construction();

}

void graph_io_stream::addEdgesToIncidentBlock(PartitionConfig &partition_config, NodeID node, NodeID target,
                                              EdgeWeight edge_weight, LongEdgeID &used_edges) {
    if (partition_config.minimal_mode) {
        if ((*partition_config.blocks_on_node_minimal)[(*partition_config.nodes_on_edge_conv)[node][0]] != -1) {
            //std::cout << "Here for node " << node << " with target " << target << std::endl;
            used_edges++;
            target = (*partition_config.blocks_on_node_minimal)[(*partition_config.nodes_on_edge_conv)[node][0]];
            EdgeWeight edge_weight = 1;
            partition_config.quotient_edges_count++;
            processQuotientEdgeInBatchEdge(partition_config, node, target, edge_weight);
        }
    } else {
        if (!(*partition_config.blocks_on_node)[(*partition_config.nodes_on_edge_conv)[node][0]].empty()) {
            if (partition_config.past_subset_size == -1 || (partition_config.past_subset_size != 0 &&
                                                            (*partition_config.blocks_on_node)[(*partition_config.nodes_on_edge_conv)[node][0]].size() <=
                                                            partition_config.past_subset_size)) {
                // if all blocks incident are allowed to be treated as quotient nodes or
                // if number of blocks is less than specified past subset size
                for (google::dense_hash_set<PartitionID>::const_iterator it = (*partition_config.blocks_on_node)[(*partition_config.nodes_on_edge_conv)[node][0]].begin();
                     it !=
                     (*partition_config.blocks_on_node)[(*partition_config.nodes_on_edge_conv)[node][0]].end(); ++it) {
                    used_edges++;
                    target = *it;
                    EdgeWeight edge_weight = 1;

                    partition_config.quotient_edges_count++;
                    processQuotientEdgeInBatchEdge(partition_config, node, target, edge_weight);
                }
            } else if (partition_config.past_subset_size != 0 && partition_config.past_subset_size != -1) {
                // otherwise, random subset of blocks to be considered as quotient nodes
                std::vector <PartitionID> in;
                int curr_idx = 0;
                for (google::dense_hash_set<PartitionID>::const_iterator it = (*partition_config.blocks_on_node)[(*partition_config.nodes_on_edge_conv)[node][0]].begin();
                     it !=
                     (*partition_config.blocks_on_node)[(*partition_config.nodes_on_edge_conv)[node][0]].end(); ++it) {
                    in.push_back(*it);
                }
                std::vector<unsigned int> reservoir;
                reservoir.reserve(partition_config.past_subset_size);
                std::random_device rd;
                std::mt19937_64 gen(rd());
                for (size_t i = 0; i < in.size(); ++i) {
                    if (reservoir.size() < partition_config.past_subset_size) {
                        reservoir.push_back(in[i]);
                    } else {
                        std::uniform_int_distribution <size_t> dis(0, i);
                        size_t j = dis(gen);
                        if (j < partition_config.past_subset_size) {
                            reservoir[j] = in[i];
                        }
                    }
                }
                for (auto i: reservoir) {
                    used_edges++;
                    target = i;
                    EdgeWeight edge_weight = 1;

                    partition_config.quotient_edges_count++;
                    processQuotientEdgeInBatchEdge(partition_config, node, target, edge_weight);
                }
            } else {
                // for = 0, do nothing
            }
        }
    }
}

void graph_io_stream::streamEvaluatePartitionEdgeBatch(PartitionConfig &config, const std::string &filename,
                                                       NodeID &vertex_cut, NodeID &replicas,
                                                       double &replication_factor, double &balance) {

    // load entire graph from input text file
    graph_access G_temp;
    long nmbNodes;
    long nmbEdges;
    if (hasEnding(config.graph_filename, ".bin") || hasEnding(config.graph_filename, ".parhip")) {
        std::vector<unsigned long long> buffer(3, 0);
        std::ifstream filebin(config.graph_filename.c_str(), std::ios::binary | std::ios::in);
        if (filebin) {
            filebin.read((char *) (&buffer[0]), 3 * sizeof(unsigned long long));
        }

        unsigned long long version = buffer[0];
        nmbNodes = static_cast<NodeID>(buffer[1]);
        nmbEdges = static_cast<NodeID>(buffer[2]) / 2;

        unsigned long long nodes_in_batch = static_cast<unsigned long long>(nmbNodes);
        unsigned long long *vertex_offsets = new unsigned long long[nodes_in_batch + 1];
        filebin.seekg(3 * (sizeof(unsigned long long)));
        filebin.read((char *) (vertex_offsets), (nodes_in_batch + 1) * sizeof(unsigned long long));

        unsigned long long edge_start_pos = vertex_offsets[0];
        unsigned long long num_reads = vertex_offsets[nodes_in_batch] - vertex_offsets[0];
        unsigned long long num_edges_to_read = num_reads / sizeof(unsigned long long);
        unsigned long long *edges = new unsigned long long[num_edges_to_read]; // we also need the next vertex offset
        filebin.seekg(edge_start_pos);
        filebin.read((char *) (edges), (num_edges_to_read) * sizeof(unsigned long long));

        NodeID node;
        G_temp.start_construction(nmbNodes, nmbEdges * 2);
        unsigned long long pos = 0;
        for (unsigned long long i = 0; i < nodes_in_batch; ++i) {
            node = static_cast<NodeID>(i);
            G_temp.new_node();
            G_temp.setPartitionIndex(node, 0);
            G_temp.setNodeWeight(node, 1);
            unsigned long long degree = (vertex_offsets[i + 1] - vertex_offsets[i]) / sizeof(unsigned long long);
            for (unsigned long long j = 0; j < degree; j++, pos++) {
                auto target = static_cast<NodeID>(edges[pos]);
                const EdgeID e = G_temp.new_edge(node, target);
                G_temp.setEdgeWeight(e, 1);
            }
        }
        G_temp.finish_construction();
        delete[] vertex_offsets;
        delete[] edges;
        filebin.close();
    } else {
        std::vector <std::vector<LongNodeID>> *input;
        std::vector <std::string> *lines;
        lines = new std::vector<std::string>(1);
        LongNodeID node_counter = 0;
        buffered_input *ss2 = NULL;
        std::string line;
        std::ifstream in(filename.c_str());
        if (!in) {
            std::cerr << "Error opening " << filename << std::endl;
            // return 1;
        }
        int ew = 0;
        std::getline(in, (*lines)[0]);
        while ((*lines)[0][0] == '%')
            std::getline(in, (*lines)[0]); // a comment in the file

        std::stringstream ss((*lines)[0]);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;
        bool read_ew = false;
        bool read_nw = false;
        if (ew == 1) {
            read_ew = true;
        } else if (ew == 11) {
            read_ew = true;
            read_nw = true;
        } else if (ew == 10) {
            read_nw = true;
        }
        NodeID target;
        NodeWeight total_nodeweight = 0;
        EdgeWeight total_edgeweight = 0;
        EdgeID edge = 0;

        G_temp.start_construction(nmbNodes, nmbEdges * 2);

        while (std::getline(in, (*lines)[0])) {
            if ((*lines)[0][0] == '%')
                continue; // a comment in the file
            NodeID node = node_counter++;
            G_temp.new_node();
            G_temp.setPartitionIndex(node, 0);
            G_temp.setNodeWeight(node, 1);
            input = new std::vector <std::vector<LongNodeID>>(1);
            ss2 = new buffered_input(lines);
            ss2->simple_scan_line((*input)[0]);
            std::vector <LongNodeID> &line_numbers = (*input)[0];
            LongNodeID col_counter = 0;

            NodeWeight weight = 1;
            if (read_nw) {
                weight = line_numbers[col_counter++];
                total_nodeweight += weight;
            }
            while (col_counter < line_numbers.size()) {
                target = line_numbers[col_counter++];
                target = target - 1;
                EdgeWeight edge_weight = 1;
                if (read_ew) {
                    edge_weight = line_numbers[col_counter++];
                }
                total_edgeweight += edge_weight;
                const EdgeID e = G_temp.new_edge(node, target);
                G_temp.setEdgeWeight(e, 1);
            }
            (*lines)[0].clear();
            delete ss2;
            delete input;
            if (in.eof()) {
                break;
            }
        }
        G_temp.finish_construction();
        delete lines;
    }

    // load partition
    if (config.evaluate_mode) {
        if (config.stream_nodes_assign != NULL) {
            delete config.stream_nodes_assign;
            config.stream_nodes_assign = new std::vector<PartitionID>(nmbEdges, INVALID_PARTITION);
        } else {
            config.stream_nodes_assign = new std::vector<PartitionID>(nmbEdges, INVALID_PARTITION);
        }
        std::string line;
        std::ifstream part_file(config.filename_output);
        if (!part_file) {
            std::cerr << "Error opening partition ID file." << filename << std::endl;
            exit;
        }
        for (int i = 0; i < nmbEdges; i++) {
            // fetch current line
            std::getline(part_file, line);
            if (line[0] == '%') { // Comment
                continue;
            }
            (*config.stream_nodes_assign)[i] = (PartitionID) atol(line.c_str());
        }
        part_file.close();
    }

    std::vector <EdgeID> m_reverse_edge;
    if (!config.evaluate_mode)
        graph_io_stream::find_reverse_edges(G_temp, m_reverse_edge);
    else {
        m_reverse_edge.resize(G_temp.number_of_edges());
        std::vector <std::vector<EdgeID>> temp_rev_storage(G_temp.number_of_nodes());
        NodeID back_idx = 0;
        for (NodeID u = 0; u < G_temp.number_of_nodes(); ++u) {
            back_idx = 0;
            for (EdgeID e = G_temp.get_first_edge(u); e < G_temp.get_first_invalid_edge(u); ++e) {
                NodeID j = G_temp.getEdgeTarget(e);
                if (u < j) { // if left to right edge, push edge ID to index j
                    temp_rev_storage[j].push_back(e);
                } else { // if right to left edge, fetch edge ID from index u (which was
                    // j earlier) of corresponding left to right edge
                    m_reverse_edge[temp_rev_storage[u][back_idx]] = e;
                    m_reverse_edge[e] = temp_rev_storage[u][back_idx];
                    back_idx++;
                }
            }
        }
    }

    std::vector <PartitionID> edge_partition(nmbEdges * 2, -1);
    NodeID edgeCount = 0;
    NodeID remaining_nodes = nmbNodes;
    NodeID batchSize = MIN(config.stream_buffer_len, remaining_nodes);
    NodeID lower_node = 0;
    NodeID upper_node = lower_node + batchSize;
    NodeID incremental_edge_ID = 0;

    while (remaining_nodes != 0) {
        for (NodeID u = lower_node; u < upper_node; u++) {
            for (EdgeID e = G_temp.get_first_edge(u); e < G_temp.get_first_invalid_edge(u); ++e) {
                NodeID v = G_temp.getEdgeTarget(e);
                if (v >= upper_node)
                    continue;
                if (v < lower_node) {
                    edge_partition[e] = (*config.stream_nodes_assign)[edgeCount];
                    edge_partition[m_reverse_edge[e]] = (*config.stream_nodes_assign)[edgeCount];
                    edgeCount++;

                } else {
                    if (u < v) {
                        edge_partition[e] = (*config.stream_nodes_assign)[edgeCount];
                        edge_partition[m_reverse_edge[e]] = (*config.stream_nodes_assign)[edgeCount];
                        edgeCount++;
                    }
                }
            }
        }
        lower_node = upper_node;
        remaining_nodes = remaining_nodes - batchSize;
        batchSize = MIN(config.stream_buffer_len, remaining_nodes);
        upper_node = lower_node + batchSize;
    }

    std::vector <NodeID> v_ei(config.k, 0);
    std::vector <NodeID> block_weights(config.k, 0);

    for (NodeID u = 0; u < G_temp.number_of_nodes(); ++u) {
        if (G_temp.getNodeDegree(u) == 0)
            continue;

        std::vector<bool> counted(config.k);
        NodeID incident_part_count = 0;

        for (EdgeID e = G_temp.get_first_edge(u); e < G_temp.get_first_invalid_edge(u); ++e) {
            PartitionID part = edge_partition[e];
            if (u < G_temp.getEdgeTarget(e))
                block_weights[part]++;
            if (!counted[part]) {
                counted[part] = true;
                v_ei[part]++;
                ++incident_part_count;
                ++replicas;
            }
        }
        --replicas;
        if (incident_part_count > 1)
            vertex_cut++;
    }

    double sum = 0;
    double max = -1;
    double total_weight = 0;
    for (int i = 0; i < v_ei.size(); i++) {
        sum += v_ei[i];
        if (block_weights[i] > max) {
            max = block_weights[i];
        }
        if (block_weights[i] > config.stream_total_upperbound) {
            std::cout << "Partition is imbalanced" << std::endl;
        }
        total_weight += block_weights[i];
    }
    replication_factor = sum / (double) G_temp.number_of_nodes();
    double balance_part_weight = ceil(total_weight / (double) config.k);
    balance = max / balance_part_weight;
}

void
graph_io_stream::streamEvaluateEdgePartition(PartitionConfig &config, const std::string &filename, NodeID &vertex_cut,
                                             NodeID &replicas, double &replication_factor, double &balance) {
    // evaluate using (n*k) vector instead of loading entire graph
    NodeID nmbNodes;
    EdgeID nmbEdges;
    std::string bin_ending(".bin");
    std::string parhip_ending(".parhip");
    std::ifstream *in;
    if (hasEnding(config.graph_filename, bin_ending) || hasEnding(config.graph_filename, parhip_ending)) {
        std::vector<unsigned long long> buffer(3, 0);
        in = new std::ifstream(config.graph_filename.c_str(), std::ios::binary | std::ios::in);;
        if ((*in)) {
            (*in).read((char *) (&buffer[0]), 3 * sizeof(unsigned long long));
        }
        unsigned long long version = buffer[0];
        nmbNodes = static_cast<NodeID>(buffer[1]);
        nmbEdges = static_cast<NodeID>(buffer[2]) / 2;

    } else {
        in = new std::ifstream(config.graph_filename.c_str());
        if (!(*(in))) {
            std::cerr << "Error opening " << config.graph_filename << std::endl;
            exit(1);
        }
        std::vector <std::string> *lines;
        int ew = 0;

        lines = new std::vector<std::string>(1);
        std::getline(*(in), (*lines)[0]);

        // skip comments
        while ((*lines)[0][0] == '%') {
            std::getline(*(in), (*lines)[0]);
        }

        std::stringstream ss((*lines)[0]);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;

        delete lines;
    }

    // load partition
    if (config.evaluate_mode) {
        if (config.stream_nodes_assign != NULL) {
            delete config.stream_nodes_assign;
            config.stream_nodes_assign = new std::vector<PartitionID>(config.total_stream_edges, INVALID_PARTITION);
        } else {
            config.stream_nodes_assign = new std::vector<PartitionID>(config.total_stream_edges, INVALID_PARTITION);
        }
        std::string line;
        std::ifstream part_file(config.filename_output);
        if (!part_file) {
            std::cerr << "Error opening partition ID file." << filename << std::endl;
            exit;
        }
        for (int i = 0; i < nmbEdges; i++) {
            // fetch current line
            std::getline(part_file, line);
            if (line[0] == '%') { // Comment
                continue;
            }
            (*config.stream_nodes_assign)[i] = (PartitionID) atol(line.c_str());
        }
        part_file.close();
    }

    NodeID remaining_nodes = nmbNodes;
    NodeID batch_size = MIN(config.stream_buffer_len, remaining_nodes);
    NodeID lower_node = 0;
    NodeID upper_node = lower_node + batch_size;
    NodeID node, target;
    EdgeID incremental_edge_counter = 0;
    unsigned long long start_pos = 3 * sizeof(unsigned long long);
    std::vector <std::vector<PartitionID>> part_ids_on_nodes(nmbNodes, std::vector<PartitionID>(config.k, 0));

    if (hasEnding(config.graph_filename, bin_ending) || hasEnding(config.graph_filename, parhip_ending)) {
        while (remaining_nodes != 0) {
            unsigned long long nodes_in_batch = static_cast<unsigned long long>(batch_size);
            unsigned long long *vertex_offsets = new unsigned long long[nodes_in_batch + 1];
            (*in).seekg(start_pos);
            (*in).read((char *) (vertex_offsets), (nodes_in_batch + 1) * sizeof(unsigned long long));
            unsigned long long next_pos = start_pos + (nodes_in_batch) * sizeof(unsigned long long);

            unsigned long long edge_start_pos = vertex_offsets[0];
            unsigned long long num_reads = vertex_offsets[nodes_in_batch] - vertex_offsets[0];
            unsigned long long num_edges_to_read = num_reads / sizeof(unsigned long long);
            unsigned long long *edges = new unsigned long long[num_edges_to_read]; // we also need the next vertex offset
            (*in).seekg(edge_start_pos);
            (*in).read((char *) (edges), (num_edges_to_read) * sizeof(unsigned long long));
            start_pos = next_pos;

            unsigned long long pos = 0;
            for (unsigned long long i = 0; i < nodes_in_batch; ++i) {
                node = static_cast<NodeID>(i);
                unsigned long long degree = (vertex_offsets[i + 1] - vertex_offsets[i]) / sizeof(unsigned long long);
                for (unsigned long long j = 0; j < degree; j++, pos++) {
                    auto target = static_cast<NodeID>(edges[pos]);
                    if (target >= upper_node) {
                        // edge to future batch
                        continue;
                    } else if (target < lower_node) {
                        // edge to past batch
                        part_ids_on_nodes[node + lower_node][(*config.stream_nodes_assign)[incremental_edge_counter]]++;
                        part_ids_on_nodes[target][(*config.stream_nodes_assign)[incremental_edge_counter]]++;
                        incremental_edge_counter++;
                    } else {
                        // edge to current batch
                        if (node + lower_node < target) {
                            part_ids_on_nodes[node +
                                              lower_node][(*config.stream_nodes_assign)[incremental_edge_counter]]++;
                            part_ids_on_nodes[target][(*config.stream_nodes_assign)[incremental_edge_counter]]++;
                            incremental_edge_counter++;
                        }
                    }
                }
            }
            delete vertex_offsets;
            delete edges;
            if (incremental_edge_counter != 0) {
                lower_node = upper_node;
                remaining_nodes = remaining_nodes - batch_size;
                batch_size = MIN(config.stream_buffer_len, remaining_nodes);
                upper_node = lower_node + batch_size;
            }
        }
    } else {
        NodeID node_counter = 0;
        while (remaining_nodes != 0) {
            std::vector <std::string> *lines;
            lines = new std::vector<std::string>(1);
            buffered_input *ss2 = NULL;
            while (node_counter < batch_size) {
                std::getline(*(in), (*lines)[0]);
                if ((*lines)[0][0] == '%') { // a comment in the file
                    continue;
                }
                std::vector <LongNodeID> input;
                ss2 = new buffered_input(lines);
                ss2->simple_scan_line(input);
                (*lines)[0].clear();
                delete ss2;

                NodeID col_counter = 0;
                node = node_counter;
                while (col_counter < input.size()) {
                    target = input[col_counter++] - 1;
                    if (target >= upper_node) {
                        // edge to future batch
                        continue;
                    } else if (target < lower_node) {
                        // edge to past batch
                        part_ids_on_nodes[node + lower_node][(*config.stream_nodes_assign)[incremental_edge_counter]]++;
                        part_ids_on_nodes[target][(*config.stream_nodes_assign)[incremental_edge_counter]]++;
                        incremental_edge_counter++;
                    } else {
                        // edge to current batch
                        if (node + lower_node < target) {
                            part_ids_on_nodes[node +
                                              lower_node][(*config.stream_nodes_assign)[incremental_edge_counter]]++;
                            part_ids_on_nodes[target][(*config.stream_nodes_assign)[incremental_edge_counter]]++;
                            incremental_edge_counter++;
                        }
                    }
                }
                node_counter++;
            }
            delete lines;
            node_counter = 0;
            lower_node = upper_node;
            remaining_nodes = remaining_nodes - batch_size;
            batch_size = MIN(config.stream_buffer_len, remaining_nodes);
            upper_node = lower_node + batch_size;
        }
    }
    delete in;

    // iterate over all nodes
    std::vector <NodeID> v_ei(config.k, 0);
    std::vector <NodeID> block_weights(config.k, 0);

    for (auto &part_ids_on_node: part_ids_on_nodes) {
        bool set = false;
        NodeID blocks_on_node_count = 0;
        PartitionID cur_part = 0;
        for (auto &j: part_ids_on_node) {
            if (j > 0) {
                block_weights[cur_part] += j;
                v_ei[cur_part]++;
                replicas++;
                set = true;
                blocks_on_node_count++;
            }
            cur_part++;
        }
        if (set) {
            --replicas;
        }
        if (blocks_on_node_count > 1) {
            ++vertex_cut;
        }
    }

    double sum = 0;
    double max = -1;
    double total_weight = 0;
    for (int i = 0; i < v_ei.size(); i++) {
        sum += v_ei[i];
        if (block_weights[i] / 2 > max) {
            max = block_weights[i] / 2;
        }
        if (block_weights[i] / 2 > config.stream_total_upperbound) {
            std::cout << "Partition is imbalanced" << std::endl;
        }
        total_weight += block_weights[i] / 2;
    }
    replication_factor = sum / (double) nmbNodes;
    double balance_part_weight = ceil(total_weight / (double) config.k);
    balance = max / balance_part_weight;
}

void graph_io_stream::buildGraphModel(graph_access &G, PartitionConfig &partition_config,
                                      std::vector <EdgeID> &contracted_edge_graph_id,
                                      std::vector <EdgeID> &rev_edge, EdgeID S_nodes, NodeID S_edges,
                                      graph_access *G_temp) {
    // Directly constructs graph model G from contracted split graph (S) of graph
    // G_temp.

    partition_config.nmbNodes = S_nodes / 2;
    partition_config.prev_batch_edge_ID = S_nodes / 2;
    partition_config.remaining_stream_edges = S_edges / 2;

    if (partition_config.batch_alpha) {
        partition_config.fennel_alpha = partition_config.remaining_stream_edges *
                                        std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
                                        (std::pow(partition_config.remaining_stream_nodes,
                                                  partition_config.fennel_gamma));

        partition_config.fennel_alpha_gamma = partition_config.fennel_alpha * partition_config.fennel_gamma;
    }

    NodeID node_counter = 0;
    EdgeID edge_counter = 0;
    NodeID node = 0;
    LongNodeID target;
    NodeWeight total_nodeweight = 0;
    LongEdgeID used_edges = 0;
    LongEdgeID nmbEdges = 2 * partition_config.remaining_stream_edges;
    NodeWeight weight;
    std::vector < std::vector < std::pair < NodeID, EdgeWeight>>> all_edges;
    std::vector <NodeWeight> all_nodes;
    all_edges.resize(partition_config.nmbNodes + partition_config.quotient_nodes);
    all_nodes.resize(partition_config.nmbNodes + partition_config.quotient_nodes);
    partition_config.lower_global_node =
            partition_config.total_stream_nodecounter + 1; // Bounds below start from 1 instead of 0
    partition_config.upper_global_node = partition_config.total_stream_nodecounter + partition_config.nmbNodes;
    partition_config.curr_batch++;

    if (nmbEdges > std::numeric_limits<EdgeWeight>::max() ||
        partition_config.nmbNodes > std::numeric_limits<LongNodeID>::max()) {
#ifdef MODE64BITEDGES
        std::cerr << "The graph is too large. Currently only 64bits supported!"
                  << std::endl;
#else
        std::cerr << "The graph is too large. Currently only 32bits supported!" << std::endl;
#endif
        exit(0);
    }

    partition_config.edge_block_nodes = new std::vector <std::vector<std::pair < NodeID, NodeWeight>> >
                                        (partition_config.k, std::vector < std::pair < NodeID, NodeWeight >> ());

    NodeID in_deg = 0;
    NodeID out_deg = 0;
    NodeID current_edge_node;

    for (NodeID u = 0; u < G_temp->number_of_nodes(); ++u) {
        NodeID first_edge_u = G_temp->get_first_edge(u);
        for (EdgeID e = G_temp->get_first_edge(u); e < G_temp->get_first_invalid_edge(u); ++e) {
            NodeID v = G_temp->getEdgeTarget(e);
            if (u < v) {
                NodeID first_edge_v = G_temp->get_first_edge(v);
                in_deg = G_temp->getNodeDegree(u);
                out_deg = G_temp->getNodeDegree(v);
                current_edge_node = contracted_edge_graph_id[e];
                EdgeWeight edge_weight = 1;
                node_counter++;
                weight = 1;
                node = current_edge_node;

                all_nodes[node] = weight;

                if (in_deg <= 1) {
                    // do nothing
                } else if (in_deg > 1) {
                    EdgeID u_next = (in_deg + e - first_edge_u + 1) % in_deg + first_edge_u;
                    EdgeID u_prev = (in_deg + e - first_edge_u - 1) % in_deg + first_edge_u;

                    target = contracted_edge_graph_id[u_next] + 1 + partition_config.last_edge_count;
                    used_edges += ((NodeID)(target - partition_config.lower_global_node) <
                                   node); // used_edges only counts arcs to previous nodes
                    edge_counter += insertRegularEdgeInBatchEdge(partition_config, all_edges, node, target,
                                                                 edge_weight);
                    if (u_next != u_prev) {
                        target = contracted_edge_graph_id[u_prev] + 1 + partition_config.last_edge_count;
                        used_edges += ((NodeID)(target - partition_config.lower_global_node) <
                                       node); // used_edges only counts arcs to previous nodes
                        edge_counter += insertRegularEdgeInBatchEdge(partition_config, all_edges, node, target,
                                                                     edge_weight);
                    }
                }

                if (out_deg <= 1) {
                    // do nothing
                } else if (out_deg > 1) {
                    EdgeID u_next = (out_deg + rev_edge[e] - first_edge_v + 1) % out_deg + first_edge_v;
                    EdgeID u_prev = (out_deg + rev_edge[e] - first_edge_v - 1) % out_deg + first_edge_v;

                    target = contracted_edge_graph_id[u_next] + 1 + partition_config.last_edge_count;
                    used_edges += ((NodeID)(target - partition_config.lower_global_node) <
                                   node); // used_edges only counts arcs to previous nodes
                    edge_counter += insertRegularEdgeInBatchEdge(partition_config, all_edges, node, target,
                                                                 edge_weight);
                    if (u_next != u_prev) {
                        target = contracted_edge_graph_id[u_prev] + 1 + partition_config.last_edge_count;
                        used_edges += ((NodeID)(target - partition_config.lower_global_node) <
                                       node); // used_edges only counts arcs to previous nodes
                        edge_counter += insertRegularEdgeInBatchEdge(partition_config, all_edges, node, target,
                                                                     edge_weight);
                    }
                }

                addEdgesToIncidentBlock(partition_config, node, target, 1, used_edges);
                if ((partition_config.minimal_mode && partition_config.quotient_edges_count > 1) ||
                    (partition_config.past_subset_size != -1 &&
                     partition_config.quotient_edges_count > partition_config.past_subset_size)) {
                    std::cout << "Added: " << partition_config.quotient_edges_count << std::endl;
                }
                partition_config.quotient_edges_count = 0;
            }
        }
    }

    delete G_temp;
    insertQuotientNodesInBatchEdge(partition_config, all_nodes, node_counter);
    edge_counter += insertQuotientEdgesInBatchEdge(partition_config, all_edges);

    createGraphForBatchEdge(partition_config, G, node_counter, edge_counter, all_edges, all_nodes);

    delete partition_config.edge_block_nodes;

    partition_config.total_stream_nodecounter += partition_config.nmbNodes;
    partition_config.total_stream_nodeweight += total_nodeweight;
    partition_config.remaining_stream_nodes -= partition_config.nmbNodes;
    partition_config.remaining_stream_edges -= used_edges;

    if (node_counter != (NodeID)
        partition_config.nmbNodes + partition_config.quotient_nodes) {
        std::cerr << "number of specified nodes mismatch" << std::endl;
        std::cerr << (partition_config.nmbNodes + partition_config.quotient_nodes) << " " << node_counter << std::endl;
        exit(0);
    }
}

void graph_io_stream::updateBatchVariables(PartitionConfig &partition_config) {
    // Update values required to map nodes between global node IDs to per-batch
    // node IDs
    if (partition_config.total_stream_nodecounter == 0) {
        partition_config.lower_global_node_conv = 0;
    } else {
        partition_config.lower_global_node_conv = partition_config.upper_global_node_conv;
    }
    partition_config.upper_global_node_conv = partition_config.lower_global_node_conv + partition_config.nmbNodes;
    partition_config.last_edge_count = partition_config.incremental_edge_ID;
    partition_config.prev_batch_edge_ID = partition_config.last_edge_count;
    partition_config.remaining_stream_graph_nodes -= partition_config.nmbNodes;
}

void graph_io_stream::constructBatchIOGraph(PartitionConfig &partition_config,
                                            std::vector <std::vector<Edge>> &subgraph_edges,
                                            std::vector <NodeID> &mapping, std::vector <EdgeID> &rev_edge,
                                            std::vector <EdgeID> &contracted_edge_graph_id,
                                            NodeID &number_of_deg1_vertices, NodeID &number_of_deg2_vertices,
                                            NodeID &G_temp_no_nodes, NodeID &G_temp_no_edges,
                                            graph_access &G_temp) {
    // store graph built from current batch nodes and edges in G_temp
    // G_temp will later be turned into the CSPAC per-batch model
    // the function also builds rev_edge and contracted_edge_graph_id
    number_of_deg1_vertices = 0;
    number_of_deg2_vertices = 0;
    NodeID edge_to_node_key = 0;
    NodeID u_OG, v_OG;
    contracted_edge_graph_id.resize(G_temp_no_edges, -1);
    EdgeID edge_index = 0;

    if (partition_config.nodes_on_edge_conv != NULL) {
        delete partition_config.nodes_on_edge_conv;
        //partition_config.nodes_on_edge_conv->clear();
    }
    partition_config.nodes_on_edge_conv = new std::vector <std::vector<NodeID>>(G_temp_no_edges / 2,
                                                                                std::vector<NodeID>(2));

    rev_edge.resize(G_temp_no_edges);

    G_temp.start_construction_light(G_temp_no_nodes, G_temp_no_edges);
    for (NodeID u = 0; u < G_temp_no_nodes; ++u) {
        if (subgraph_edges[u].size() == 1) {
            ++number_of_deg1_vertices;
        } else if (subgraph_edges[u].size() == 2) {
            ++number_of_deg2_vertices;
        }
        u_OG = mapping[u];

        G_temp.new_node();
        G_temp.setNodeWeight(u, 1);
        for (auto &j: subgraph_edges[u]) {
            NodeID target_idx = j.target;
            v_OG = mapping[target_idx];
            const EdgeID e = G_temp.new_edge(u, target_idx);
            if (u < target_idx) { // if left to right edge, push edge ID to index j
                j.weight = e;
                edge_to_node_key = partition_config.incremental_edge_ID++;
                contracted_edge_graph_id[e] = edge_to_node_key - partition_config.last_edge_count;
            } else { // if right to left edge, fetch edge ID from index u (which
                // was j earlier) of corresponding left to right edge
                EdgeWeight rev_edge_ID = subgraph_edges[target_idx][j.weight].weight;
                rev_edge[rev_edge_ID] = e;
                rev_edge[e] = rev_edge_ID;
                contracted_edge_graph_id[e] = contracted_edge_graph_id[rev_edge[e]];
            }

            if (u_OG < v_OG) {
                (*partition_config.nodes_on_edge_conv)[contracted_edge_graph_id[e]][0] = u_OG;
                (*partition_config.nodes_on_edge_conv)[contracted_edge_graph_id[e]][1] = v_OG;
            }
            edge_index++;
            G_temp.setEdgeWeight(e, 1);
        }
    }

    G_temp.finish_construction_light();
}

void graph_io_stream::readInputAsGraph(PartitionConfig &partition_config, std::vector <EdgeID> &rev_edge,
                                       std::vector <EdgeID> &contracted_edge_graph_id,
                                       NodeID &number_of_deg1_vertices, NodeID &number_of_deg2_vertices,
                                       graph_access &G_temp) {
    // Read graph IO and store in G_temp

    updateBatchVariables(partition_config);
    timer fit;

    std::vector <std::vector<Edge>> G_temp_edges;
    G_temp_edges.resize(partition_config.nmbNodes); // stores edges in G_temp
    std::vector <NodeID> mapping(partition_config.nmbNodes, -1);
    LongNodeID target;

    NodeID target_idx = -1; // store mapping of target nodes in current batch
    LongNodeID node_counter = 0;
    NodeID G_temp_no_nodes = partition_config.nmbNodes;
    NodeID G_temp_no_edges = 0;
    NodeID map_index = partition_config.nmbNodes;
    NodeID node = 0;
    google::dense_hash_map <NodeID, NodeID> prev_batch_mapping; // remember which nodes from current to previous batch
    // have been mapped
    prev_batch_mapping.set_empty_key(-1);

    fit.restart();
    std::vector <std::string> *lines;
    lines = new std::vector<std::string>(1);
    buffered_input *ss2 = NULL;
    partition_config.read_graph_time += fit.elapsed();

    while (node_counter < partition_config.nmbNodes) {
        fit.restart();
        std::getline(*(partition_config.stream_in), (*lines)[0]);
        if ((*lines)[0][0] == '%') { // a comment in the file
            continue;
        }
        std::vector <LongNodeID> input;
        ss2 = new buffered_input(lines);
        ss2->simple_scan_line(input);
        (*lines)[0].clear();
        delete ss2;
        partition_config.read_graph_time += fit.elapsed();

        LongNodeID col_counter = 0;
        node = (NodeID)
        node_counter;
        mapping[node] = node + partition_config.lower_global_node_conv;

        while (col_counter < input.size()) {
            target = input[col_counter++] - 1;

            if (target >= partition_config.upper_global_node_conv) {
                // ignore edge to future batch
                continue;
            } else if (target < partition_config.lower_global_node_conv) {
                // edge to previous batch

                if (prev_batch_mapping.find(target) == prev_batch_mapping.end()) {
                    // if target has not been mapped, assign it to a mapping
                    prev_batch_mapping[target] = map_index;
                    target_idx = map_index;
                    G_temp_edges.emplace_back();
                    mapping.push_back(target);
                    map_index++;
                    G_temp_no_nodes++;
                } else {
                    target_idx = prev_batch_mapping[target];
                }

                // add forward and backward edge to previous batch
                Edge forward_edge;
                forward_edge.target = target_idx;
                forward_edge.weight = -1;
                G_temp_edges[node].emplace_back(forward_edge);

                Edge backward_edge;
                backward_edge.target = node;
                backward_edge.weight = G_temp_edges[node].size() - 1;
                G_temp_edges[target_idx].emplace_back(backward_edge);

                G_temp_no_edges += 2;
            } else {
                // edge to current batch

                target_idx = target - partition_config.lower_global_node_conv;
                if (node < target_idx) {
                    Edge forward_edge;
                    forward_edge.target = target_idx;
                    forward_edge.weight = -1;
                    G_temp_edges[node].emplace_back(forward_edge);

                    Edge backward_edge;
                    backward_edge.target = node;
                    backward_edge.weight = G_temp_edges[node].size() - 1;
                    G_temp_edges[target_idx].emplace_back(backward_edge);

                    G_temp_no_edges += 2;
                }
            }
        }

        node_counter++;
    }
    delete lines;

    timer aet;
    aet.restart();

    // build graph from G_temp_edges
    constructBatchIOGraph(partition_config, G_temp_edges, mapping, rev_edge, contracted_edge_graph_id,
                          number_of_deg1_vertices, number_of_deg2_vertices, G_temp_no_nodes, G_temp_no_edges, G_temp);

    partition_config.assign_edge_ID_time += aet.elapsed();
}

void graph_io_stream::readInputAsGraphBinary(PartitionConfig &partition_config, std::vector <EdgeID> &rev_edge,
                                             std::vector <EdgeID> &contracted_edge_graph_id,
                                             NodeID &number_of_deg1_vertices, NodeID &number_of_deg2_vertices,
                                             graph_access &G_temp) {
    // Read graph IO and store in G_temp

    updateBatchVariables(partition_config);
    timer fit;

    std::vector <std::vector<Edge>> G_temp_edges;
    G_temp_edges.resize(partition_config.nmbNodes); // stores edges in G_temp
    std::vector <NodeID> mapping(partition_config.nmbNodes, -1);

    NodeID target_idx = -1; // store mapping of target nodes in current batch
    NodeID G_temp_no_nodes = partition_config.nmbNodes;
    NodeID G_temp_no_edges = 0;
    NodeID map_index = partition_config.nmbNodes;
    NodeID node = 0;
    google::dense_hash_map <NodeID, NodeID> prev_batch_mapping; // remember which nodes from current to previous batch
    // have been mapped
    prev_batch_mapping.set_empty_key(-1);

    fit.restart();

    unsigned long long nodes_in_batch = static_cast<unsigned long long>(partition_config.nmbNodes);
    unsigned long long *vertex_offsets = new unsigned long long[nodes_in_batch + 1];
    (*partition_config.stream_in).seekg(partition_config.start_pos);
    (*partition_config.stream_in).read((char *) (vertex_offsets), (nodes_in_batch + 1) * sizeof(unsigned long long));
    unsigned long long next_pos = partition_config.start_pos + (nodes_in_batch) * sizeof(unsigned long long);

    unsigned long long edge_start_pos = vertex_offsets[0];
    unsigned long long num_reads = vertex_offsets[nodes_in_batch] - vertex_offsets[0];
    unsigned long long num_edges_to_read = num_reads / sizeof(unsigned long long);
    unsigned long long *edges = new unsigned long long[num_edges_to_read]; // we also need the next vertex offset
    (*partition_config.stream_in).seekg(edge_start_pos);
    (*partition_config.stream_in).read((char *) (edges), (num_edges_to_read) * sizeof(unsigned long long));

    partition_config.read_graph_time += fit.elapsed();

    partition_config.start_pos = next_pos;

    unsigned long long pos = 0;
    for (unsigned long long i = 0; i < nodes_in_batch; ++i) {
        node = static_cast<NodeID>(i);
        unsigned long long degree = (vertex_offsets[i + 1] - vertex_offsets[i]) / sizeof(unsigned long long);
        G_temp_edges[node].reserve(degree);
        mapping[node] = node + partition_config.lower_global_node_conv;
        for (unsigned long long j = 0; j < degree; j++, pos++) {
            auto target = static_cast<NodeID>(edges[pos]);

            if (target >= partition_config.upper_global_node_conv) {
                // ignore edge to future batch
                continue;
            } else if (target < partition_config.lower_global_node_conv) {
                // edge to previous batch

                if (prev_batch_mapping.find(target) == prev_batch_mapping.end()) {
                    // if target has not been mapped, assign it to a mapping
                    prev_batch_mapping[target] = map_index;
                    target_idx = map_index;
                    G_temp_edges.emplace_back();
                    mapping.push_back(target);
                    map_index++;
                    G_temp_no_nodes++;
                } else {
                    target_idx = prev_batch_mapping[target];
                }

                // add forward and backward edge to previous batch
                Edge forward_edge{};
                forward_edge.target = target_idx;
                forward_edge.weight = -1;
                G_temp_edges[node].emplace_back(forward_edge);

                Edge backward_edge{};
                backward_edge.target = node;
                backward_edge.weight = G_temp_edges[node].size() - 1;
                G_temp_edges[target_idx].emplace_back(backward_edge);

                G_temp_no_edges += 2;
            } else {
                // edge to current batch

                target_idx = target - partition_config.lower_global_node_conv;
                if (node < target_idx) {
                    Edge forward_edge{};
                    forward_edge.target = target_idx;
                    forward_edge.weight = -1;
                    G_temp_edges[node].emplace_back(forward_edge);

                    Edge backward_edge{};
                    backward_edge.target = node;
                    backward_edge.weight = G_temp_edges[node].size() - 1;
                    G_temp_edges[target_idx].emplace_back(backward_edge);

                    G_temp_no_edges += 2;
                }
            }
        }

    }
    delete vertex_offsets;
    delete edges;

    timer aet;
    aet.restart();

    // build graph from G_temp_edges
    constructBatchIOGraph(partition_config, G_temp_edges, mapping, rev_edge, contracted_edge_graph_id,
                          number_of_deg1_vertices, number_of_deg2_vertices, G_temp_no_nodes, G_temp_no_edges, G_temp);

    partition_config.assign_edge_ID_time += aet.elapsed();
}

void graph_io_stream::constructBatchModel(PartitionConfig &partition_config, graph_access &G) {
    //******construct graph model from input******

    graph_access *G_temp = new graph_access(); // stores graph constructed from nodes and edges in current buffer
    std::vector <EdgeID> rev_edge; // resized later to number of edges in current batch
    // stores edge ID of edge v->u in index of edge u->v where v>u
    std::vector <EdgeID> contracted_edge_graph_id; // resized later to number of edges in current batch
    // stores for each edge of G_temp the corresponding node ID
    // in the constructed graph model
    NodeID number_of_deg1_vertices, number_of_deg2_vertices;

    // readInputAsGraph builds G_temp, rev_edge, and contracted_edge_graph_id
    if (hasEnding(partition_config.graph_filename, ".bin") || hasEnding(partition_config.graph_filename, ".parhip")) {
        readInputAsGraphBinary(partition_config, rev_edge, contracted_edge_graph_id, number_of_deg1_vertices,
                               number_of_deg2_vertices, *G_temp);
    } else {
        readInputAsGraph(partition_config, rev_edge, contracted_edge_graph_id, number_of_deg1_vertices,
                         number_of_deg2_vertices, *G_temp);
    }


    //******construct contracted split graph from graph******

    // number of nodes and edges of split graph
    NodeID split_n = G_temp->number_of_edges(); // two times the number of edges
    EdgeID split_m = 3 * G_temp->number_of_edges() - 2 * (number_of_deg1_vertices + number_of_deg2_vertices);

    // update dynamic alpha if using that mode
    if (partition_config.dynamic_alpha) {
        partition_config.fennel_edges -= (2 * (number_of_deg2_vertices + number_of_deg1_vertices) + split_n);
        partition_config.fennel_alpha =
                partition_config.fennel_edges * std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
                (std::pow(partition_config.remaining_stream_nodes, partition_config.fennel_gamma));
        partition_config.fennel_alpha_gamma = partition_config.fennel_alpha * partition_config.fennel_gamma;
    }
    timer gct;
    gct.restart();


    buildGraphModel(G, partition_config, contracted_edge_graph_id, rev_edge, split_n, split_m - split_n, G_temp);

    partition_config.graph_model_time += gct.elapsed();

}

void graph_io_stream::writeStreamOutput(PartitionConfig &config, graph_access &G_local) {
    timer sot;
    sot.restart();
    for (NodeID node = 0, end = config.nmbNodes; node < end; node++) {
        PartitionID block = G_local.getPartitionIndex(node);
        (*config.stream_out) << block << "\n";
    }
    config.stream_output_time += sot.elapsed();
}