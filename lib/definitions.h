/******************************************************************************
 * definitions.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef DEFINITIONS_H_CHR
#define DEFINITIONS_H_CHR

#include <limits>
#include <queue>
#include <vector>
#include <unordered_map>
#include <list>

#include "limits.h"
#include "macros_assertions.h"
#include "stdio.h"


// allows us to disable most of the output during partitioning
#ifdef KAFFPAOUTPUT
        #define PRINT(x) x
#else
        #define PRINT(x) do {} while (false);
#endif

/**********************************************
 * Constants
 * ********************************************/
//Types needed for the graph ds
typedef unsigned int 	NodeID;
typedef double 		EdgeRatingType;
typedef unsigned int 	PathID;
typedef unsigned int 	PartitionID;
#ifdef MODE64BITEDGES
typedef int64_t 	EdgeWeight;
typedef uint64_t 	NodeWeight;
typedef uint64_t 	EdgeID;
#else
typedef int 		EdgeWeight;
typedef unsigned int 	NodeWeight;
typedef unsigned int 	EdgeID;
#endif
typedef EdgeWeight 	Gain;
typedef int 		Color;
typedef unsigned int 	Count;
typedef std::vector<NodeID> boundary_starting_nodes;
typedef long FlowType;

const EdgeID UNDEFINED_EDGE            = std::numeric_limits<EdgeID>::max();
const EdgeWeight UNDEFINED_EDGEWEIGHT  = std::numeric_limits<EdgeWeight>::max();
const NodeID UNDEFINED_NODE            = std::numeric_limits<NodeID>::max();
const NodeID UNASSIGNED                = std::numeric_limits<NodeID>::max();
const NodeID ASSIGNED                  = std::numeric_limits<NodeID>::max()-1;
const PartitionID INVALID_PARTITION    = std::numeric_limits<PartitionID>::max();
const PartitionID BOUNDARY_STRIPE_NODE = std::numeric_limits<PartitionID>::max();
const int NOTINQUEUE 		       = std::numeric_limits<int>::max();
const int ROOT 			       = 0;

//for the gpa algorithm
struct edge_source_pair {
        EdgeID e;
        NodeID source;       
};

struct source_target_pair {
        NodeID source;       
        NodeID target;       
};

//matching array has size (no_of_nodes), so for entry in this table we get the matched neighbor
typedef std::vector<NodeID> CoarseMapping;
typedef std::vector<NodeID> Matching;
typedef std::vector<NodeID> NodePermutationMap;

typedef double ImbalanceType;
//Coarsening
typedef enum {
        EXPANSIONSTAR, 
        EXPANSIONSTAR2, 
 	WEIGHT, 
 	REALWEIGHT, 
	PSEUDOGEOM, 
	EXPANSIONSTAR2ALGDIST, 
        SEPARATOR_MULTX,
        SEPARATOR_ADDX,
        SEPARATOR_MAX,
        SEPARATOR_LOG,
        SEPARATOR_R1,
        SEPARATOR_R2,
        SEPARATOR_R3,
        SEPARATOR_R4,
        SEPARATOR_R5,
        SEPARATOR_R6,
        SEPARATOR_R7,
        SEPARATOR_R8
} EdgeRating;

typedef enum {
        PERMUTATION_QUALITY_NONE, 
	PERMUTATION_QUALITY_FAST,  
	PERMUTATION_QUALITY_GOOD
} PermutationQuality;

typedef enum {
        MATCHING_RANDOM, 
	MATCHING_GPA, 
	MATCHING_RANDOM_GPA,
        CLUSTER_COARSENING
} MatchingType;

typedef enum {
	INITIAL_PARTITIONING_RECPARTITION, 
	INITIAL_PARTITIONING_BIPARTITION,
	INITIAL_PARTITIONING_MULTIBFS,
	INITIAL_PARTITIONING_FENNEL
} InitialPartitioningType;

typedef enum {
        REFINEMENT_SCHEDULING_FAST, 
	REFINEMENT_SCHEDULING_ACTIVE_BLOCKS, 
	REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY
} RefinementSchedulingAlgorithm;

typedef enum {
        REFINEMENT_TYPE_FM, 
	REFINEMENT_TYPE_FM_FLOW, 
	REFINEMENT_TYPE_FLOW
} RefinementType;

typedef enum {
        STOP_RULE_SIMPLE, 
	STOP_RULE_MULTIPLE_K, 
	STOP_RULE_STRONG,
        STOP_RULE_MULTIBFS 
} StopRule;

typedef enum {
        BIPARTITION_BFS, 
	BIPARTITION_FM
} BipartitionAlgorithm ;

typedef enum {
        KWAY_SIMPLE_STOP_RULE, 
	KWAY_ADAPTIVE_STOP_RULE
} KWayStopRule;

typedef enum {
        COIN_RNDTIE, 
	COIN_DIFFTIE, 
	NOCOIN_RNDTIE, 
	NOCOIN_DIFFTIE 
} MLSRule;

typedef enum {
        CYCLE_REFINEMENT_ALGORITHM_PLAYFIELD, 
        CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL, 
	CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL_PLUS
} CycleRefinementAlgorithm;

typedef enum {
        RANDOM_NODEORDERING, 
        DEGREE_NODEORDERING,
	NATURAL_NODEORDERING
} NodeOrderingType; typedef enum {
        NSQUARE, 
        NSQUAREPRUNED, 
        COMMUNICATIONGRAPH
} LsNeighborhoodType;

typedef enum {
        MAP_CONST_RANDOM, 
        MAP_CONST_IDENTITY,
        MAP_CONST_OLDGROWING,
        MAP_CONST_OLDGROWING_FASTER,
        MAP_CONST_OLDGROWING_MATRIX,
        MAP_CONST_FASTHIERARCHY_BOTTOMUP,
        MAP_CONST_FASTHIERARCHY_TOPDOWN
} ConstructionAlgorithm;

typedef enum {
        DIST_CONST_RANDOM, 
        DIST_CONST_IDENTITY,
        DIST_CONST_HIERARCHY,
        DIST_CONST_HIERARCHY_ONLINE
} DistanceConstructionAlgorithm;

typedef enum {
        PRE_CONFIG_MAPPING_FAST, 
        PRE_CONFIG_MAPPING_ECO,
        PRE_CONFIG_MAPPING_STRONG
} PreConfigMapping;


//////////////////////////////////////////////////////
//////////////// Streaming algorithms ////////////////
////////////////////// (begin) ///////////////////////
//////////////////////////////////////////////////////

// One-pass streaming algorithms
typedef enum {
        ONEPASS_BALANCED, 
        ONEPASS_CHUNKING, 
        ONEPASS_HASHING,
	ONEPASS_GREEDY,
	ONEPASS_LDG,
	ONEPASS_FENNEL,
	ONEPASS_FRACTIONAL_GREEDY 
} OnePassStreamAlgorithm;

// Dynamic modification of Fennel objective function within label propagation
typedef enum {
        FENNELADP_ORIGINAL, 
	FENNELADP_DOUBLE,
        FENNELADP_LINEAR, 
        FENNELADP_QUADRATIC, 
        FENNELADP_MID_LINEAR, 
        FENNELADP_MID_QUADRATIC,
        FENNELADP_MID_CONSTANT,
	FENNELADP_EDGE_CUT
} FennelDynamicSetup;

// Order to compute Fennel initial solution inside each batch
typedef enum {
        FENN_ORDER_UNCHANGED, 
	FENN_ORDER_ASC_DEGREE,
        FENN_ORDER_DESC_DEGREE 
} FennelBatchOrder;

// Procedure for ghost neighbors
typedef enum {
        GHOST_CONTRACT_ALL, 
	GHOST_KEEP_ALL,
        GHOST_KEEP_THRESHOLD_CONTRACT_REST
} GhostNeighborsProcedure;

// Datatypes and constants for huge graphs
typedef uint64_t LongNodeID;
typedef uint64_t LongEdgeID;
typedef int32_t	ShortEdgeWeight;
const LongNodeID UNDEFINED_LONGNODE = std::numeric_limits<LongNodeID>::max();
const LongEdgeID UNDEFINED_LONGEDGE = std::numeric_limits<LongEdgeID>::max();

//////////////////////////////////////////////////////
//////////////// Streaming algorithms ////////////////
/////////////////////// (end) ////////////////////////
//////////////////////////////////////////////////////

// Conversion of graphs
typedef enum {
	EDGE_STREAM_TO_METISEXTERNAL,
        EDGE_STREAM_TO_METISBUFFERED, 
        EDGE_STREAM_TO_METIS, 
	EDGE_STREAM_TO_HMETIS,
        METIS_TO_HMETIS 
} GraphTranslateFormat;

 
//////////////////////////////////////////////////////
////////////////// Hashing Function //////////////////
////////////////////// (begin) ///////////////////////
//////////////////////////////////////////////////////

// default values recommended by http://isthe.com/chongo/tech/comp/fnv/
const uint32_t Prime = 0x01000193; //   16777619
const uint32_t Seed  = 0x811C9DC5; // 2166136261
// hash a single byte
inline uint32_t fnv0a(unsigned char oneByte, uint32_t hash = Seed) {
	return (oneByte ^ hash) * Prime;
}
// hash a 32 bit integer (four bytes)
inline uint32_t fnv1a(uint32_t fourBytes, uint32_t hash = Seed) {
	const unsigned char* ptr = (const unsigned char*) &fourBytes;
	hash = fnv0a(*ptr++, hash);
	hash = fnv0a(*ptr++, hash);
	hash = fnv0a(*ptr++, hash);
	return fnv0a(*ptr  , hash);
}
// hash a 64 bit integer (eight bytes)
inline uint32_t fnv2a(uint64_t eightBytes, uint32_t hash = Seed) {
	const unsigned char* ptr = (const unsigned char*) &eightBytes;
	hash = fnv0a(*ptr++, hash);
	hash = fnv0a(*ptr++, hash);
	hash = fnv0a(*ptr++, hash);
	hash = fnv0a(*ptr++, hash);
	hash = fnv0a(*ptr++, hash);
	hash = fnv0a(*ptr++, hash);
	hash = fnv0a(*ptr++, hash);
	return fnv0a(*ptr  , hash);
}
//////////////////////////////////////////////////////
////////////////// Hashing Function //////////////////
/////////////////////// (end) ////////////////////////
//////////////////////////////////////////////////////

#endif

