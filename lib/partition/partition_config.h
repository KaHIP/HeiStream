/******************************************************************************
 * partition_config.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARTITION_CONFIG_DI1ES4T0
#define PARTITION_CONFIG_DI1ES4T0

#include <fstream>
#include <map>
#include <memory>
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include <unordered_set>

#include "algorithms/node_partitioning/execution/batch_id_manager.h"
#include "data_structure/buffered_map.h"
#include "data_structure/graph_access.h"
#include "definitions.h"

typedef struct {
    PartitionID block;
    int gain;
    int degree;
} DELTA;

class matrix;
namespace partition {
namespace state {
struct NodePartitionerPassState;
struct EdgePartitionerPassState;
}  // namespace state
}  // namespace partition

// Global runtime configuration for streaming algorithms.
struct Config {
    Config() {}

    // NOTE:
    // `heistream` and `heistream_edge` are streaming executables, but they still
    // reuse portions of the multilevel KaHIP stack during batch/model solving.
    // The sections above stream-specific fields are retained for that runtime
    // compatibility and should be treated as legacy compatibility knobs.

    //============================================================
    //=======================MATCHING=============================
    //============================================================
    bool edge_rating_tiebreaking;

    EdgeRating edge_rating;

    PermutationQuality permutation_quality;

    MatchingType matching_type;

    bool match_islands;

    bool first_level_random_matching;

    bool rate_first_level_inner_outer;

    NodeWeight max_vertex_weight;

    NodeWeight largest_graph_weight;

    NodeWeight work_load;

    unsigned aggressive_random_levels;

    bool disable_max_vertex_weight_constraint;

    //============================================================
    //===================INITIAL PARTITIONING=====================
    //============================================================
    unsigned int initial_partitioning_repetitions;

    unsigned int minipreps;

    bool refined_bubbling;

    InitialPartitioningType initial_partitioning_type;

    bool initial_partition_optimize;

    BipartitionAlgorithm bipartition_algorithm;

    bool initial_partitioning;

    int bipartition_tries;

    int bipartition_post_fm_limits;

    int bipartition_post_ml_limits;

    //============================================================
    //====================REFINEMENT PARAMETERS===================
    //============================================================
    bool corner_refinement_enabled;

    bool use_bucket_queues;

    RefinementType refinement_type;

    PermutationQuality permutation_during_refinement;

    ImbalanceType imbalance;

    unsigned bubbling_iterations;

    unsigned kway_rounds;

    bool quotient_graph_refinement_disabled;

    KWayStopRule kway_stop_rule;

    double kway_adaptive_limits_alpha;

    double kway_adaptive_limits_beta;

    unsigned max_flow_iterations;

    unsigned local_multitry_rounds;

    unsigned local_multitry_fm_alpha;

    bool graph_allready_partitioned;

    unsigned int fm_search_limit;

    unsigned int kway_fm_search_limit;

    NodeWeight upper_bound_partition;

    double bank_account_factor;

    RefinementSchedulingAlgorithm refinement_scheduling_algorithm;

    bool most_balanced_minimum_cuts;

    bool most_balanced_minimum_cuts_node_sep;

    unsigned toposort_iterations;

    bool softrebalance;

    bool rebalance;

    double flow_region_factor;

    bool gpa_grow_paths_between_blocks;

    //=======================================
    //==========GLOBAL SEARCH PARAMETERS=====
    //=======================================
    unsigned global_cycle_iterations;

    bool use_wcycles;

    bool use_fullmultigrid;

    unsigned level_split;

    bool no_new_initial_partitioning;

    bool omit_given_partitioning;

    StopRule stop_rule;

    int num_vert_stop_factor;

    bool no_change_convergence;

    //=======================================
    //===PERFECTLY BALANCED PARTITIONING ====
    //=======================================
    bool remove_negative_cycles;

    bool kaba_include_removal_of_paths;

    bool kaba_enable_zero_weight_cycles;

    double kabaE_internal_bal;

    CycleRefinementAlgorithm cycle_refinement_algorithm;

    int kaba_internal_no_aug_steps_aug;

    unsigned kaba_packing_iterations;

    bool kaba_flip_packings;

    MLSRule kaba_lsearch_p;  // more localized search pseudo directed

    bool kaffpa_perfectly_balanced_refinement;

    unsigned kaba_unsucc_iterations;

    //=======================================
    //============PAR_PSEUDOMH / MH =========
    //=======================================
    double time_limit;

    double epsilon;

    unsigned no_unsuc_reps;

    unsigned local_partitioning_repetitions;

    bool mh_plain_repetitions;

    bool mh_easy_construction;

    bool mh_enable_gal_combine;

    bool mh_no_mh;

    bool mh_print_log;

    int mh_flip_coin;

    int mh_initial_population_fraction;

    bool mh_disable_cross_combine;

    bool mh_cross_combine_original_k;

    bool mh_disable_nc_combine;

    bool mh_disable_combine;

    bool mh_enable_quickstart;

    bool mh_disable_diversify_islands;

    bool mh_diversify;

    bool mh_diversify_best;

    bool mh_enable_tournament_selection;

    bool mh_optimize_communication_volume;

    unsigned mh_num_ncs_to_compute;

    unsigned mh_pool_size;

    bool combine;  // in this case the second index is filled and edges between both partitions are
                   // not contracted

    unsigned initial_partition_optimize_fm_limits;

    unsigned initial_partition_optimize_multitry_fm_alpha;

    unsigned initial_partition_optimize_multitry_rounds;

    unsigned walshaw_mh_repetitions;

    unsigned scaleing_factor;

    bool scale_back;

    bool suppress_partitioner_output;

    unsigned maxT;

    unsigned maxIter;
    //=======================================
    //===============BUFFOON=================
    //=======================================
    bool disable_hard_rebalance;

    bool buffoon;

    bool kabapE;

    bool mh_penalty_for_unconnected;
    //=======================================
    //===============MISC====================
    //=======================================
    std::string input_partition;

    int seed;

    bool fast;

    bool eco;

    bool strong;

    bool kaffpaE;

    bool balance_edges;

    // number of blocks the graph should be partitioned in
    PartitionID k;

    bool compute_vertex_separator;

    bool only_first_level;

    bool use_balance_singletons;

    int amg_iterations;

    std::string graph_filename;

    std::string filename_output;

    bool kaffpa_perfectly_balance;

    bool mode_node_separators;

    //=======================================
    //===========SNW PARTITIONING============
    //=======================================
    NodeOrderingType node_ordering;

    int cluster_coarsening_factor;

    bool ensemble_clusterings;

    int label_iterations;

    int label_iterations_refinement;

    int number_of_clusterings;

    bool label_propagation_refinement;

    double balance_factor;

    bool cluster_coarsening_during_ip;

    bool set_upperbound;

    int repetitions;

    //=======================================
    //===========NODE SEPARATOR==============
    //=======================================
    int max_flow_improv_steps;

    int max_initial_ns_tries;

    double region_factor_node_separators;

    bool sep_flows_disabled;

    bool sep_fm_disabled;

    bool sep_loc_fm_disabled;

    int sep_loc_fm_no_snodes;

    bool sep_greedy_disabled;

    int sep_fm_unsucc_steps;

    int sep_loc_fm_unsucc_steps;

    int sep_num_fm_reps;

    int sep_num_loc_fm_reps;

    int sep_num_vert_stop;

    bool sep_full_boundary_ip;

    bool faster_ns;

    EdgeRating sep_edge_rating_during_ip;

    //=======================================
    //=========LABEL PROPAGATION=============
    //=======================================
    NodeWeight cluster_upperbound;

    //=======================================
    //=========INITIAL PARTITIONING==========
    //=======================================

    // variables controling the size of the blocks during
    // multilevel recursive bisection
    // (for the case where k is not a power of 2)
    std::vector<int> target_weights;

    bool initial_bipartitioning;

    int grow_target;

    //=======================================
    //===============QAP=====================
    //=======================================

    int communication_neighborhood_dist;

    LsNeighborhoodType ls_neighborhood;

    ConstructionAlgorithm construction_algorithm;

    DistanceConstructionAlgorithm distance_construction_algorithm;

    std::vector<int> group_sizes;

    std::vector<int> distances;

    int search_space_s;

    PreConfigMapping preconfiguration_mapping;

    int max_recursion_levels_construction;

    bool enable_mapping;

    //=======================================
    //===========integrated_mapping==========
    //=======================================

    bool integrated_mapping;
    bool multisection;
    bool qap_label_propagation_refinement;
    bool qap_blabel_propagation_refinement;
    bool qap_alabel_propagation_refinement;
    bool qap_multitry_kway_fm;
    bool qap_bmultitry_kway_fm;
    bool qap_kway_fm;
    bool qap_bkway_fm;
    bool qap_quotient_ref;
    bool qap_bquotient_ref;
    bool qap_0quotient_ref;
    bool bipartition_gp_local_search;
    bool skip_map_ls;
    bool suppress_output;
    bool no_change_convergence_map;
    bool full_matrix;
    matrix* D;
    std::vector<NodeID>* perm_rank;

    //=======================================
    //============= Delta gains =============
    //=======================================

    std::vector<std::pair<int, std::vector<DELTA*>>>* delta;
    std::vector<bool>* has_gains;
    bool use_delta_gains;
    bool quotient_more_mem;
    int* ref_layer;
    bool skip_delta_gains;

    //=======================================
    //======= Binary Online Distance ========
    //=======================================

    std::vector<std::vector<int>>* bin_id;
    bool use_bin_id;
    std::vector<unsigned int>* compact_bin_id;
    bool use_compact_bin_id;
    int bit_sec_len;
    int label_iterations_refinement_map;

    //=======================================
    //======== Adaptative Balancing =========
    //=======================================

    bool adapt_bal;
    double glob_block_upperbound;
    std::vector<int> interval_sizes;

    //=======================================
    //========== Stream Partition ===========
    //=======================================
    // Primary runtime configuration for `heistream`.

    // BuffCut / buffered streaming (node)
    bool run_parallel;
    float param_dbl1;
    bool param_enbld1;
    std::vector<std::vector<LongNodeID>>* batch_unpartitioned_neighbors;
    bool ghost_neighbors_enabled;
    EdgeWeight default_weight_non_ghost;
    float ghost_weight;
    float inv_ghost_weight;  // Inverse ghost weight for buffer score updates
    bool sep_batch_marker;
    LongNodeID num_ghost_nodes;
    bool restream_include_high_degree_nodes;
    BatchExtractionStrategy batch_extraction_strategy;
    BatchIDManager* batch_manager;
    size_t max_active_batches;
    bool print_times;
    bool part_adj_directly;
    size_t max_input_q_size;
    BPQStorageType bpq_storage_type;
    BufferScoreType buffer_score_type;
    LongNodeID d_max;
    float haa_beta;
    float haa_theta;
    float cbs_theta;
    LongNodeID max_block_weight;
    LongNodeID number_of_nodes;
    std::vector<NodeID>* local_to_global_map;
    unsigned bb_ratio;
    LongNodeID max_buffer_size;
    unsigned bq_disc_factor;

    bool stream_input;
    LongNodeID batch_size;
    LongNodeID total_nodes;
    LongEdgeID total_edges;
    bool write_log;
    LongNodeID stream_assigned_nodes;
    LongNodeID stream_n_nodes;
    LongNodeID nmbNodes;
    std::vector<std::vector<EdgeWeight>>* degree_nodeBlock;
    LongNodeID stream_total_upperbound;
    double fennel_gamma;
    double fennel_alpha;
    double fennel_alpha_gamma;
    bool use_fennel_objective;  // maps global blocks to current stream blocks
    int fennel_dynamics;
    bool ram_stream;
    int quotient_nodes;
    int lhs_nodes;
    bool stream_initial_bisections;
    int n_batches;
    double stream_global_epsilon;
    bool stream_output_progress;
    double batch_inbalance;
    bool skip_outer_ls;

    // Initial partition via growing multiple BFS trees
    bool initial_part_multi_bfs;
    int multibfs_tries;

    // Initial partitioning via Fennel on the coarsest level
    int initial_part_fennel_tries;

    // Ghost neighbors
    bool stream_allow_ghostnodes;
    LongNodeID ghost_nodes;
    LongNodeID ghost_nodes_threshold;
    bool double_non_ghost_edges;

    // Restreaming and partial restreaming
    int num_streams_passes;
    bool restream_vcycle;

    int xxx;
    double* t1;
    double* t2;
    double* t3;

    //=======================================
    //======= Stream Edge Partition ========
    //=======================================
    LongNodeID total_stream_edges;
    bool edge_partition;
    bool benchmark;
    bool dynamic_alpha;
    bool batch_alpha;
    bool minimal_mode;
    bool light_evaluator;
    bool use_queue;
    bool async_mode;
    bool evaluate_mode;
    bool evaluate;
    int past_subset_size;
    int quotient_edges_count;
    NodeID num_split_edges;
    LongEdgeID fennel_edges;
    LongNodeID lower_global_store;
    double reps;
    //=======================================
    //========== SpMxV partitioning =========
    //=======================================

    int matrix_m;
    int matrix_n;
    int matrix_nnz;

    //=======================================
    //===============Shared Mem OMP==========
    //=======================================
    bool enable_omp;

    // Grouped stream-runtime views to reduce direct coupling to the full config.
    struct StreamInputConfig {
        std::string& graph_filename;
        bool& ram_stream;
        LongNodeID& batch_size;
        int& num_streams_passes;
    };
    struct ConstStreamInputConfig {
        const std::string& graph_filename;
        const bool& ram_stream;
        const LongNodeID& batch_size;
        const int& num_streams_passes;
    };

    struct StreamExecutionConfig {
        bool& evaluate;
        bool& benchmark;
        bool& run_parallel;
    };
    struct ConstStreamExecutionConfig {
        const bool& evaluate;
        const bool& benchmark;
        const bool& run_parallel;
    };

    // Generic stream controls shared by all algorithms.
    struct CommonConfig {
        int& seed;
        bool& evaluate;
        bool& benchmark;
        bool& suppress_output;
        bool& write_log;
        std::string& graph_filename;
        LongNodeID& batch_size;
        bool& ram_stream;
        int& num_streams_passes;
    };
    struct ConstCommonConfig {
        const int& seed;
        const bool& evaluate;
        const bool& benchmark;
        const bool& suppress_output;
        const bool& write_log;
        const std::string& graph_filename;
        const LongNodeID& batch_size;
        const bool& ram_stream;
        const int& num_streams_passes;
    };

    struct StreamOutputConfig {
        std::string& filename_output;
        bool& stream_output_progress;
        bool& write_log;
    };
    struct ConstStreamOutputConfig {
        const std::string& filename_output;
        const bool& stream_output_progress;
        const bool& write_log;
    };

    // Partitioning-domain settings shared by node/edge partitioning algorithms.
    struct PartitionAlgorithmConfig {
        PartitionID& k;
        double& imbalance;
        bool& balance_edges;
        InitialPartitioningType& initial_partitioning_type;
        RefinementType& refinement_type;
        int& label_iterations;
        int& label_iterations_refinement;
        bool& use_fennel_objective;
        int& fennel_dynamics;
    };
    struct ConstPartitionAlgorithmConfig {
        const PartitionID& k;
        const double& imbalance;
        const bool& balance_edges;
        const InitialPartitioningType& initial_partitioning_type;
        const RefinementType& refinement_type;
        const int& label_iterations;
        const int& label_iterations_refinement;
        const bool& use_fennel_objective;
        const int& fennel_dynamics;
    };

    struct NodePartitioningTuningConfig {
        LongNodeID& max_buffer_size;
        unsigned& bq_disc_factor;
        BufferScoreType& buffer_score_type;
        LongNodeID& d_max;
        float& haa_beta;
        float& haa_theta;
        float& cbs_theta;
        size_t& max_input_q_size;
        bool& restream_include_high_degree_nodes;
    };
    struct ConstNodePartitioningTuningConfig {
        const LongNodeID& max_buffer_size;
        const unsigned& bq_disc_factor;
        const BufferScoreType& buffer_score_type;
        const LongNodeID& d_max;
        const float& haa_beta;
        const float& haa_theta;
        const float& cbs_theta;
        const size_t& max_input_q_size;
        const bool& restream_include_high_degree_nodes;
    };

    struct EdgePartitioningTuningConfig {
        bool& minimal_mode;
        bool& dynamic_alpha;
        bool& batch_alpha;
        bool& use_queue;
        int& past_subset_size;
        NodeID& num_split_edges;
    };
    struct ConstEdgePartitioningTuningConfig {
        const bool& minimal_mode;
        const bool& dynamic_alpha;
        const bool& batch_alpha;
        const bool& use_queue;
        const int& past_subset_size;
        const NodeID& num_split_edges;
    };

    StreamInputConfig stream_input_config() {
        return StreamInputConfig{graph_filename, ram_stream, batch_size, num_streams_passes};
    }
    ConstStreamInputConfig stream_input_config() const {
        return ConstStreamInputConfig{graph_filename, ram_stream, batch_size, num_streams_passes};
    }

    StreamExecutionConfig stream_execution_config() {
        return StreamExecutionConfig{evaluate, benchmark, run_parallel};
    }
    ConstStreamExecutionConfig stream_execution_config() const {
        return ConstStreamExecutionConfig{evaluate, benchmark, run_parallel};
    }
    CommonConfig common_config() {
        return CommonConfig{seed,           evaluate,   benchmark,  suppress_output,   write_log,
                            graph_filename, batch_size, ram_stream, num_streams_passes};
    }
    ConstCommonConfig common_config() const {
        return ConstCommonConfig{
            seed,           evaluate,   benchmark,  suppress_output,   write_log,
            graph_filename, batch_size, ram_stream, num_streams_passes};
    }

    StreamOutputConfig stream_output_config() {
        return StreamOutputConfig{filename_output, stream_output_progress, write_log};
    }
    ConstStreamOutputConfig stream_output_config() const {
        return ConstStreamOutputConfig{filename_output, stream_output_progress, write_log};
    }

    PartitionAlgorithmConfig partition_algorithm_config() {
        return PartitionAlgorithmConfig{k,
                                        imbalance,
                                        balance_edges,
                                        initial_partitioning_type,
                                        refinement_type,
                                        label_iterations,
                                        label_iterations_refinement,
                                        use_fennel_objective,
                                        fennel_dynamics};
    }
    ConstPartitionAlgorithmConfig partition_algorithm_config() const {
        return ConstPartitionAlgorithmConfig{k,
                                             imbalance,
                                             balance_edges,
                                             initial_partitioning_type,
                                             refinement_type,
                                             label_iterations,
                                             label_iterations_refinement,
                                             use_fennel_objective,
                                             fennel_dynamics};
    }
    // Backward-compatible alias while callers migrate.
    PartitionAlgorithmConfig partition_config() {
        return partition_algorithm_config();
    }
    ConstPartitionAlgorithmConfig partition_config() const {
        return partition_algorithm_config();
    }

    NodePartitioningTuningConfig node_tuning_config() {
        return NodePartitioningTuningConfig{
            max_buffer_size, bq_disc_factor,   buffer_score_type,
            d_max,           haa_beta,         haa_theta,
            cbs_theta,       max_input_q_size, restream_include_high_degree_nodes};
    }
    ConstNodePartitioningTuningConfig node_tuning_config() const {
        return ConstNodePartitioningTuningConfig{
            max_buffer_size, bq_disc_factor,   buffer_score_type,
            d_max,           haa_beta,         haa_theta,
            cbs_theta,       max_input_q_size, restream_include_high_degree_nodes};
    }

    EdgePartitioningTuningConfig edge_tuning_config() {
        return EdgePartitioningTuningConfig{minimal_mode, dynamic_alpha,    batch_alpha,
                                            use_queue,    past_subset_size, num_split_edges};
    }
    ConstEdgePartitioningTuningConfig edge_tuning_config() const {
        return ConstEdgePartitioningTuningConfig{minimal_mode, dynamic_alpha,    batch_alpha,
                                                 use_queue,    past_subset_size, num_split_edges};
    }

    void LogDump(FILE* out) const {}
};

#endif /* end of include guard: PARTITION_CONFIG_DI1ES4T0 */
