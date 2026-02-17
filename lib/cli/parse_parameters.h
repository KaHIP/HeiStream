/******************************************************************************
 * parse_parameters.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARSE_PARAMETERS_GPJMGSM8
#define PARSE_PARAMETERS_GPJMGSM8

#include <omp.h>
#include <regex.h>

#include <cstring>
#include <limits>
#include <sstream>

#include "algorithms/node_partitioning/mode/node_partitioning_mode_resolver.h"
#include "cli/configuration.h"
#include "cli/stream_parameter_apply.h"

int parse_parameters(int argn, char** argv, Config& partition_config, std::string& graph_filename,
                     bool& is_graph_weighted, bool& suppress_program_output, bool& recursive) {
    const char* progname = argv[0];

    // Setup argtable parameters.
    struct arg_lit* help = arg_lit0(nullptr, "help", "Print help.");
    struct arg_lit* edge_rating_tiebreaking =
        arg_lit0(nullptr, "edge_rating_tiebreaking", "Enable random edgerating tiebreaking.");
    struct arg_lit* match_islands =
        arg_lit0(nullptr, "match_islands", "Enable matching of islands during gpa algorithm.");
    struct arg_lit* only_first_level =
        arg_lit0(nullptr, "only_first_level",
                 "Disable Multilevel Approach. Only perform on the first level. (Currently only "
                 "initial partitioning).");
    struct arg_lit* graph_weighted =
        arg_lit0(nullptr, "weighted", "Read the graph as weighted graph.");
    struct arg_lit* enable_corner_refinement =
        arg_lit0(nullptr, "enable_corner_refinement", "Enables corner refinement.");
    struct arg_lit* disable_qgraph_refinement =
        arg_lit0(nullptr, "disable_qgraph_refinement", "Disables qgraph refinement.");
    struct arg_lit* use_fullmultigrid = arg_lit0(
        nullptr, "use_fullmultigrid", "Enable full multigrid (wcycles have to be enabled).");
    struct arg_lit* use_vcycle = arg_lit0(nullptr, "use_vcycle", "Enable vcycle .");
    struct arg_lit* compute_vertex_separator =
        arg_lit0(nullptr, "compute_vertex_separator", "Compute vertex separator.");
    struct arg_lit* first_level_random_matching = arg_lit0(
        nullptr, "first_level_random_matching", "The first level will be matched randomly.");
    struct arg_lit* rate_first_level_inner_outer =
        arg_lit0(nullptr, "rate_first_level_inner_outer",
                 "The edge rating for the first level is inner outer.");
    struct arg_lit* use_bucket_queues =
        arg_lit0(nullptr, "use_bucket_queues", "Use bucket priority queues during refinement.");
    struct arg_lit* use_wcycles = arg_lit0(nullptr, "use_wcycle", "Enables wcycles.");
    struct arg_lit* disable_refined_bubbling = arg_lit0(
        nullptr, "disable_refined_bubbling",
        "Disables refinement during initial partitioning using bubbling (Default: enabled).");
    struct arg_lit* enable_convergence = arg_lit0(nullptr, "enable_convergence",
                                                  "Enables convergence mode, i.e. every step is "
                                                  "running until no change.(Default: disabled).");
    struct arg_lit* enable_omp = arg_lit0(nullptr, "enable_omp", "Enable parallel omp.");
    struct arg_lit* wcycle_no_new_initial_partitioning =
        arg_lit0(nullptr, "wcycle_no_new_initial_partitioning",
                 "Using this option, the graph is initially partitioned only the first time we are "
                 "at the deepest level.");
    struct arg_str* filename =
        arg_strn(nullptr, nullptr, "FILE", 1, 1, "Path to graph file to partition.");
    struct arg_str* filename_output =
        arg_str0(nullptr, "output_filename", nullptr,
                 "Specify the name of the output file (that contains the partition).");
    struct arg_int* user_seed = arg_int0(nullptr, "seed", nullptr, "Seed to use for the PRNG.");
    struct arg_int* k = arg_int1(nullptr, "k", nullptr, "Number of blocks to partition the graph.");
    struct arg_rex* edge_rating =
        arg_rex0(nullptr, "edge_rating",
                 "^(weight|realweight|expansionstar|expansionstar2|expansionstar2deg|punch|"
                 "expansionstar2algdist|expansionstar2algdist2|algdist|algdist2|sepmultx|sepaddx|"
                 "sepmax|seplog|r1|r2|r3|r4|r5|r6|r7|r8)$",
                 "RATING", REG_EXTENDED,
                 "Edge rating to use. One of {weight, expansionstar, expansionstar2, punch, "
                 "sepmultx, sepaddx, sepmax, seplog, "
                 " expansionstar2deg}. Default: weight");
    struct arg_rex* refinement_type =
        arg_rex0(nullptr, "refinement_type", "^(fm|fm_flow|flow)$", "TYPE", REG_EXTENDED,
                 "Refinementvariant to use. One of {fm, fm_flow, flow}. Default: fm");
    struct arg_rex* matching_type =
        arg_rex0(nullptr, "matching", "^(cluster|random|hem|shem|regions|gpa|randomgpa|localmax)$",
                 "TYPE", REG_EXTENDED,
                 "Type of matchings to use during coarsening. One of {random, hem,"
                 " shem, regions, gpa, randomgpa, localmax}.");
    struct arg_int* mh_pool_size =
        arg_int0(nullptr, "mh_pool_size", nullptr, "MetaHeuristic Pool Size.");
    struct arg_lit* mh_plain_repetitions = arg_lit0(nullptr, "mh_plain_repetitions", "");
    struct arg_lit* mh_penalty_for_unconnected =
        arg_lit0(nullptr, "mh_penalty_for_unconnected",
                 "Add a penalty on the objective function if the computed partition contains "
                 "blocks that are not connected.");
    struct arg_lit* mh_disable_nc_combine = arg_lit0(nullptr, "mh_disable_nc_combine", "");
    struct arg_lit* mh_disable_cross_combine = arg_lit0(nullptr, "mh_disable_cross_combine", "");
    struct arg_lit* mh_disable_combine = arg_lit0(nullptr, "mh_disable_combine", "");
    struct arg_lit* mh_enable_quickstart =
        arg_lit0(nullptr, "mh_enable_quickstart", "Enables the quickstart option.");
    struct arg_lit* mh_disable_diversify_islands =
        arg_lit0(nullptr, "mh_disable_diversify_islands", "");
    struct arg_lit* mh_disable_diversify = arg_lit0(nullptr, "mh_disable_diversify", "");
    struct arg_lit* mh_diversify_best =
        arg_lit0(nullptr, "mh_diversify_best",
                 "Uses best individuum instead of random during diversification.");
    struct arg_lit* mh_enable_tournament_selection = arg_lit0(
        nullptr, "mh_enable_tournament_selection",
        "Enables the tournament selection roule instead of choosing two random inidiviuums.");
    struct arg_lit* mh_cross_combine_original_k =
        arg_lit0(nullptr, "mh_cross_combine_original_k", "");
    struct arg_lit* mh_optimize_communication_volume =
        arg_lit0(nullptr, "mh_optimize_communication_volume",
                 "Fitness function is modified to optimize communication volume instead of the "
                 "number of cut edges.");
    struct arg_lit* disable_balance_singletons =
        arg_lit0(nullptr, "disable_balance_singletons", "");
    struct arg_lit* gpa_grow_internal =
        arg_lit0(nullptr, "gpa_grow_internal",
                 "If the graph is allready partitions the paths are grown only block internally.");
    struct arg_int* initial_partitioning_repetitions =
        arg_int0(nullptr, "initial_partitioning_repetitions", nullptr,
                 "Number of initial partitioning repetitons. Default: 5.");
    struct arg_int* minipreps = arg_int0(nullptr, "minipreps", nullptr, "Default: 10.");
    struct arg_int* aggressive_random_levels =
        arg_int0(nullptr, "aggressive_random_levels", nullptr,
                 "In case matching is randomgpa, this is the number of levels that should be "
                 "matched using random matching. Default: 3.");
    struct arg_dbl* imbalance =
        arg_dbl0(nullptr, "imbalance", nullptr, "Desired balance. Default: 3 (%).");
    struct arg_rex* initial_partition =
        arg_rex0(nullptr, "initial_partitioner",
                 "^(metis|scotch|hybrid|bubbling|squeez|metaheuristic|recursive)$", "PARTITIONER",
                 REG_EXTENDED,
                 "Type of matchings to use during coarsening. One of {metis, scotch, bubbling, "
                 "hybrid, recursive).");
    struct arg_lit* initial_partition_optimize = arg_lit0(
        nullptr, "initial_partition_optimize", "Enables postoptimization of initial partition.");
    struct arg_rex* bipartition_algorithm =
        arg_rex0(nullptr, "bipartition_algorithm", "^(bfs|fm|squeezing)$", "TYPE", REG_EXTENDED,
                 "Type of bipartition algorithm to use in case of recursive partitioning. One of "
                 " {bfs, fm, squeezing}.");
    struct arg_rex* permutation_quality =
        arg_rex0(nullptr, "permutation_quality", "^(none|fast|good|cacheefficient)$", "QUALITY",
                 REG_EXTENDED,
                 "The quality of permutations to use. One of {none, fast,"
                 " good, cacheefficient}.");
    struct arg_rex* permutation_during_refinement = arg_rex0(
        nullptr, "permutation_during_refinement", "^(none|fast|good)$", "QUALITY", REG_EXTENDED,
        "The quality of permutations to use during 2way fm refinement. One of {none, fast,"
        " good}.");
    struct arg_int* fm_search_limit =
        arg_int0(nullptr, "fm_search_limit", nullptr,
                 "Search limit for 2way fm local search: Default 1 (%).");
    struct arg_int* bipartition_post_fm_limit =
        arg_int0(nullptr, "bipartition_post_fm_limit", nullptr,
                 "Search limit for the fm search after a bipartition has been created. :");
    struct arg_int* bipartition_post_ml_limit = arg_int0(
        nullptr, "bipartition_post_ml_limit", nullptr,
        "Search limit for the multilevel fm search after a bipartition has been created. :");
    struct arg_int* bipartition_tries =
        arg_int0(nullptr, "bipartition_tries", nullptr,
                 "Number of tries to find a bipartition (during recursive intial partitioning).");
    struct arg_rex* refinement_scheduling_algorithm = arg_rex0(
        nullptr, "refinement_scheduling_algorithm", "^(fast|active_blocks|active_blocks_kway)$",
        "QUALITY", REG_EXTENDED, " One of {fast, active_blocks, active_blocks_kway}.");
    struct arg_dbl* bank_account_factor =
        arg_dbl0(nullptr, "bank_account_factor", nullptr,
                 "The bank account factor for the scheduler. Default 1.5 (%).");
    struct arg_dbl* flow_region_factor =
        arg_dbl0(nullptr, "flow_region_factor", nullptr,
                 "If using flow, then the regions found are sized flow_region_factor * imbalance. "
                 "Default: 4 (%).");
    struct arg_dbl* kway_adaptive_limits_alpha =
        arg_dbl0(nullptr, "kway_adaptive_limits_alpha", nullptr,
                 "This is the factor alpha used for the adaptive stopping criteria. Default: 1.0");
    struct arg_rex* stop_rule =
        arg_rex0(nullptr, "stop_rule", "^(simple|multiplek|strong)$", "VARIANT", REG_EXTENDED,
                 "Stop rule to use. One of {simple, multiplek, strong}. Default: simple");
    struct arg_int* num_vert_stop_factor = arg_int0(nullptr, "num_vert_stop_factor", nullptr,
                                                    "x*k (for multiple_k stop rule). Default 20.");
    struct arg_rex* kway_search_stop_rule = arg_rex0(
        nullptr, "kway_stop_rule", "^(simple|adaptive)$", "VARIANT", REG_EXTENDED,
        "Stop rule to use during kway_refinement. One of {simple, adaptive}. Default: simple");
    struct arg_int* bubbling_iterations =
        arg_int0(nullptr, "bubbling_iterations", nullptr,
                 "Number of bubbling iterations to perform: Default 1 .");
    struct arg_int* kway_rounds =
        arg_int0(nullptr, "kway_rounds", nullptr,
                 "Number of kway refinement rounds to perform: Default 1 .");
    struct arg_int* kway_fm_limits = arg_int0(nullptr, "kway_fm_search_limit", nullptr,
                                              "Search limit for kway fm local search: Default 1 .");
    struct arg_int* global_cycle_iterations = arg_int0(nullptr, "global_cycle_iterations", nullptr,
                                                       "Number of V-cycle iterations: Default 2.");
    struct arg_int* level_split = arg_int0(nullptr, "level_split", nullptr,
                                           "Number of trial tree levels (1 means on each level a "
                                           "two trials are performed). Default: 2.");
    struct arg_int* toposort_iterations = arg_int0(nullptr, "toposort_iterations", nullptr,
                                                   "Number of topo sort iterations). Default: 4.");
    struct arg_lit* most_balanced_flows =
        arg_lit0(nullptr, "most_balanced_flows", "(Default: disabled)");
    struct arg_str* input_partition =
        arg_str0(nullptr, "input_partition", nullptr, "Input partition to use.");
    struct arg_lit* recursive_bipartitioning =
        arg_lit0(nullptr, "recursive_bipartitioning",
                 "Use recursive bipartitioning instead of kway methods.");
    struct arg_lit* suppress_output =
        arg_lit0(nullptr, "suppress_output", "(Default: output enabled)");
    struct arg_lit* disable_max_vertex_weight_constraint =
        arg_lit0(nullptr, "disable_max_vertex_weight_constraint",
                 "Disables the max vertex weight constraint during the contraction.");
    struct arg_int* local_multitry_fm_alpha = arg_int0(
        nullptr, "local_multitry_fm_alpha", nullptr, "Search limit factor alpha for multitry fm.");
    struct arg_int* local_multitry_rounds = arg_int0(nullptr, "local_multitry_rounds", nullptr,
                                                     "Number of rounds for local multitry fm.");
    struct arg_int* initial_partition_optimize_fm_limits =
        arg_int0(nullptr, "initial_partition_optimize_fm_limits", nullptr,
                 "Initial Partition Optimize FM limits. (Default: 20)");
    struct arg_int* initial_partition_optimize_multitry_fm_alpha =
        arg_int0(nullptr, "initial_partition_optimize_multitry_fm_limits", nullptr,
                 "Initial Partition Optimize Multitry FM limits. (Default: 20)");
    struct arg_int* initial_partition_optimize_multitry_rounds =
        arg_int0(nullptr, "initial_partition_optimize_multitry_rounds", nullptr, "(Default: 100)");

#ifdef MODE_KAFFPA
    struct arg_rex* preconfiguration = arg_rex1(
        nullptr, "preconfiguration", "^(strong|eco|fast|fsocial|esocial|ssocial)$", "VARIANT",
        REG_EXTENDED,
        "Use a preconfiguration. (Default: eco) [strong|eco|fast|fsocial|esocial|ssocial].");
#else
    struct arg_rex* preconfiguration = arg_rex0(
        nullptr, "preconfiguration", "^(strong|eco|fast|fsocial|esocial|ssocial)$", "VARIANT",
        REG_EXTENDED,
        "Use a preconfiguration. (Default: strong) [strong|eco|fast|fsocial|esocial|ssocial].");
#endif

    struct arg_dbl* time_limit =
        arg_dbl0(nullptr, "time_limit", nullptr, "Time limit in s. Default 0s .");
    struct arg_int* unsuccessful_reps =
        arg_int0(nullptr, "unsuccessful_reps", nullptr, "Unsuccessful reps to fresh start.");
    struct arg_int* local_partitioning_repetitions = arg_int0(
        nullptr, "local_partitioning_repetitions", nullptr, "Number of local repetitions.");
    struct arg_int* amg_iterations =
        arg_int0(nullptr, "amg_iterations", nullptr, "Number of amg iterations.");
    struct arg_int* mh_flip_coin = arg_int0(
        nullptr, "mh_flip_coin", nullptr,
        "Control the ratio of mutation and crossovers. c/10 Mutation and (10-c)/10 crossovers.");
    struct arg_int* mh_initial_population_fraction =
        arg_int0(nullptr, "mh_initial_population_fraction", nullptr,
                 "Control the initial population fraction parameter (Default: 1000).");
    struct arg_lit* mh_print_log =
        arg_lit0(nullptr, "mh_print_log", "Each PE prints a logfile (timestamp, edgecut).");
    struct arg_lit* mh_sequential_mode =
        arg_lit0(nullptr, "mh_sequential_mode",
                 "Disables all MH algorithms. Use KaFFPa in a parallel setting.");
    struct arg_rex* kaba_neg_cycle_algorithm =
        arg_rex0(nullptr, "kaba_neg_cycle_algorithm",
                 "^(ultramodel|randomcycle|playfield|ultramodelplus)$", "VARIANT", REG_EXTENDED,
                 "Balanced refinement operator to use. On of randomcycle, ultramodel, playfield, "
                 "ultramodelplus");
    struct arg_dbl* kabaE_internal_bal =
        arg_dbl0(nullptr, "kabaE_internal_bal", nullptr,
                 "Control the internal balance paramter for kaffpaE (Default: 0.01) (1 percent)");
    struct arg_int* kaba_internal_no_aug_steps_aug =
        arg_int0(nullptr, "kaba_internal_no_aug_steps_aug", nullptr,
                 "Internal number of steps in the augmented models of negative cycle detection.");
    struct arg_int* kaba_packing_iterations =
        arg_int0(nullptr, "kaba_packing_iterations", nullptr, "Number of packing iterations.");
    struct arg_int* kaba_unsucc_iterations =
        arg_int0(nullptr, "kaba_unsucc_iterations", nullptr,
                 "Number of unsucc iterations until a rebalancing step is performed.");
    struct arg_lit* kaba_flip_packings = arg_lit0(
        nullptr, "kaba_flip_packings", "Enable flip packing mode (if ultramodelplus is used).");
    struct arg_rex* kaba_lsearch_p =
        arg_rex0(nullptr, "kaba_lsearch_p", "^(coindiff|coinrnd|nocoindiff|nocoinrnd)$", "VARIANT",
                 REG_EXTENDED, "Make more localized search in ultraplus model.");
    struct arg_lit* kaffpa_perfectly_balanced_refinement =
        arg_lit0(nullptr, "kaffpa_perfectly_balanced_refinement",
                 "Enable perfectly balanced refinement during ML KaFFPa.");
    struct arg_lit* kaba_disable_zero_weight_cycles = arg_lit0(
        nullptr, "kaba_disable_zero_weight_cycles", "Disable zero weight cycle diversification.");
    struct arg_lit* enforce_balance =
        arg_lit0(nullptr, "enforce_balance",
                 "Uses eps+1 to create a partition, and then runs rebalancing and negative cycle "
                 "detection to output a partition that fulfills the eps-balance constraint.");
    struct arg_lit* mh_enable_tabu_search =
        arg_lit0(nullptr, "mh_enable_tabu_search",
                 "Enables our version of combine operation by block matching; +tabusearch and all "
                 "our refinement algorithms.");
    struct arg_lit* mh_enable_kabapE =
        arg_lit0(nullptr, "mh_enable_kabapE", "Enable combine operator KaBaPE");
    struct arg_int* maxT = arg_int0(nullptr, "maxT", nullptr, "maxT parameter for Tabu Search");
    struct arg_int* maxIter =
        arg_int0(nullptr, "maxIter", nullptr, "maxIter parameter for Tabu Search");
    struct arg_lit* balance_edges =
        arg_lit0(nullptr, "balance_edges", "Turn on balancing of edges among blocks.");

    struct arg_int* cluster_upperbound =
        arg_int0(nullptr, "cluster_upperbound", nullptr,
                 "Set a size-constraint on the size of a cluster. Default: none");
    struct arg_int* label_propagation_iterations =
        arg_int0(nullptr, "label_propagation_iterations", nullptr,
                 "Set the number of label propgation iterations. Default: 10.");
    struct arg_int* qap_label_iterations =
        arg_int0(nullptr, "qap_label_iterations", nullptr,
                 "Set the number of label propgation iterations. Default: 10.");

    struct arg_int* max_initial_ns_tries =
        arg_int0(nullptr, "max_initial_ns_tries", nullptr,
                 "Number of NS tries during initial partitioning.");
    struct arg_int* max_flow_improv_steps =
        arg_int0(nullptr, "max_flow_improv_steps", nullptr,
                 "Maximum number of tries to improve a node separator using flows.");
    struct arg_lit* most_balanced_flows_node_sep =
        arg_lit0(nullptr, "most_balanced_flows_node_sep", "(Default: disabled)");
    struct arg_dbl* region_factor_node_separators =
        arg_dbl0(nullptr, "region_factor_node_separators", nullptr,
                 "Region factor for flow problems to obtain node separators.");
    struct arg_lit* sep_flows_disabled =
        arg_lit0(nullptr, "sep_flows_disabled", "(Default: disabled)");
    struct arg_lit* sep_fm_disabled = arg_lit0(nullptr, "sep_fm_disabled", "(Default: disabled)");
    struct arg_lit* sep_loc_fm_disabled =
        arg_lit0(nullptr, "sep_loc_fm_disabled", "(Default: disabled)");
    struct arg_lit* sep_greedy_disabled =
        arg_lit0(nullptr, "sep_greedy_disabled", "(Default: disabled)");
    struct arg_lit* sep_full_boundary_ip =
        arg_lit0(nullptr, "sep_full_boundary_ip", "(Default: disabled)");
    struct arg_lit* sep_faster_ns = arg_lit0(nullptr, "sep_faster_ns", "(Default: disabled)");
    struct arg_int* sep_fm_unsucc_steps =
        arg_int0(nullptr, "sep_fm_unsucc_steps", nullptr,
                 "Maximum number of steps till last improvement in FM algorithm.");
    struct arg_int* sep_num_fm_reps =
        arg_int0(nullptr, "sep_num_fm_reps", nullptr,
                 "Number of FM repetitions during uncoarsening on each level.");
    struct arg_int* sep_loc_fm_unsucc_steps =
        arg_int0(nullptr, "sep_loc_fm_unsucc_steps", nullptr,
                 "Maximum number of steps till last improvement in FM algorithm.");
    struct arg_int* sep_num_loc_fm_reps =
        arg_int0(nullptr, "sep_num_loc_fm_reps", nullptr,
                 "Number of FM repetitions during uncoarsening on each level.");
    struct arg_int* sep_loc_fm_no_snodes =
        arg_int0(nullptr, "sep_loc_fm_no_snodes", nullptr,
                 "Number of FM repetitions during uncoarsening on each level.");
    struct arg_int* sep_num_vert_stop = arg_int0(nullptr, "sep_num_vert_stop", nullptr,
                                                 "Number of vertices to stop coarsening at.");
    struct arg_rex* sep_edge_rating_during_ip =
        arg_rex0(nullptr, "sep_edge_rating_during_ip",
                 "^(weight|expansionstar|expansionstar2|expansionstar2deg|punch|"
                 "expansionstar2algdist|expansionstar2algdist2|algdist|algdist2|sepmultx|sepaddx|"
                 "sepmax|seplog|r1|r2|r3|r4|r5|r6|r7|r8)$",
                 "RATING", REG_EXTENDED,
                 "Edge rating to use. One of {weight, expansionstar, expansionstar2, punch, "
                 "sepmultx, sepaddx, sepmax, seplog, "
                 " expansionstar2deg}. Default: weight");

    // mapping stuff
    //
    //
    struct arg_lit* integrated_mapping =
        arg_lit0(nullptr, "integrated_mapping",
                 "Enable integrated mapping algorithms to map quotient graph onto processor graph "
                 "defined by hierarchy and distance options. (Default: disabled)");
    struct arg_lit* multisection = arg_lit0(nullptr, "multisection",
                                            "Enable multisectioning through the hierarchy instead "
                                            "of bipartitioning. (Default: disabled)");
    struct arg_lit* qap_label_propagation =
        arg_lit0(nullptr, "qap_label_propagation",
                 "Enable mapping label propagation local search to work with integrated mapping. "
                 "(Default: disabled)");
    struct arg_lit* qap_blabel_propagation =
        arg_lit0(nullptr, "qap_blabel_propagation",
                 "Enable mapping (before all) label propagation local search to work with "
                 "integrated mapping. (Default: disabled)");
    struct arg_lit* qap_alabel_propagation =
        arg_lit0(nullptr, "qap_alabel_propagation",
                 "Enable mapping (after all) label propagation local search to work with "
                 "integrated mapping. (Default: disabled)");
    struct arg_lit* qap_multitry_fm =
        arg_lit0(nullptr, "qap_multitry_fm",
                 "Enable mapping localized local search (mapping multitry fm) to work with "
                 "integrated mapping. (Default: disabled)");
    struct arg_lit* qap_bmultitry_fm =
        arg_lit0(nullptr, "qap_bmultitry_fm",
                 "Enable mapping multitry fm before mapping label propagation for integrated "
                 "mapping. (Default: disabled)");
    struct arg_lit* qap_kway_fm =
        arg_lit0(nullptr, "qap_kway_fm",
                 "Enable mapping kway fm to work with integrated mapping. (Default: disabled)");
    struct arg_lit* qap_bkway_fm = arg_lit0(nullptr, "qap_bkway_fm",
                                            "Enable mapping kway fm to work with integrated "
                                            "mapping before bmultitry. (Default: disabled)");
    struct arg_lit* qap_quotient_ref = arg_lit0(nullptr, "qap_quotient_ref",
                                                "Enable mapping quotient graph refinement to work "
                                                "with integrated mapping. (Default: disabled)");
    struct arg_lit* qap_bquotient_ref =
        arg_lit0(nullptr, "qap_bquotient_ref",
                 "Enable mapping (before) quotient graph refinement to work with integrated "
                 "mapping. (Default: disabled)");
    struct arg_lit* qap_0quotient_ref =
        arg_lit0(nullptr, "qap_0quotient_ref",
                 "Enable mapping (first) quotient graph refinement to work with integrated "
                 "mapping. (Default: disabled)");
    struct arg_lit* quotient_more_mem =
        arg_lit0(nullptr, "quotient_more_mem",
                 "Enable mapping more memory storage during QAP quotient graph refinement step. "
                 "(Default: disabled)");
    struct arg_lit* disable_bipartition_gp_local_search = arg_lit0(
        nullptr, "disable_bipartition_gp_local_search",
        "Disable graph partitioning--based local searches during the bipartitioning for the "
        "initial partitioning. Only works with integrated mapping. (Default: enabled)");
    struct arg_lit* enable_mapping =
        arg_lit0(nullptr, "enable_mapping",
                 "Enable mapping algorithms to map quotient graph onto processor graph defined by "
                 "hierarchy and distance options. (Default: disabled)");
    struct arg_str* hierarchy_parameter_string =
        arg_str0(nullptr, "hierarchy_parameter_string", nullptr,
                 "Specify as 4:8:8 for 4 cores per PE, 8 PEs per rack, ... and so forth.");
    struct arg_str* distance_parameter_string =
        arg_str0(nullptr, "distance_parameter_string", nullptr,
                 "Specify as 1:10:100 if cores on the same chip have distance 1, PEs in the same "
                 "rack have distance 10, ... and so forth.");
    struct arg_lit* online_distances = arg_lit0(
        nullptr, "online_distances",
        "Do not store processor distances in a matrix, but do recomputation. (Default: disabled)");
    struct arg_str* map_construction_algorithm = arg_str0(
        nullptr, "map_construction_algorithm", nullptr, "Specify mapping construction algorithm.");
    struct arg_lit* skip_map_ls =
        arg_lit0(nullptr, "skip_map_ls", "Skip mapping local search. (Default: disabled)");
    struct arg_lit* delta_gains =
        arg_lit0(nullptr, "delta_gains",
                 "Use delta gains to accelerate map local searches. (Default: disabled)");
    struct arg_lit* use_bin_id =
        arg_lit0(nullptr, "use_bin_id",
                 "Hierarchically split IDs for faster online distances. (Default: disabled)");
    struct arg_lit* use_compact_bin_id = arg_lit0(
        nullptr, "use_compact_bin_id", "Use binary IDs for online distances. (Default: disabled)");
    struct arg_lit* full_matrix =
        arg_lit0(nullptr, "full_matrix", "Store full topology matrix. (Default: disabled)");
    struct arg_lit* enable_convergence_map =
        arg_lit0(nullptr, "enable_convergence_map",
                 "Enables convergence mode, i.e. every step is running until no change.(Default: "
                 "disabled).");
    struct arg_lit* adapt_bal = arg_lit0(nullptr, "adapt_bal",
                                         "Use adaptive balancing to improve balancing during "
                                         "inital construction. (Default: disabled)");

    // Stream Partition
    struct arg_int* stream_buffer = arg_int0(nullptr, "stream_buffer", nullptr,
                                             "Legacy alias for --batch_size (nodes per batch).");
    struct arg_lit* run_parallel = arg_lit0(
        nullptr, "run-parallel",
        "Run the parallel priority-buffer pipeline for node streaming. (Default: disabled).");
    struct arg_int* max_buffer_size =
        arg_int0(nullptr, "buffer_size", nullptr,
                 "Maximum bucket priority queue (buffer) size. Set >1 to enable priority "
                 "buffering. Default 0 (disabled).");
    struct arg_int* batch_size = arg_int0(nullptr, "batch_size", nullptr,
                                          "Batch size (nodes per MLP batch): Default 16384.");
    struct arg_int* bq_disc_factor = arg_int0(nullptr, "bq_disc_factor", nullptr,
                                              "Discretization factor for bucket pq: Default 1000.");
    struct arg_str* buffer_score =
        arg_str0(nullptr, "b_score", nullptr,
                 "Buffer score used in PQ (haa|cbs|cbsq|anr|cms|nss|gts). Default: haa.");
    struct arg_dbl* haa_beta =
        arg_dbl0(nullptr, "haa_beta", nullptr, "Beta parameter for HAA buffer score. Default 2.0.");
    struct arg_dbl* haa_theta =
        arg_dbl0(nullptr, "theta", nullptr, "Theta parameter for HAA buffer score. Default 0.75.");
    struct arg_lit* ghost_neighbors_enabled =
        arg_lit0(nullptr, "enable_ghost", "Enable ghost nodes in MLP. (Default: disabled).");
    struct arg_lit* sep_batch_marker =
        arg_lit0(nullptr, "sep_batch_marker",
                 "Use separate batch marker for batch id marking. (Default: disabled).");
    struct arg_lit* restream_include_high_degree_nodes =
        arg_lit0(nullptr, "include_highdeg_nodes",
                 "Include high degree nodes in MLP when restreaming. (Default: disabled).");
    struct arg_int* d_max = arg_int0(nullptr, "d_max", nullptr,
                                     "Maximum degree to be inserted into queue. Default 10000.");
    struct arg_int* bb_ratio =
        arg_int0(nullptr, "bb_ratio", nullptr, "Buffer to batch size ratio. Default: not active.");
    struct arg_lit* use_fennel_objective =
        arg_lit0(nullptr, "use_fennel_objective",
                 "Use Fennel objetcive function in the local search. (Default: disabled)");
    struct arg_str* fennel_dynamics =
        arg_str0(nullptr, "fennel_dynamics", nullptr,
                 "Dynamic behavior of Fennel objective in local search "
                 "(original|double|linear|quadratic|midlinear|midquadratic|midconstant|edgecut). "
                 "(Default: original)");
    struct arg_lit* ram_stream =
        arg_lit0(nullptr, "ram_stream", "Stream from RAM instead of HD. (Default: disabled)");
    struct arg_lit* write_log =
        arg_lit0(nullptr, "write_log",
                 "Log experimental output in a flatbuffer .bin file. (Default: disabled)");
    struct arg_str* evaluate =
        arg_str0(nullptr, "evaluate", nullptr,
                 "Evaluate final partition and print run summaries (true|false). Default: true.");
    struct arg_lit* stream_initial_bisections =
        arg_lit0(nullptr, "stream_initial_bisections",
                 "Compute initial solution at the coarsest level for each buffer. (Default: "
                 "compute preliminary Fennel)");
    struct arg_lit* stream_output_progress =
        arg_lit0(nullptr, "stream_output_progress",
                 "Output global partition after each batch is partitioned. (Default: disabled)");
    struct arg_lit* stream_allow_ghostnodes =
        arg_lit0(nullptr, "stream_allow_ghostnodes",
                 "Improve partitioning using information about unvisited neighbors of the nodes in "
                 "each batch. (Default: disabled)");
    struct arg_str* ghost_nodes_procedure =
        arg_str0(nullptr, "ghost_nodes_procedure", nullptr,
                 "procedure regarding ghost neighbors (contract|keep|keepthresholdcontract). "
                 "(Default: contract)");
    struct arg_dbl* ghost_nodes_threshold =
        arg_dbl0(nullptr, "ghost_nodes_threshold", nullptr,
                 "Number of uncontracted ghost nodes allowed per batch. (Default: 1024).");
    struct arg_dbl* num_streams_passes =
        arg_dbl0(nullptr, "num_streams_passes", nullptr,
                 "Number of times input graph should be streamed. (Default: 1).");
    struct arg_lit* restream_vcycle =
        arg_lit0(nullptr, "restream_vcycle",
                 "Don't recompute recursive bisections in restreamed graphs. (Default: disabled)");
    struct arg_dbl* batch_inbalance =
        arg_dbl0(nullptr, "batch_inbalance", nullptr,
                 "Relative inbalance allowed in a batch. Default: 20 (%).");
    struct arg_lit* initial_part_multi_bfs =
        arg_lit0(nullptr, "initial_part_multi_bfs",
                 "Initial partition w/ multiple BFS trees initialized w/ artificial nodes. "
                 "(Default: initial partition via recursive bisection)");
    struct arg_lit* initial_part_fennel = arg_lit0(
        nullptr, "initial_part_fennel",
        "Initial partition w/ Fennel. (Default: initial partition via recursive bisection)");
    struct arg_lit* skip_outer_ls = arg_lit0(
        nullptr, "skip_outer_ls", "Skip outer local search procedures. (Default: disabled)");
    struct arg_str* stream_label_rounds =
        arg_str0(nullptr, "stream_label_rounds", nullptr,
                 "Amount of label propagation rounds (minimal|normal|high). (Default: normal)");
    struct arg_lit* automatic_buffer_len =
        arg_lit0(nullptr, "automatic_buffer_len",
                 "Automatically choose buffer size for fastest performance. (Default: disabled)");
    struct arg_int* xxx =
        arg_int0(nullptr, "xxx", nullptr, "tuning factor for size of coarsest graph. Default 4.");

    // Stream Edge Partition
    struct arg_lit* benchmark = arg_lit0(nullptr, "benchmark",
                                         "Do not output partition IDs to benchmark time and memory "
                                         "consumption. (Default: disabled)");
    struct arg_str* use_queue = arg_str0(
        nullptr, "use_queue", nullptr,
        "Use queue candidate mode for Fennel (true|false|1|0). (Default: true)");
    struct arg_lit* dynamic_alpha =
        arg_lit0(nullptr, "dynamic_alpha", "Dynamically update fennel alpha. (Default: disabled)");
    struct arg_lit* batch_alpha =
        arg_lit0(nullptr, "batch_alpha",
                 "Dynamically update fennel alpha based on batch nodes and "
                 "edges. (Default: disabled)");
    struct arg_lit* minimal_mode =
        arg_lit0(nullptr, "minimal_mode",
                 "Build graph model with only one latest quotient block per "
                 "node. (Default: disabled)");
    struct arg_lit* light_evaluator =
        arg_lit0(nullptr, "light_evaluator",
                 "Use light evaluator for huge graphs at small k. (Default: disabled)");
    struct arg_int* num_split_edges =
        arg_int0(nullptr, "num_split_edges", nullptr,
                 "Provide number of edges of entire split graph. (Default: -1)");
    struct arg_int* past_subset_size =
        arg_int0(nullptr, "past_subset_size", nullptr,
                 "Number of random past blocks to connect to in graph model. "
                 "(Default: -1 = All)");
    struct arg_dbl* reps = arg_dbl0(nullptr, "reps", nullptr,
                                    "Number of times batches should be coarsened. (Default: 1).");
    struct arg_int* label_propagation_iterations_refinement =
        arg_int0(nullptr, "label_propagation_iterations_refinement", nullptr,
                 "Set the number of label propagation iterations for refinement. "
                 "Default: 5.");

    // translation of graphs

    struct arg_end* end = arg_end(100);

    // Define argtable.
    void* argtable[] = {help,
                        filename,
                        user_seed,
#ifdef MODE_DEVEL
                        k,
                        graph_weighted,
                        imbalance,
                        edge_rating_tiebreaking,
                        matching_type,
                        edge_rating,
                        rate_first_level_inner_outer,
                        first_level_random_matching,
                        aggressive_random_levels,
                        gpa_grow_internal,
                        match_islands,
                        stop_rule,
                        num_vert_stop_factor,
                        initial_partition,
                        initial_partitioning_repetitions,
                        disable_refined_bubbling,
                        bubbling_iterations,
                        initial_partition_optimize,
                        bipartition_post_fm_limit,
                        bipartition_post_ml_limit,
                        bipartition_tries,
                        bipartition_algorithm,
                        permutation_quality,
                        permutation_during_refinement,
                        enforce_balance,
                        refinement_scheduling_algorithm,
                        bank_account_factor,
                        refinement_type,
                        fm_search_limit,
                        flow_region_factor,
                        most_balanced_flows,
                        toposort_iterations,
                        kway_rounds,
                        kway_search_stop_rule,
                        kway_fm_limits,
                        kway_adaptive_limits_alpha,
                        enable_corner_refinement,
                        disable_qgraph_refinement,
                        local_multitry_fm_alpha,
                        local_multitry_rounds,
                        global_cycle_iterations,
                        use_wcycles,
                        wcycle_no_new_initial_partitioning,
                        use_fullmultigrid,
                        use_vcycle,
                        level_split,
                        enable_convergence,
                        compute_vertex_separator,
                        suppress_output,
                        input_partition,
                        preconfiguration,
                        only_first_level,
                        disable_max_vertex_weight_constraint,
                        recursive_bipartitioning,
                        use_bucket_queues,
                        time_limit,
                        unsuccessful_reps,
                        local_partitioning_repetitions,
                        mh_pool_size,
                        mh_plain_repetitions,
                        mh_disable_nc_combine,
                        mh_disable_cross_combine,
                        mh_enable_tournament_selection,
                        mh_disable_combine,
                        mh_enable_quickstart,
                        mh_disable_diversify_islands,
                        mh_flip_coin,
                        mh_initial_population_fraction,
                        mh_print_log,
                        mh_sequential_mode,
                        mh_optimize_communication_volume,
                        mh_enable_tabu_search,
                        mh_disable_diversify,
                        mh_diversify_best,
                        mh_cross_combine_original_k,
                        disable_balance_singletons,
                        initial_partition_optimize_fm_limits,
                        initial_partition_optimize_multitry_fm_alpha,
                        initial_partition_optimize_multitry_rounds,
                        enable_omp,
                        amg_iterations,
                        kaba_neg_cycle_algorithm,
                        kabaE_internal_bal,
                        kaba_internal_no_aug_steps_aug,
                        kaba_packing_iterations,
                        kaba_flip_packings,
                        kaba_lsearch_p,
                        kaffpa_perfectly_balanced_refinement,
                        kaba_unsucc_iterations,
                        kaba_disable_zero_weight_cycles,
                        maxT,
                        maxIter,
                        minipreps,
                        mh_penalty_for_unconnected,
                        mh_enable_kabapE,

#elif defined MODE_STREAMPARTITION_DEV
                        k,
                        imbalance,
                        preconfiguration,
                        time_limit,
                        enforce_balance,
                        balance_edges,
                        integrated_mapping,
                        multisection,
                        qap_label_propagation,
                        qap_blabel_propagation,
                        qap_alabel_propagation,
                        qap_multitry_fm,
                        qap_bmultitry_fm,
                        qap_kway_fm,
                        qap_bkway_fm,
                        qap_quotient_ref,
                        qap_bquotient_ref,
                        qap_0quotient_ref,
                        quotient_more_mem,
                        disable_bipartition_gp_local_search,
                        enable_mapping,
                        hierarchy_parameter_string,
                        distance_parameter_string,
                        online_distances,
                        filename_output,
                        map_construction_algorithm,
                        skip_map_ls,
                        delta_gains,
                        use_bin_id,
                        use_compact_bin_id,
                        full_matrix,
                        matching_type,
                        global_cycle_iterations,
                        suppress_output,
                        kway_fm_limits,
                        enable_convergence_map,
                        qap_label_iterations,
                        adapt_bal,
                        stream_buffer,
                        run_parallel,
                        use_fennel_objective,
                        ram_stream,
                        write_log,
                        evaluate,
                        stream_output_progress,
                        fennel_dynamics,
                        ghost_nodes_procedure,
                        stream_initial_bisections,
                        stream_allow_ghostnodes,
                        ghost_nodes_threshold,
                        num_streams_passes,
                        restream_vcycle,
                        batch_inbalance,
                        initial_part_multi_bfs,
                        initial_part_fennel,
                        skip_outer_ls,
                        stream_label_rounds,
                        automatic_buffer_len,
                        xxx,
                        max_buffer_size,
                        batch_size,
                        bq_disc_factor,
                        buffer_score,
                        haa_beta,
                        haa_theta,
                        ghost_neighbors_enabled,
                        sep_batch_marker,
                        restream_include_high_degree_nodes,
                        d_max,
                        bb_ratio,
                        use_queue,
                        dynamic_alpha,
                        batch_alpha,
                        minimal_mode,
                        num_split_edges,
                        past_subset_size,
                        reps,
#elif defined MODE_STREAMPARTITION
                        k,
                        imbalance,
                        filename_output,
                        stream_buffer,
                        run_parallel,
                        ram_stream,
                        write_log,
                        evaluate,
                        stream_output_progress,
                        stream_allow_ghostnodes,
                        num_streams_passes,
                        balance_edges,
                        benchmark,
                        light_evaluator,
                        label_propagation_iterations,
                        label_propagation_iterations_refinement,
                        max_buffer_size,
                        batch_size,
                        bq_disc_factor,
                        buffer_score,
                        haa_beta,
                        haa_theta,
                        ghost_neighbors_enabled,
                        sep_batch_marker,
                        restream_include_high_degree_nodes,
                        d_max,
                        bb_ratio,
                        use_queue,
                        dynamic_alpha,
                        batch_alpha,
                        minimal_mode,
                        num_split_edges,
                        past_subset_size,
                        reps,
                        xxx,
#elif defined MODE_SPMXV_MULTILEVELMAPPING
                        k,
                        imbalance,
                        preconfiguration,
                        time_limit,
                        enforce_balance,
                        balance_edges,
                        integrated_mapping,
                        multisection,
                        qap_label_propagation,
                        qap_blabel_propagation,
                        qap_alabel_propagation,
                        qap_multitry_fm,
                        qap_bmultitry_fm,
                        qap_kway_fm,
                        qap_bkway_fm,
                        qap_quotient_ref,
                        qap_bquotient_ref,
                        qap_0quotient_ref,
                        quotient_more_mem,
                        disable_bipartition_gp_local_search,
                        enable_mapping,
                        hierarchy_parameter_string,
                        distance_parameter_string,
                        online_distances,
                        filename_output,
                        map_construction_algorithm,
                        skip_map_ls,
                        delta_gains,
                        use_bin_id,
                        use_compact_bin_id,
                        full_matrix,
                        matching_type,
                        global_cycle_iterations,
                        suppress_output,
                        kway_fm_limits,
                        enable_convergence_map,
                        qap_label_iterations,
#elif defined MODE_MULTILEVELMAPPING
                        k,
                        imbalance,
                        preconfiguration,
                        time_limit,
                        enforce_balance,
                        balance_edges,
                        integrated_mapping,
                        multisection,
                        qap_label_propagation,
                        qap_blabel_propagation,
                        qap_alabel_propagation,
                        qap_multitry_fm,
                        qap_bmultitry_fm,
                        qap_kway_fm,
                        qap_bkway_fm,
                        qap_quotient_ref,
                        qap_bquotient_ref,
                        qap_0quotient_ref,
                        quotient_more_mem,
                        disable_bipartition_gp_local_search,
                        enable_mapping,
                        hierarchy_parameter_string,
                        distance_parameter_string,
                        online_distances,
                        filename_output,
                        map_construction_algorithm,
                        skip_map_ls,
                        delta_gains,
                        use_bin_id,
                        use_compact_bin_id,
                        full_matrix,
                        matching_type,
                        global_cycle_iterations,
                        suppress_output,
                        kway_fm_limits,
                        enable_convergence_map,
                        qap_label_iterations,
                        adapt_bal,
#elif defined MODE_KAFFPA
                        k,
                        imbalance,
                        preconfiguration,
                        time_limit,
                        enforce_balance,
                        balance_edges,
                        enable_mapping,
                        hierarchy_parameter_string,
                        distance_parameter_string,
                        online_distances,
                        filename_output,
#elif defined MODE_GRAPH_TRANSLATOR
                        filename_output,
#elif defined MODE_EVALUATOR
                        k,
                        preconfiguration,
                        input_partition,
                        k,
                        enable_mapping,
                        hierarchy_parameter_string,
                        distance_parameter_string,
                        online_distances,
                        use_bin_id,
                        use_compact_bin_id,
#elif defined MODE_NODESEP
                        // k,
                        imbalance, preconfiguration, filename_output,
    // time_limit,
    // edge_rating,
    // max_flow_improv_steps,
    // max_initial_ns_tries,
    // region_factor_node_separators,
    // global_cycle_iterations,
    // most_balanced_flows_node_sep,
    // sep_flows_disabled,
    // sep_fm_disabled,
    // sep_loc_fm_disabled,
    // sep_greedy_disabled,
    // sep_fm_unsucc_steps,
    // sep_num_fm_reps,
    // sep_loc_fm_unsucc_steps,
    // sep_num_loc_fm_reps,
    // sep_loc_fm_no_snodes,
    // sep_num_vert_stop,
    // sep_full_boundary_ip,
    // sep_edge_rating_during_ip,
    // sep_faster_ns,
#elif defined MODE_PARTITIONTOVERTEXSEPARATOR
                        k,    input_partition, filename_output,
#elif defined MODE_IMPROVEVERTEXSEPARATOR
                        input_partition,
                        filename_output,
#elif defined MODE_KAFFPAE
                        k,
                        imbalance,
                        preconfiguration,
                        time_limit,
                        mh_enable_quickstart,
                        mh_print_log,
                        mh_optimize_communication_volume,
                        mh_enable_tabu_search,
                        maxT,
                        maxIter,
                        mh_enable_kabapE,
                        kabaE_internal_bal,
                        balance_edges,
                        input_partition,
                        filename_output,
#elif defined MODE_LABELPROPAGATION
                        cluster_upperbound,
                        label_propagation_iterations,
                        filename_output,
#endif
                        end};
    // Parse arguments.
    int nerrors = arg_parse(argn, argv, argtable);

    // Catch case that help was requested.
    if (help->count > 0) {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout, argtable, "\n");
        arg_print_glossary(stdout, argtable, "  %-40s %s\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return 1;
    }

    if (nerrors > 0) {
        arg_print_errors(stderr, end, progname);
        printf("Try '%s --help' for more information.\n", progname);
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return 1;
    }

    if (k->count > 0) {
        partition_config.k = k->ival[0];
    }

    if (filename->count > 0) {
        graph_filename = filename->sval[0];
    }

    recursive = false;

    configuration cfg;
    cfg.standard(partition_config);
#ifdef MODE_MULTILEVELMAPPING
    cfg.eco(partition_config);
#else
    cfg.strong(partition_config);
#endif
#ifdef MODE_KAFFPA
    cfg.eco(partition_config);
#else
    cfg.strong(partition_config);
#endif

#ifdef MODE_NODESEP
    cfg.eco_separator(partition_config);
#endif

#ifdef MODE_STREAMPARTITION
    cfg.stream_partition(partition_config);
#endif

    if (preconfiguration->count > 0) {
#ifdef MODE_NODESEP
        if (strcmp("strong", preconfiguration->sval[0]) == 0) {
            cfg.strong_separator(partition_config);
        } else if (strcmp("eco", preconfiguration->sval[0]) == 0) {
            cfg.eco_separator(partition_config);
        } else if (strcmp("fast", preconfiguration->sval[0]) == 0) {
            cfg.fast_separator(partition_config);
        } else if (strcmp("fsocial", preconfiguration->sval[0]) == 0) {
            std::cout << "fsocial not supported yet" << '\n';
            exit(0);
        } else if (strcmp("esocial", preconfiguration->sval[0]) == 0) {
            std::cout << "esocial not supported yet" << '\n';
            exit(0);
        } else if (strcmp("ssocial", preconfiguration->sval[0]) == 0) {
            std::cout << "ssocial not supported yet" << '\n';
            exit(0);
        } else {
            fprintf(stderr, "Invalid preconfiguration variant: \"%s\"\n",
                    preconfiguration->sval[0]);
            exit(0);
        }
#else
        if (strcmp("strong", preconfiguration->sval[0]) == 0) {
            cfg.strong(partition_config);
        } else if (strcmp("eco", preconfiguration->sval[0]) == 0) {
            cfg.eco(partition_config);
        } else if (strcmp("fast", preconfiguration->sval[0]) == 0) {
            cfg.fast(partition_config);
        } else if (strcmp("fsocial", preconfiguration->sval[0]) == 0) {
            cfg.fastsocial(partition_config);
        } else if (strcmp("esocial", preconfiguration->sval[0]) == 0) {
            cfg.ecosocial(partition_config);
        } else if (strcmp("ssocial", preconfiguration->sval[0]) == 0) {
            cfg.strongsocial(partition_config);
        } else {
            fprintf(stderr, "Invalid preconfiguration variant: \"%s\"\n",
                    preconfiguration->sval[0]);
            exit(0);
        }
#endif
    }

    if (enable_mapping->count > 0) {
        partition_config.enable_mapping = true;
        if (!hierarchy_parameter_string->count) {
            std::cout
                << "Please specify the hierarchy using the --hierarchy_parameter_string option."
                << '\n';
            exit(0);
        }

        if (!distance_parameter_string->count) {
            std::cout
                << "Please specify the distances using the --distance_parameter_string option."
                << '\n';
            exit(0);
        }
    }

    if (hierarchy_parameter_string->count) {
        std::istringstream f(hierarchy_parameter_string->sval[0]);
        std::string s;
        partition_config.group_sizes.clear();
        while (getline(f, s, ':')) {
            partition_config.group_sizes.push_back(stoi(s));
        }

        PartitionID old_k = partition_config.k;
        partition_config.k = 1;  // recompute k
        for (unsigned int i = 0; i < partition_config.group_sizes.size(); i++) {
            partition_config.k *= partition_config.group_sizes[i];
        }
        if (old_k != partition_config.k) {
            std::cout
                << "number of processors defined through specified hierarchy does not match k!"
                << '\n';
            std::cout << "please specify k as " << partition_config.k << '\n';
            exit(0);
        }
    }

    if (distance_parameter_string->count) {
        std::istringstream f(distance_parameter_string->sval[0]);
        std::string s;
        partition_config.distances.clear();
        while (getline(f, s, ':')) {
            partition_config.distances.push_back(stoi(s));
        }
    }

    if (online_distances->count > 0) {
        partition_config.distance_construction_algorithm = DIST_CONST_HIERARCHY_ONLINE;
    }

    if (filename_output->count > 0) {
        partition_config.filename_output = filename_output->sval[0];
    }

    if (initial_partition_optimize->count > 0) {
        partition_config.initial_partition_optimize = true;
    }

    if (disable_balance_singletons->count > 0) {
        partition_config.use_balance_singletons = false;
    }

    if (mh_disable_nc_combine->count > 0) {
        partition_config.mh_disable_nc_combine = true;
    }

    if (mh_disable_cross_combine->count > 0) {
        partition_config.mh_disable_cross_combine = true;
    }

    if (imbalance->count > 0) {
        partition_config.epsilon = imbalance->dval[0];
    }

    if (mh_disable_combine->count > 0) {
        partition_config.mh_disable_combine = true;
    }

    if (balance_edges->count > 0) {
        partition_config.balance_edges = true;
    }

    if (mh_optimize_communication_volume->count > 0) {
        partition_config.mh_optimize_communication_volume = true;
    }

    if (mh_enable_tournament_selection->count > 0) {
        partition_config.mh_enable_tournament_selection = true;
    }

    if (amg_iterations->count > 0) {
        partition_config.amg_iterations = amg_iterations->ival[0];
    }

    if (max_initial_ns_tries->count > 0) {
        partition_config.max_initial_ns_tries = max_initial_ns_tries->ival[0];
    }

    if (max_flow_improv_steps->count > 0) {
        partition_config.max_flow_improv_steps = max_flow_improv_steps->ival[0];
    }

    if (region_factor_node_separators->count > 0) {
        partition_config.region_factor_node_separators = region_factor_node_separators->dval[0];
    }

    if (kabaE_internal_bal->count > 0) {
        partition_config.kabaE_internal_bal = kabaE_internal_bal->dval[0];
    }

    if (kaba_internal_no_aug_steps_aug->count > 0) {
        partition_config.kaba_internal_no_aug_steps_aug = kaba_internal_no_aug_steps_aug->ival[0];
    }

    if (kaffpa_perfectly_balanced_refinement->count > 0) {
        partition_config.kaffpa_perfectly_balanced_refinement = true;
    }

    if (kaba_disable_zero_weight_cycles->count > 0) {
        partition_config.kaba_enable_zero_weight_cycles = false;
    }

    if (kaba_unsucc_iterations->count > 0) {
        partition_config.kaba_unsucc_iterations = kaba_unsucc_iterations->ival[0];
    }

    if (kaba_flip_packings->count > 0) {
        partition_config.kaba_flip_packings = true;
    }

    if (sep_flows_disabled->count > 0) {
        partition_config.sep_flows_disabled = true;
    }

    if (sep_faster_ns->count > 0) {
        partition_config.faster_ns = true;
    }

    if (sep_loc_fm_no_snodes->count > 0) {
        partition_config.sep_loc_fm_no_snodes = sep_loc_fm_no_snodes->ival[0];
    }

    if (sep_num_vert_stop->count > 0) {
        partition_config.sep_num_vert_stop = sep_num_vert_stop->ival[0];
    }

    if (sep_fm_unsucc_steps->count > 0) {
        partition_config.sep_fm_unsucc_steps = sep_fm_unsucc_steps->ival[0];
    }
    if (sep_fm_unsucc_steps->count > 0) {
        partition_config.sep_fm_unsucc_steps = sep_fm_unsucc_steps->ival[0];
    }

    if (sep_num_fm_reps->count > 0) {
        partition_config.sep_num_fm_reps = sep_num_fm_reps->ival[0];
    }

    if (sep_loc_fm_unsucc_steps->count > 0) {
        partition_config.sep_loc_fm_unsucc_steps = sep_loc_fm_unsucc_steps->ival[0];
    }

    if (sep_num_loc_fm_reps->count > 0) {
        partition_config.sep_num_loc_fm_reps = sep_num_loc_fm_reps->ival[0];
    }

    if (sep_fm_disabled->count > 0) {
        partition_config.sep_fm_disabled = true;
    }

    if (sep_loc_fm_disabled->count > 0) {
        partition_config.sep_loc_fm_disabled = true;
    }

    if (sep_greedy_disabled->count > 0) {
        partition_config.sep_greedy_disabled = true;
    }

    if (sep_full_boundary_ip->count > 0) {
        partition_config.sep_full_boundary_ip = true;
    }

    if (kaba_lsearch_p->count) {
        if (strcmp("coindiff", kaba_lsearch_p->sval[0]) == 0) {
            partition_config.kaba_lsearch_p = COIN_DIFFTIE;
        } else if (strcmp("nocoindiff", kaba_lsearch_p->sval[0]) == 0) {
            partition_config.kaba_lsearch_p = NOCOIN_DIFFTIE;
        } else if (strcmp("coinrnd", kaba_lsearch_p->sval[0]) == 0) {
            partition_config.kaba_lsearch_p = COIN_RNDTIE;
        } else if (strcmp("nocoinrnd", kaba_lsearch_p->sval[0]) == 0) {
            partition_config.kaba_lsearch_p = NOCOIN_RNDTIE;
        } else {
            fprintf(stderr, "Invalid combine variant: \"%s\"\n", kaba_lsearch_p->sval[0]);
            exit(0);
        }
    }

    if (maxT->count > 0) {
        partition_config.maxT = maxT->ival[0];
    }

    if (maxIter->count > 0) {
        partition_config.maxIter = maxIter->ival[0];
    }

    if (mh_enable_tabu_search->count > 0) {
        partition_config.mh_enable_gal_combine = true;
    }

    if (kaba_packing_iterations->count > 0) {
        partition_config.kaba_packing_iterations = kaba_packing_iterations->ival[0];
    }

    if (mh_flip_coin->count > 0) {
        partition_config.mh_flip_coin = mh_flip_coin->ival[0];
    }

    if (mh_initial_population_fraction->count > 0) {
        partition_config.mh_initial_population_fraction = mh_initial_population_fraction->ival[0];
    }

    if (minipreps->count > 0) {
        partition_config.minipreps = minipreps->ival[0];
    }

    if (mh_enable_quickstart->count > 0) {
        partition_config.mh_enable_quickstart = true;
    }

    if (mh_disable_diversify_islands->count > 0) {
        partition_config.mh_disable_diversify_islands = true;
    }

    if (gpa_grow_internal->count > 0) {
        partition_config.gpa_grow_paths_between_blocks = false;
    }

    if (suppress_output->count > 0) {
        suppress_program_output = true;
        partition_config.suppress_output = true;
    }

    if (mh_print_log->count > 0) {
        partition_config.mh_print_log = true;
    }

    if (use_bucket_queues->count > 0) {
        partition_config.use_bucket_queues = true;
    }

    if (recursive_bipartitioning->count > 0) {
        recursive = true;
    }

    if (time_limit->count > 0) {
        partition_config.time_limit = time_limit->dval[0];
    }

    if (unsuccessful_reps->count > 0) {
        partition_config.no_unsuc_reps = unsuccessful_reps->ival[0];
    }

    if (mh_pool_size->count > 0) {
        partition_config.mh_pool_size = mh_pool_size->ival[0];
    }

    if (mh_penalty_for_unconnected->count > 0) {
        partition_config.mh_penalty_for_unconnected = true;
    }

    if (mh_enable_kabapE->count > 0) {
        partition_config.kabapE = true;
    }

    if (initial_partition_optimize_multitry_fm_alpha->count > 0) {
        partition_config.initial_partition_optimize_multitry_fm_alpha =
            initial_partition_optimize_multitry_fm_alpha->ival[0];
    }

    if (initial_partition_optimize_multitry_rounds->count > 0) {
        partition_config.initial_partition_optimize_multitry_rounds =
            initial_partition_optimize_multitry_rounds->ival[0];
    }

    if (initial_partition_optimize_fm_limits->count > 0) {
        partition_config.initial_partition_optimize_fm_limits =
            initial_partition_optimize_fm_limits->ival[0];
    }

    if (mh_disable_diversify->count > 0) {
        partition_config.mh_diversify = false;
    }

    if (mh_diversify_best->count > 0) {
        partition_config.mh_diversify_best = true;
    }

    if (enforce_balance->count > 0) {
        partition_config.kaffpa_perfectly_balance = true;
    }

    if (mh_plain_repetitions->count > 0) {
        partition_config.mh_plain_repetitions = true;
    }

    if (local_partitioning_repetitions->count > 0) {
        partition_config.local_partitioning_repetitions = local_partitioning_repetitions->ival[0];
    }

    if (only_first_level->count > 0) {
        partition_config.only_first_level = true;
    }

    if (mh_cross_combine_original_k->count > 0) {
        partition_config.mh_cross_combine_original_k = true;
    }

    if (mh_sequential_mode->count > 0) {
        partition_config.mh_no_mh = true;
    }

    if (enable_omp->count > 0) {
        partition_config.enable_omp = true;
    }

    if (compute_vertex_separator->count > 0) {
        partition_config.compute_vertex_separator = true;
    }

    if (most_balanced_flows->count > 0) {
        partition_config.most_balanced_minimum_cuts = true;
    }

    if (most_balanced_flows_node_sep->count > 0) {
        partition_config.most_balanced_minimum_cuts_node_sep = true;
    }

    if (use_wcycles->count > 0) {
        partition_config.use_wcycles = true;
    }

    if (enable_convergence->count > 0) {
        partition_config.no_change_convergence = true;
    }

    if (use_fullmultigrid->count > 0) {
        partition_config.use_fullmultigrid = true;
    }

    if (use_vcycle->count > 0) {
        partition_config.use_fullmultigrid = false;
        partition_config.use_wcycles = false;
    }

    if (toposort_iterations->count > 0) {
        partition_config.toposort_iterations = toposort_iterations->ival[0];
    }

    if (bipartition_tries->count > 0) {
        partition_config.bipartition_tries = bipartition_tries->ival[0];
    }

    if (bipartition_post_fm_limit->count > 0) {
        partition_config.bipartition_post_fm_limits = bipartition_post_fm_limit->ival[0];
    }

    if (bipartition_post_ml_limit->count > 0) {
        partition_config.bipartition_post_ml_limits = bipartition_post_ml_limit->ival[0];
    }

    if (disable_max_vertex_weight_constraint->count > 0) {
        partition_config.disable_max_vertex_weight_constraint = true;
    }

    if (num_vert_stop_factor->count > 0) {
        partition_config.num_vert_stop_factor = num_vert_stop_factor->ival[0];
    }

    if (local_multitry_rounds->count > 0) {
        partition_config.local_multitry_rounds = local_multitry_rounds->ival[0];
    }

    if (local_multitry_fm_alpha->count > 0) {
        partition_config.local_multitry_fm_alpha = local_multitry_fm_alpha->ival[0];
    }

    if (wcycle_no_new_initial_partitioning->count > 0) {
        partition_config.no_new_initial_partitioning = true;
    }

    if (graph_weighted->count > 0) {
        is_graph_weighted = true;
    }

    if (disable_refined_bubbling->count > 0) {
        partition_config.refined_bubbling = false;
    }

    if (input_partition->count > 0) {
        partition_config.input_partition = input_partition->sval[0];
    }

    if (global_cycle_iterations->count > 0) {
        partition_config.global_cycle_iterations = global_cycle_iterations->ival[0];
    }

    if (level_split->count > 0) {
        partition_config.level_split = level_split->ival[0];
    }

    if (disable_qgraph_refinement->count > 0) {
        partition_config.quotient_graph_refinement_disabled = true;
    }

    if (bubbling_iterations->count > 0) {
        partition_config.bubbling_iterations = bubbling_iterations->ival[0];
    }

    if (kway_fm_limits->count > 0) {
        partition_config.kway_fm_search_limit = kway_fm_limits->ival[0];
    }

    if (kway_rounds->count > 0) {
        partition_config.kway_rounds = kway_rounds->ival[0];
    }

    if (enable_corner_refinement->count > 0) {
        partition_config.corner_refinement_enabled = true;
    }

    if (match_islands->count > 0) {
        partition_config.match_islands = true;
    }

    if (aggressive_random_levels->count > 0) {
        partition_config.aggressive_random_levels = aggressive_random_levels->ival[0];
    }

    if (rate_first_level_inner_outer->count > 0) {
        partition_config.rate_first_level_inner_outer = true;
    }

    if (user_seed->count > 0) {
        partition_config.seed = user_seed->ival[0];
    }

    if (fm_search_limit->count > 0) {
        partition_config.fm_search_limit = fm_search_limit->ival[0];
    }

    if (bank_account_factor->count > 0) {
        partition_config.bank_account_factor = bank_account_factor->dval[0];
    }

    if (flow_region_factor->count > 0) {
        partition_config.flow_region_factor = flow_region_factor->dval[0];
    }

    if (kway_adaptive_limits_alpha->count > 0) {
        partition_config.kway_adaptive_limits_alpha = kway_adaptive_limits_alpha->dval[0];
    }

    if (imbalance->count > 0) {
        partition_config.imbalance = imbalance->dval[0];
        partition_config.stream_global_epsilon = (partition_config.imbalance) / 100.;
    }

    if (initial_partitioning_repetitions->count > 0) {
        partition_config.initial_partitioning_repetitions =
            initial_partitioning_repetitions->ival[0];
    }

    if (edge_rating_tiebreaking->count > 0) {
        partition_config.edge_rating_tiebreaking = true;
    }

    if (first_level_random_matching->count > 0) {
        partition_config.first_level_random_matching = true;
    } else {
        partition_config.first_level_random_matching = false;
    }

    if (kaba_neg_cycle_algorithm->count > 0) {
        if (strcmp("playfield", kaba_neg_cycle_algorithm->sval[0]) == 0) {
            partition_config.cycle_refinement_algorithm = CYCLE_REFINEMENT_ALGORITHM_PLAYFIELD;
        } else if (strcmp("ultramodel", kaba_neg_cycle_algorithm->sval[0]) == 0) {
            partition_config.cycle_refinement_algorithm = CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL;
        } else if (strcmp("ultramodelplus", kaba_neg_cycle_algorithm->sval[0]) == 0) {
            partition_config.cycle_refinement_algorithm =
                CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL_PLUS;
        } else {
            fprintf(stderr, "Invalid balanced refinement operator: \"%s\"\n",
                    kaba_neg_cycle_algorithm->sval[0]);
            exit(0);
        }
    }

    if (sep_edge_rating_during_ip->count > 0) {
        if (strcmp("expansionstar", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = EXPANSIONSTAR;
        } else if (strcmp("expansionstar2", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = EXPANSIONSTAR2;
        } else if (strcmp("expansionstar2algdist", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = EXPANSIONSTAR2ALGDIST;
        } else if (strcmp("geom", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = PSEUDOGEOM;
        } else if (strcmp("sepaddx", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_ADDX;
        } else if (strcmp("sepmultx", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_MULTX;
        } else if (strcmp("sepmax", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_MAX;
        } else if (strcmp("seplog", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_LOG;
        } else if (strcmp("r1", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_R1;
        } else if (strcmp("r2", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_R2;
        } else if (strcmp("r3", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_R3;
        } else if (strcmp("r4", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_R4;
        } else if (strcmp("r5", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_R5;
        } else if (strcmp("r6", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_R6;
        } else if (strcmp("r7", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_R7;
        } else if (strcmp("r8", sep_edge_rating_during_ip->sval[0]) == 0) {
            partition_config.sep_edge_rating_during_ip = SEPARATOR_R8;
        } else {
            fprintf(stderr, "Invalid edge rating variant: \"%s\"\n",
                    sep_edge_rating_during_ip->sval[0]);
            exit(0);
        }
    }

    if (edge_rating->count > 0) {
        if (strcmp("expansionstar", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = EXPANSIONSTAR;
        } else if (strcmp("expansionstar2", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = EXPANSIONSTAR2;
        } else if (strcmp("weight", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = WEIGHT;
        } else if (strcmp("realweight", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = REALWEIGHT;
        } else if (strcmp("expansionstar2algdist", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = EXPANSIONSTAR2ALGDIST;
        } else if (strcmp("geom", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = PSEUDOGEOM;
        } else if (strcmp("sepaddx", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_ADDX;
        } else if (strcmp("sepmultx", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_MULTX;
        } else if (strcmp("sepmax", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_MAX;
        } else if (strcmp("seplog", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_LOG;
        } else if (strcmp("r1", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_R1;
        } else if (strcmp("r2", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_R2;
        } else if (strcmp("r3", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_R3;
        } else if (strcmp("r4", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_R4;
        } else if (strcmp("r5", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_R5;
        } else if (strcmp("r6", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_R6;
        } else if (strcmp("r7", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_R7;
        } else if (strcmp("r8", edge_rating->sval[0]) == 0) {
            partition_config.edge_rating = SEPARATOR_R8;
        } else {
            fprintf(stderr, "Invalid edge rating variant: \"%s\"\n", edge_rating->sval[0]);
            exit(0);
        }
    }

    if (bipartition_algorithm->count > 0) {
        if (strcmp("bfs", bipartition_algorithm->sval[0]) == 0) {
            partition_config.bipartition_algorithm = BIPARTITION_BFS;
        } else if (strcmp("fm", bipartition_algorithm->sval[0]) == 0) {
            partition_config.bipartition_algorithm = BIPARTITION_FM;
        } else {
            fprintf(stderr, "Invalid bipartition algorthim: \"%s\"\n",
                    bipartition_algorithm->sval[0]);
            exit(0);
        }
    }

    if (refinement_scheduling_algorithm->count > 0) {
        if (strcmp("fast", refinement_scheduling_algorithm->sval[0]) == 0) {
            partition_config.refinement_scheduling_algorithm = REFINEMENT_SCHEDULING_FAST;
        } else if (strcmp("active_blocks", refinement_scheduling_algorithm->sval[0]) == 0) {
            partition_config.refinement_scheduling_algorithm = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS;
        } else if (strcmp("active_blocks_kway", refinement_scheduling_algorithm->sval[0]) == 0) {
            partition_config.refinement_scheduling_algorithm =
                REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY;
        } else {
            fprintf(stderr, "Invalid refinement scheduling variant: \"%s\"\n",
                    refinement_scheduling_algorithm->sval[0]);
            exit(0);
        }
    }

    if (stop_rule->count > 0) {
        if (strcmp("simple", stop_rule->sval[0]) == 0) {
            partition_config.stop_rule = STOP_RULE_SIMPLE;
        } else if (strcmp("multiplek", stop_rule->sval[0]) == 0) {
            partition_config.stop_rule = STOP_RULE_MULTIPLE_K;
        } else if (strcmp("strong", stop_rule->sval[0]) == 0) {
            partition_config.stop_rule = STOP_RULE_STRONG;
        } else {
            fprintf(stderr, "Invalid stop rule: \"%s\"\n", stop_rule->sval[0]);
            exit(0);
        }
    }

    if (kway_search_stop_rule->count > 0) {
        if (strcmp("simple", kway_search_stop_rule->sval[0]) == 0) {
            partition_config.kway_stop_rule = KWAY_SIMPLE_STOP_RULE;
        } else if (strcmp("adaptive", kway_search_stop_rule->sval[0]) == 0) {
            partition_config.kway_stop_rule = KWAY_ADAPTIVE_STOP_RULE;
        } else {
            fprintf(stderr, "Invalid kway stop rule: \"%s\"\n", kway_search_stop_rule->sval[0]);
            exit(0);
        }
    }

    if (permutation_quality->count > 0) {
        if (strcmp("none", permutation_quality->sval[0]) == 0) {
            partition_config.permutation_quality = PERMUTATION_QUALITY_NONE;
        } else if (strcmp("fast", permutation_quality->sval[0]) == 0) {
            partition_config.permutation_quality = PERMUTATION_QUALITY_FAST;
        } else if (strcmp("good", permutation_quality->sval[0]) == 0) {
            partition_config.permutation_quality = PERMUTATION_QUALITY_GOOD;
        } else {
            fprintf(stderr, "Invalid permutation quality variant: \"%s\"\n",
                    permutation_quality->sval[0]);
            exit(0);
        }
    }

    if (permutation_during_refinement->count > 0) {
        if (strcmp("none", permutation_during_refinement->sval[0]) == 0) {
            partition_config.permutation_during_refinement = PERMUTATION_QUALITY_NONE;
        } else if (strcmp("fast", permutation_during_refinement->sval[0]) == 0) {
            partition_config.permutation_during_refinement = PERMUTATION_QUALITY_FAST;
        } else if (strcmp("good", permutation_during_refinement->sval[0]) == 0) {
            partition_config.permutation_during_refinement = PERMUTATION_QUALITY_GOOD;
        } else {
            fprintf(stderr, "Invalid permutation quality during refinement variant: \"%s\"\n",
                    permutation_during_refinement->sval[0]);
            exit(0);
        }
    }

    if (matching_type->count > 0) {
        if (strcmp("random", matching_type->sval[0]) == 0) {
            partition_config.matching_type = MATCHING_RANDOM;
        } else if (strcmp("gpa", matching_type->sval[0]) == 0) {
            partition_config.matching_type = MATCHING_GPA;
        } else if (strcmp("randomgpa", matching_type->sval[0]) == 0) {
            partition_config.matching_type = MATCHING_RANDOM_GPA;
        } else if (strcmp("cluster", matching_type->sval[0]) == 0) {
            partition_config.matching_type = CLUSTER_COARSENING;
        } else {
            fprintf(stderr, "Invalid matching variant: \"%s\"\n", matching_type->sval[0]);
            exit(0);
        }
    }

    if (refinement_type->count > 0) {
        if (strcmp("fm", refinement_type->sval[0]) == 0) {
            partition_config.refinement_type = REFINEMENT_TYPE_FM;
        } else if (strcmp("fm_flow", refinement_type->sval[0]) == 0) {
            partition_config.refinement_type = REFINEMENT_TYPE_FM_FLOW;
        } else if (strcmp("flow", refinement_type->sval[0]) == 0) {
            partition_config.refinement_type = REFINEMENT_TYPE_FLOW;
        } else {
            fprintf(stderr, "Invalid refinement type variant: \"%s\"\n", refinement_type->sval[0]);
            exit(0);
        }
    }

    if (initial_partition->count > 0) {
        if (strcmp("recursive", initial_partition->sval[0]) == 0) {
            partition_config.initial_partitioning_type = INITIAL_PARTITIONING_RECPARTITION;
        } else {
            fprintf(stderr, "Invalid initial partition variant: \"%s\"\n",
                    initial_partition->sval[0]);
            exit(0);
        }
    }

    if (label_propagation_iterations->count > 0) {
        partition_config.label_iterations = label_propagation_iterations->ival[0];
    }

    if (cluster_upperbound->count > 0) {
        partition_config.cluster_upperbound = cluster_upperbound->ival[0];
    } else {
        partition_config.cluster_upperbound = std::numeric_limits<NodeWeight>::max() / 2;
    }

    if (integrated_mapping->count > 0) {
        partition_config.integrated_mapping = true;
        cfg.integrated_mapping(partition_config);
        if (!hierarchy_parameter_string->count) {
            std::cout
                << "Please specify the hierarchy using the --hierarchy_parameter_string option."
                << '\n';
            exit(0);
        }

        if (!distance_parameter_string->count) {
            std::cout
                << "Please specify the distances using the --distance_parameter_string option."
                << '\n';
            exit(0);
        }
    }

    if (multisection->count > 0) {
        partition_config.multisection = true;
        if (!integrated_mapping->count && !enable_mapping->count) {
            std::cout << "--multisection only works together with either of these flags: "
                      << "--integrated_mapping , --enable_mapping" << '\n';
            exit(0);
        }
    }

    if (qap_label_propagation->count > 0) {
        partition_config.qap_label_propagation_refinement = true;
        if (!integrated_mapping->count) {
            std::cout
                << "--qap_label_propagation only works together with the --integrated_mapping flag."
                << '\n';
            exit(0);
        }
    }

    if (qap_blabel_propagation->count > 0) {
        partition_config.qap_blabel_propagation_refinement = true;
        if (!integrated_mapping->count) {
            std::cout << "--qap_blabel_propagation only works together with the "
                         "--integrated_mapping flag."
                      << '\n';
            exit(0);
        }
    }

    if (qap_alabel_propagation->count > 0) {
        partition_config.qap_alabel_propagation_refinement = true;
        if (!integrated_mapping->count) {
            std::cout << "--qap_alabel_propagation only works together with the "
                         "--integrated_mapping flag."
                      << '\n';
            exit(0);
        }
    }

    if (qap_multitry_fm->count > 0) {
        partition_config.qap_multitry_kway_fm = true;
        if (!integrated_mapping->count) {
            std::cout << "--qap_multitry_fm only works together with the --integrated_mapping flag."
                      << '\n';
            exit(0);
        }
    }

    if (qap_bmultitry_fm->count > 0) {
        partition_config.qap_bmultitry_kway_fm = true;
        if (!integrated_mapping->count) {
            std::cout
                << "--qap_bmultitry_fm only works together with the --integrated_mapping flag."
                << '\n';
            exit(0);
        }
    }

    if (qap_kway_fm->count > 0) {
        partition_config.qap_kway_fm = true;
        if (!integrated_mapping->count) {
            std::cout << "--qap_kway_fm only works together with the --integrated_mapping flag."
                      << '\n';
            exit(0);
        }
    }

    if (qap_bkway_fm->count > 0) {
        partition_config.qap_bkway_fm = true;
        if (!integrated_mapping->count) {
            std::cout << "--qap_bkway_fm only works together with the --integrated_mapping flag."
                      << '\n';
            exit(0);
        }
    }

    if (enable_convergence_map->count > 0) {
        partition_config.no_change_convergence_map = true;
        if (!integrated_mapping->count) {
            std::cout << "--enable_convergence_map only works together with the "
                         "--integrated_mapping flag."
                      << '\n';
            exit(0);
        }
    }

    if (qap_quotient_ref->count > 0) {
        partition_config.qap_quotient_ref = true;
    }

    if (qap_bquotient_ref->count > 0) {
        partition_config.qap_bquotient_ref = true;
    }

    if (qap_0quotient_ref->count > 0) {
        partition_config.qap_bquotient_ref = true;
        partition_config.qap_0quotient_ref = true;
    }

    if (quotient_more_mem->count > 0) {
        partition_config.quotient_more_mem = true;
        if (!integrated_mapping->count && !use_fennel_objective->count) {
            std::cout << "--qap_quotient_ref only works together with --integrated_mapping or "
                         "--use_fennel_objective."
                      << '\n';
            exit(0);
        }
    }

    if (disable_bipartition_gp_local_search->count > 0) {
        partition_config.bipartition_gp_local_search = false;
        if (!integrated_mapping->count) {
            std::cout << "--disable_bipartition_gp_local_search only works together with the "
                         "--integrated_mapping flag."
                      << '\n';
            exit(0);
        }
    }

    if (map_construction_algorithm->count) {
        if (strcmp("oldgrowingfaster", map_construction_algorithm->sval[0]) == 0) {
            partition_config.construction_algorithm = MAP_CONST_OLDGROWING_FASTER;
        } else if (strcmp("random", map_construction_algorithm->sval[0]) == 0) {
            partition_config.construction_algorithm = MAP_CONST_RANDOM;
        } else if (strcmp("identity", map_construction_algorithm->sval[0]) == 0) {
            partition_config.construction_algorithm = MAP_CONST_IDENTITY;
        } else if (strcmp("fasthierarchybottomup", map_construction_algorithm->sval[0]) == 0) {
            partition_config.construction_algorithm = MAP_CONST_FASTHIERARCHY_BOTTOMUP;
        } else if (strcmp("fasthierarchytopdown", map_construction_algorithm->sval[0]) == 0) {
            partition_config.construction_algorithm = MAP_CONST_FASTHIERARCHY_TOPDOWN;
        } else if (strcmp("greedyallc", map_construction_algorithm->sval[0]) == 0) {
            partition_config.construction_algorithm = MAP_CONST_OLDGROWING;
        } else {
            std::cout << "map_construction_algorithm should be one of the following options:"
                      << '\n';
            std::cout << "oldgrowingfaster random identity fasthierarchybottomup "
                         "fasthierarchytopdown greedyallc"
                      << '\n';
            exit(0);
        }
    }

    if (skip_map_ls->count > 0) {
        partition_config.skip_map_ls = true;
        if (!enable_mapping->count && !integrated_mapping->count) {
            std::cout << "--skip_map_ls only works together with either --enable_mapping or "
                         "--integrated_mapping"
                      << '\n';
            exit(0);
        }
    }

    if (delta_gains->count > 0) {
        partition_config.use_delta_gains = true;
        if (!integrated_mapping->count) {
            std::cout << "--delta_gains only works together with --integrated_mapping" << '\n';
            exit(0);
        }
    }

    if (use_bin_id->count > 0) {
        partition_config.use_bin_id = true;
        if (!integrated_mapping->count) {
            std::cout << "--use_bin_id only works together with --integrated_mapping" << '\n';
            exit(0);
        }
        if (!online_distances->count) {
            std::cout << "--use_bin_id only works together with --online_distances" << '\n';
            exit(0);
        }
    }

    if (use_compact_bin_id->count > 0) {
        partition_config.use_compact_bin_id = true;
        if (!integrated_mapping->count) {
            std::cout << "--use_compact_bin_id only works together with --integrated_mapping"
                      << '\n';
            exit(0);
        }
        if (!online_distances->count) {
            std::cout << "--use_compact_bin_id only works together with --online_distances" << '\n';
            exit(0);
        }
        if (use_bin_id->count > 0) {
            std::cout
                << "--use_bin_id only and --use_compact_bin_id cannot be used at the same time"
                << '\n';
            exit(0);
        }
    }

    if (full_matrix->count > 0) {
        partition_config.full_matrix = true;
        if (!integrated_mapping->count) {
            std::cout << "--full_matrix only works together with --integrated_mapping" << '\n';
            exit(0);
        }
        if (online_distances->count > 0) {
            std::cout << "--full_matrix and --online_distances cannot be used at the same time"
                      << '\n';
            exit(0);
        }
    }

    if (qap_label_iterations->count > 0) {
        partition_config.label_iterations_refinement_map = qap_label_iterations->ival[0];
    }

    if (adapt_bal->count > 0) {
        partition_config.adapt_bal = true;
    }

    cli::stream_parameters::NodeStreamArgRefs node_stream_args{stream_buffer,
                                                               run_parallel,
                                                               use_queue,
                                                               batch_size,
                                                               max_buffer_size,
                                                               bq_disc_factor,
                                                               buffer_score,
                                                               haa_beta,
                                                               haa_theta,
                                                               ghost_neighbors_enabled,
                                                               sep_batch_marker,
                                                               bb_ratio,
                                                               restream_include_high_degree_nodes,
                                                               d_max,
                                                               use_fennel_objective,
                                                               ram_stream,
                                                               write_log,
                                                               evaluate,
                                                               stream_output_progress,
                                                               stream_allow_ghostnodes,
                                                               stream_initial_bisections,
                                                               fennel_dynamics,
                                                               ghost_nodes_procedure,
                                                               ghost_nodes_threshold,
                                                               num_streams_passes,
                                                               batch_inbalance,
                                                               restream_vcycle,
                                                               initial_part_multi_bfs,
                                                               skip_outer_ls,
                                                               initial_part_fennel,
                                                               stream_label_rounds,
                                                               automatic_buffer_len,
                                                               xxx};
    cli::stream_parameters::apply_node_stream_options(node_stream_args, partition_config);

    // stream edge partition
    cli::stream_parameters::EdgeStreamArgRefs edge_stream_args{
        use_queue,
        dynamic_alpha,
        batch_alpha,
        minimal_mode,
        num_split_edges,
        past_subset_size,
        reps,
        benchmark,
        light_evaluator,
        label_propagation_iterations_refinement};
    cli::stream_parameters::apply_edge_stream_options(edge_stream_args, partition_config);

    return 0;
}

// Mode-specific wrappers to keep entry points explicit.
inline int parse_parameters_node(int argn, char** argv, Config& partition_config,
                                 std::string& graph_filename, bool& is_graph_weighted,
                                 bool& suppress_program_output, bool& recursive) {
    return parse_parameters(argn, argv, partition_config, graph_filename, is_graph_weighted,
                            suppress_program_output, recursive);
}

inline int parse_parameters_edge(int argn, char** argv, Config& partition_config,
                                 std::string& graph_filename, bool& is_graph_weighted,
                                 bool& suppress_program_output, bool& recursive) {
    partition_config.edge_partition = true;
    const int ret = parse_parameters(argn, argv, partition_config, graph_filename,
                                     is_graph_weighted, suppress_program_output, recursive);
    if (ret != 0) {
        return ret;
    }

    bool batch_size_overridden = false;
    bool run_parallel_requested = false;
    bool buffer_size_requested = false;
    for (int i = 1; i < argn; ++i) {
        const char* arg = argv[i];
        if (std::strncmp(arg, "--batch_size", 12) == 0 ||
            std::strncmp(arg, "--stream_buffer", 15) == 0) {
            batch_size_overridden = true;
        }
        if (std::strncmp(arg, "--run-parallel", 14) == 0) {
            run_parallel_requested = true;
        }
        if (std::strncmp(arg, "--buffer_size", 13) == 0) {
            buffer_size_requested = true;
        }
    }

    if (run_parallel_requested) {
        std::cout << "[heistream_edge] --run-parallel is not implemented for edge partitioning. "
                     "Proceeding with the sequential edge pipeline."
                  << '\n';
    }

    if (buffer_size_requested) {
        std::cout << "[heistream_edge] --buffer_size is not implemented for edge partitioning. "
                     "The option is ignored."
                  << '\n';
    }

    if (!batch_size_overridden) {
        // Legacy edge default: stream buffer (effective batch) is 32768.
        partition_config.batch_size = 32768;
    }
    return 0;
}

#endif /* end of include guard: PARSE_PARAMETERS_GPJMGSM8 */
