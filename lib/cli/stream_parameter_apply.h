/****************************************************************************
 * stream_parameter_apply.h
 *****************************************************************************/

#ifndef STREAM_PARAMETER_APPLY_H_
#define STREAM_PARAMETER_APPLY_H_

#include <argtable3.h>

#include <cstring>
#include <iostream>
#include <limits>

#include "partition/partition_config.h"


namespace cli::stream_parameters {

struct NodeStreamArgRefs {
    arg_int* stream_buffer;
    arg_lit* run_parallel;
    arg_str* use_queue;
    arg_int* batch_size;
    arg_int* max_buffer_size;
    arg_int* bq_disc_factor;
    arg_str* buffer_score;
    arg_dbl* haa_beta;
    arg_dbl* haa_theta;
    arg_lit* ghost_neighbors_enabled;
    arg_lit* sep_batch_marker;
    arg_int* bb_ratio;
    arg_lit* restream_include_high_degree_nodes;
    arg_int* d_max;
    arg_lit* use_fennel_objective;
    arg_lit* ram_stream;
    arg_lit* write_log;
    arg_str* evaluate;
    arg_lit* stream_output_progress;
    arg_lit* stream_allow_ghostnodes;
    arg_lit* stream_initial_bisections;
    arg_str* fennel_dynamics;
    arg_str* ghost_nodes_procedure;
    arg_dbl* ghost_nodes_threshold;
    arg_dbl* num_streams_passes;
    arg_dbl* batch_inbalance;
    arg_lit* restream_vcycle;
    arg_lit* initial_part_multi_bfs;
    arg_lit* skip_outer_ls;
    arg_lit* initial_part_fennel;
    arg_str* stream_label_rounds;
    arg_lit* automatic_buffer_len;
    arg_int* xxx;
};

struct EdgeStreamArgRefs {
    arg_str* use_queue;
    arg_lit* dynamic_alpha;
    arg_lit* batch_alpha;
    arg_lit* minimal_mode;
    arg_int* num_split_edges;
    arg_int* past_subset_size;
    arg_dbl* reps;
    arg_lit* benchmark;
    arg_lit* light_evaluator;
    arg_int* label_propagation_iterations_refinement;
};

inline void apply_node_stream_options(const NodeStreamArgRefs& args, Config& partition_config) {
    auto parse_bool_option = [](const char* value, const char* name) {
        if (strcmp("true", value) == 0 || strcmp("1", value) == 0) {
            return true;
        }
        if (strcmp("false", value) == 0 || strcmp("0", value) == 0) {
            return false;
        }
        std::cout << name << " should be one of: true, false, 1, 0" << '\n';
        exit(0);
    };

    if (args.stream_buffer->count > 0) {
        partition_config.batch_size = static_cast<LongNodeID>(args.stream_buffer->ival[0]);
    }

    if (args.batch_size->count > 0) {
        partition_config.batch_size = static_cast<LongNodeID>(args.batch_size->ival[0]);
    }

    if (args.max_buffer_size->count > 0) {
        partition_config.max_buffer_size = args.max_buffer_size->ival[0];
    }

    if (args.bq_disc_factor->count > 0) {
        partition_config.bq_disc_factor = args.bq_disc_factor->ival[0];
    }

    if (args.haa_beta->count > 0) {
        partition_config.haa_beta = args.haa_beta->dval[0];
    }

    if (args.haa_theta->count > 0) {
        partition_config.haa_theta = args.haa_theta->dval[0];
    }

    if (args.ghost_neighbors_enabled->count > 0) {
        partition_config.ghost_neighbors_enabled = true;
        partition_config.sep_batch_marker = true;
    }

    if (args.sep_batch_marker->count > 0) {
        partition_config.sep_batch_marker = true;
    }

    if (args.bb_ratio->count > 0) {
        partition_config.bb_ratio = args.bb_ratio->ival[0];
    }

    partition_config.run_parallel = args.run_parallel->count > 0;
    if (args.use_queue->count > 0) {
        partition_config.use_queue = parse_bool_option(args.use_queue->sval[0], "--use_queue");
    }

    if (args.restream_include_high_degree_nodes->count > 0) {
        partition_config.restream_include_high_degree_nodes = true;
    }

    if (args.d_max->count > 0) {
        if (args.d_max->ival[0] == -1) {
            partition_config.d_max = std::numeric_limits<LongNodeID>::max();
        } else {
            partition_config.d_max = args.d_max->ival[0];
        }
    }

    if (args.buffer_score->count > 0) {
        if (strcmp("cbs", args.buffer_score->sval[0]) == 0) {
            partition_config.buffer_score_type = BUFFER_SCORE_CBS;
        } else if (strcmp("cbsq", args.buffer_score->sval[0]) == 0) {
            partition_config.buffer_score_type = BUFFER_SCORE_CBSQ;
        } else if (strcmp("anr", args.buffer_score->sval[0]) == 0) {
            partition_config.buffer_score_type = BUFFER_SCORE_ANR;
        } else if (strcmp("haa", args.buffer_score->sval[0]) == 0) {
            partition_config.buffer_score_type = BUFFER_SCORE_HAA;
        } else if (strcmp("cms", args.buffer_score->sval[0]) == 0) {
            partition_config.buffer_score_type = BUFFER_SCORE_CMS;
        } else if (strcmp("nss", args.buffer_score->sval[0]) == 0) {
            partition_config.buffer_score_type = BUFFER_SCORE_NSS;
        } else if (strcmp("gts", args.buffer_score->sval[0]) == 0) {
            partition_config.buffer_score_type = BUFFER_SCORE_GTS;
        } else {
            std::cout << "--b_score should be one of: haa, cbs, cbsq, anr, cms, nss, gts" << '\n';
            exit(0);
        }
    }

    if (args.use_fennel_objective->count > 0) {
        partition_config.use_fennel_objective = true;
    }

    if (args.ram_stream->count > 0) {
        partition_config.ram_stream = true;
    }

    if (args.write_log->count > 0) {
        partition_config.write_log = true;
    }

    if (args.evaluate->count > 0) {
        if (strcmp("true", args.evaluate->sval[0]) == 0 ||
            strcmp("1", args.evaluate->sval[0]) == 0) {
            partition_config.evaluate = true;
        } else if (strcmp("false", args.evaluate->sval[0]) == 0 ||
                   strcmp("0", args.evaluate->sval[0]) == 0) {
            partition_config.evaluate = false;
        } else {
            std::cout << "--evaluate should be one of: true, false, 1, 0" << '\n';
            exit(0);
        }
    }

    if (args.stream_output_progress->count > 0) {
        partition_config.stream_output_progress = true;
    }

    if (args.stream_allow_ghostnodes->count > 0) {
        partition_config.stream_allow_ghostnodes = true;
        partition_config.double_non_ghost_edges = true;
    }

    if (args.stream_initial_bisections->count > 0) {
        partition_config.stream_initial_bisections = true;
    }

    if (args.fennel_dynamics->count > 0) {
        if (strcmp("original", args.fennel_dynamics->sval[0]) == 0) {
            partition_config.fennel_dynamics = FENNELADP_ORIGINAL;
        } else if (strcmp("double", args.fennel_dynamics->sval[0]) == 0) {
            partition_config.fennel_dynamics = FENNELADP_DOUBLE;
        } else if (strcmp("linear", args.fennel_dynamics->sval[0]) == 0) {
            partition_config.fennel_dynamics = FENNELADP_LINEAR;
        } else if (strcmp("quadratic", args.fennel_dynamics->sval[0]) == 0) {
            partition_config.fennel_dynamics = FENNELADP_QUADRATIC;
        } else if (strcmp("midlinear", args.fennel_dynamics->sval[0]) == 0) {
            partition_config.fennel_dynamics = FENNELADP_MID_LINEAR;
        } else if (strcmp("midquadratic", args.fennel_dynamics->sval[0]) == 0) {
            partition_config.fennel_dynamics = FENNELADP_MID_QUADRATIC;
        } else if (strcmp("midconstant", args.fennel_dynamics->sval[0]) == 0) {
            partition_config.fennel_dynamics = FENNELADP_MID_CONSTANT;
        } else if (strcmp("edgecut", args.fennel_dynamics->sval[0]) == 0) {
            partition_config.fennel_dynamics = FENNELADP_EDGE_CUT;
        }
    }

    bool keep_threshold_contract = false;
    if (args.ghost_nodes_procedure->count > 0) {
        if (!args.stream_allow_ghostnodes->count) {
            std::cout
                << "--ghost_nodes_procedure only works together with --stream_allow_ghostnodes"
                << '\n';
            exit(0);
        }
        if (strcmp("contract", args.ghost_nodes_procedure->sval[0]) == 0) {
            partition_config.ghost_nodes_threshold = 0;
        } else if (strcmp("keep", args.ghost_nodes_procedure->sval[0]) == 0) {
            partition_config.ghost_nodes_threshold = std::numeric_limits<LongNodeID>::max();
        } else if (strcmp("keepthresholdcontract", args.ghost_nodes_procedure->sval[0]) == 0) {
            partition_config.ghost_nodes_threshold = 1024;
            keep_threshold_contract = true;
        } else {
            std::cout
                << "--ghost_nodes_procedure should be one of: contract, keep, keepthresholdcontract"
                << '\n';
            exit(0);
        }
    }

    if (args.ghost_nodes_threshold->count > 0) {
        if (!keep_threshold_contract) {
            std::cout << "--ghost_nodes_threshold should only be used together with "
                         "--ghost_nodes_procedure=keepthresholdcontract"
                      << '\n';
            exit(0);
        }
        partition_config.ghost_nodes_threshold =
            static_cast<LongNodeID>(args.ghost_nodes_threshold->dval[0]);
    }

    if (args.num_streams_passes->count > 0) {
        partition_config.restream_vcycle = true;
        partition_config.num_streams_passes = args.num_streams_passes->dval[0];
    }

    if (args.batch_inbalance->count > 0) {
        partition_config.batch_inbalance = args.batch_inbalance->dval[0];
    }

    if (args.restream_vcycle->count > 0) {
        partition_config.restream_vcycle = true;
    }

    if (args.initial_part_multi_bfs->count > 0) {
        if (!(args.stream_initial_bisections->count > 0)) {
            std::cout
                << "--initial_part_multi_bfs cannot be used without --stream_initial_bisections."
                << '\n';
            exit(0);
        }
        partition_config.initial_part_multi_bfs = true;
    }

    if (args.skip_outer_ls->count > 0) {
        partition_config.skip_outer_ls = true;
    }

    if (args.initial_part_fennel->count > 0) {
        if (!(args.stream_initial_bisections->count > 0)) {
            std::cout << "--initial_part_fennel cannot be used without --stream_initial_bisections."
                      << '\n';
            exit(0);
        }
        if (args.initial_part_multi_bfs->count > 0) {
            std::cout
                << "--initial_part_fennel and --initial_part_multi_bfs cannot be used together."
                << '\n';
            exit(0);
        }
        partition_config.initial_partitioning_type = INITIAL_PARTITIONING_FENNEL;
        partition_config.label_iterations = 5;
        partition_config.label_iterations_refinement = 5;
        partition_config.use_balance_singletons = false;
        partition_config.node_ordering = NATURAL_NODEORDERING;
        partition_config.stop_rule = STOP_RULE_MULTIBFS;
    }

    if (args.stream_label_rounds->count > 0) {
        if (strcmp("minimal", args.stream_label_rounds->sval[0]) == 0) {
            partition_config.label_iterations = 5;
            partition_config.label_iterations_refinement = 1;
        } else if (strcmp("normal", args.stream_label_rounds->sval[0]) == 0) {
            partition_config.label_iterations = 5;
            partition_config.label_iterations_refinement = 5;
        } else if (strcmp("high", args.stream_label_rounds->sval[0]) == 0) {
            partition_config.label_iterations = 5;
            partition_config.label_iterations_refinement = 15;
        }
    }

    if (args.automatic_buffer_len->count > 0) {
        if (args.stream_buffer->count > 0) {
            std::cout << "--stream_buffer cannot be used together with --automatic_buffer_len."
                      << '\n';
        }
        partition_config.batch_size = 1024 + partition_config.k * 2500;
    }

    if (args.xxx->count > 0) {
        partition_config.xxx = args.xxx->ival[0];
    }
}

inline void apply_edge_stream_options(const EdgeStreamArgRefs& args, Config& partition_config) {
    auto parse_bool_option = [](const char* value, const char* name) {
        if (strcmp("true", value) == 0 || strcmp("1", value) == 0) {
            return true;
        }
        if (strcmp("false", value) == 0 || strcmp("0", value) == 0) {
            return false;
        }
        std::cout << name << " should be one of: true, false, 1, 0" << '\n';
        exit(0);
    };

    if (partition_config.edge_partition) {
        if (args.use_queue->count > 0) {
            partition_config.use_queue = parse_bool_option(args.use_queue->sval[0], "--use_queue");
        }

        if (args.dynamic_alpha->count > 0) {
            partition_config.dynamic_alpha = true;
            partition_config.batch_alpha = false;
        }

        if (args.batch_alpha->count > 0) {
            partition_config.batch_alpha = true;
            partition_config.dynamic_alpha = false;
        }

        if (args.minimal_mode->count > 0) {
            partition_config.minimal_mode = true;
        }

        if (args.num_split_edges->count > 0) {
            partition_config.num_split_edges = args.num_split_edges->ival[0];
        }

        if (args.past_subset_size->count > 0) {
            partition_config.past_subset_size = args.past_subset_size->ival[0];
        }

        if (args.reps->count > 0) {
            partition_config.reps = args.reps->dval[0];
        }
    }

    if (args.benchmark->count > 0) {
        partition_config.benchmark = true;
    }

    if (args.light_evaluator->count > 0) {
        partition_config.light_evaluator = true;
    }

    if (args.label_propagation_iterations_refinement->count > 0) {
        partition_config.label_iterations_refinement =
            args.label_propagation_iterations_refinement->ival[0];
    }
}

} // namespace cli::stream_parameters


#endif /* STREAM_PARAMETER_APPLY_H_ */
