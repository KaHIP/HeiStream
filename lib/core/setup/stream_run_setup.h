/****************************************************************************
 * stream_run_setup.h
 *****************************************************************************/

#ifndef CORE_STREAM_RUN_SETUP_H_
#define CORE_STREAM_RUN_SETUP_H_

#include <string>

#include "partition/partition_config.h"


namespace core::setup {

/**
 * Framework-level preflight that normalizes shared stream-run state.
 *
 * @param config Mutable global runtime config.
 * @param graph_filename Input graph path.
 * @return true when framework preflight succeeds; false otherwise.
 */
bool prepare_framework_run(Config& config, const std::string& graph_filename);

} // namespace core::setup


#endif /* CORE_STREAM_RUN_SETUP_H_ */
