cmake_minimum_required(VERSION 3.10)
include(CheckCXXCompilerFlag)
include(FetchContent)
include(FetchContent)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_program(CCACHE_PROGRAM ccache)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()
project(KaHIP C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Optional profiling/timing instrumentation.
option(ENABLE_TIME_MEASUREMENTS "Enable fine-grained scoped timing instrumentation" OFF)
if(ENABLE_TIME_MEASUREMENTS)
  add_definitions(-DENABLE_TIME_MEASUREMENTS)
endif()

# if no build mode is specified build in release mode
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# tweak compiler flags
CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
  add_definitions(-funroll-loops)
endif()
# Disabled globally to keep clang-tidy toolchain-compatible.
# GCC builds remain unaffected for this project without this flag.
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
if(COMPILER_SUPPORTS_WALL)
  add_definitions(-Wall)
endif()
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  add_definitions(-march=native)
endif()
CHECK_CXX_COMPILER_FLAG(-fpermissive COMPILER_SUPPORTS_FPERMISSIVE)
if(COMPILER_SUPPORTS_FPERMISSIVE)
  add_definitions(-fpermissive)
endif()

# check dependencies
find_package(MPI REQUIRED)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP support detected")
  add_definitions(${OpenMP_CXX_FLAGS})
else()
  message(WARNING "OpenMP not available, activating workaround")
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
endif()

# readerwriterqueue
FetchContent_Declare(
  readerwriterqueue
  GIT_REPOSITORY    https://github.com/cameron314/readerwriterqueue
  GIT_TAG           master
)

FetchContent_MakeAvailable(readerwriterqueue)

# 64 Bit option
#option(64BITMODE "64 bit mode" OFF)
option(64BITMODE "64 bit mode" ON)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()

# optimized output
option(OPTIMIZED_OUTPUT "optimized output" OFF)
if(OPTIMIZED_OUTPUT)
  add_definitions("-DKAFFPAOUTPUT")
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.0.3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/flatbuffers/include)
include_directories(${readerwriterqueue_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/sparsehash/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/io)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/partition)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/tools)

# link flatbuffers to log experimental output
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/flatbuffers
        ${CMAKE_CURRENT_BINARY_DIR}/flatbuffers-build
        EXCLUDE_FROM_ALL)





set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltcmalloc -lpthread")


set(LIBKAFFPA_SOURCE_FILES
  lib/data_structure/graph_hierarchy.cpp
  lib/graph/algorithms/strongly_connected_components.cpp
  lib/graph/algorithms/topological_sort.cpp
  lib/graph/algorithms/push_relabel.cpp
  lib/io/graph_io.cpp
  lib/partition/model/model_graph_ops.cpp
  lib/io/stream_reader.cpp
  lib/io/output/flatbuffer_log_writer.cpp
  lib/io/output/partition_output_writer.cpp
  lib/io/edge_stream_reader.cpp
  lib/partition/assignment/edge_stream_assign.cpp
  lib/partition/metrics/edge_partitioning_quality_eval.cpp
  lib/algorithms/edge_partitioning/model/edge_batch_model_builder.cpp
  lib/partition/assignment/node_stream_assign.cpp
  lib/partition/metrics/node_partitioning_quality_eval.cpp
  lib/algorithms/node_partitioning/model/node_batch_model_builder.cpp
  lib/io/node_stream_reader.cpp
  lib/algorithms/node_partitioning/execution/node_execution_common.cpp
  lib/algorithms/node_partitioning/execution/node_execution_priority_buffer_parallel.cpp
  lib/algorithms/node_partitioning/execution/node_execution_priority_buffer_sequential.cpp
  lib/algorithms/node_partitioning/execution/node_execution_no_buffer.cpp
  lib/algorithms/node_partitioning/execution/node_execution_direct.cpp
  lib/algorithms/node_partitioning/setup/node_partitioning_setup.cpp
  lib/algorithms/node_partitioning/mode/node_partitioning_mode_resolver.cpp
  lib/algorithms/node_partitioning/mode/node_partitioning_mode_prepare.cpp
  lib/algorithms/node_partitioning/scheduler/node_partitioning_scheduler.cpp
  lib/algorithms/node_partitioning/evaluation/node_partitioning_evaluator.cpp
  lib/algorithms/node_partitioning/reporting/node_partitioning_reporter.cpp
  lib/algorithms/node_partitioning/node_partitioning_algorithm.cpp
  lib/algorithms/node_partitioning/factory/node_partitioning_algorithm_factory.cpp
  lib/algorithms/bootstrap/register_builtin_algorithms.cpp
  lib/core/factory/algorithm_registry.cpp
  lib/core/setup/stream_run_setup.cpp
  lib/core/pipeline/algorithm_pipeline.cpp
  lib/core/pipeline/stream_engine.cpp
  lib/algorithms/edge_partitioning/setup/edge_partitioning_setup.cpp
  lib/algorithms/edge_partitioning/mode/edge_partitioning_mode_resolver.cpp
  lib/algorithms/edge_partitioning/scheduler/edge_partitioning_scheduler.cpp
  lib/algorithms/edge_partitioning/evaluation/edge_partitioning_evaluator.cpp
  lib/algorithms/edge_partitioning/reporting/edge_partitioning_reporter.cpp
  lib/algorithms/edge_partitioning/edge_partitioning_algorithm.cpp
  lib/algorithms/edge_partitioning/factory/edge_partitioning_algorithm_factory.cpp
  lib/partition/model/model_partition_mapper.cpp
  lib/partition/stages/stream_partition_stages.cpp
  lib/partition/stages/batch_partition_execution.cpp
  lib/tools/stream_util.cpp
  lib/algorithms/edge_partitioning/execution/edge_execution.cpp
  lib/tools/quality_metrics.cpp
  lib/tools/random_functions.cpp
  lib/tools/graph_extractor.cpp
  lib/tools/misc.cpp
  lib/tools/partition_snapshooter.cpp
  lib/partition/graph_partitioner.cpp
  lib/partition/coarsening/coarsening.cpp
  lib/partition/coarsening/contraction.cpp
  lib/partition/coarsening/edge_rating/edge_ratings.cpp
  lib/partition/coarsening/matching/matching.cpp
  lib/partition/coarsening/matching/random_matching.cpp
  lib/partition/coarsening/matching/gpa/path.cpp
  lib/partition/coarsening/matching/gpa/gpa_matching.cpp
  lib/partition/coarsening/matching/gpa/path_set.cpp
  lib/partition/coarsening/clustering/node_ordering.cpp
  lib/partition/coarsening/clustering/size_constraint_label_propagation.cpp
  lib/partition/initial_partitioning/initial_partitioning.cpp
  lib/partition/initial_partitioning/initial_partitioner.cpp
  lib/partition/initial_partitioning/initial_refinement/initial_refinement.cpp
  lib/partition/initial_partitioning/init_fennel.cpp
  lib/partition/uncoarsening/uncoarsening.cpp
  lib/partition/uncoarsening/refinement/mixed_refinement.cpp
  lib/partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.cpp
  lib/partition/uncoarsening/refinement/refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/partial_boundary.cpp
  extern/argtable3-3.0.3/argtable3.c)
add_library(libkaffpa OBJECT ${LIBKAFFPA_SOURCE_FILES})

# generate targets for each binary
add_executable(heistream app/streampartition.cpp $<TARGET_OBJECTS:libkaffpa> )
target_compile_definitions(heistream PRIVATE "-DMODE_STREAMPARTITION")
target_link_libraries(heistream ${OpenMP_CXX_LIBRARIES} flatbuffers readerwriterqueue)
#target_link_libraries (streampartition /usr/lib/x86_64-linux-gnu/libtcmalloc_and_profiler.so)
install(TARGETS heistream DESTINATION bin)

add_executable(heistream_edge app/streamedgepartition.cpp $<TARGET_OBJECTS:libkaffpa> )
target_compile_definitions(heistream_edge PRIVATE "-DMODE_STREAMPARTITION")
target_link_libraries(heistream_edge ${OpenMP_CXX_LIBRARIES} flatbuffers)
#target_link_libraries(heistream_edge flatbuffers)
install(TARGETS heistream_edge DESTINATION bin)
