# Adding a New Streaming Algorithm

HeiStream is structured as a framework for streaming graph algorithms.
Algorithm-specific behavior lives under `lib/algorithms/<name>/`, while shared
infrastructure stays in `lib/core`, `lib/io`, and `lib/partition`.

## Design Principles

1. Keep `lib/core` framework-only.
- `core` owns lifecycle (`prepare_run` -> `run` -> `finalize_run`), registry,
  orchestration, and timing primitives.
- `core` must not implement algorithm-specific logic.

2. Keep algorithm logic in `lib/algorithms/<algorithm_name>/`.
- Each algorithm owns setup, execution scheduling, mode resolution, model,
  evaluation, reporting, and assembly.

3. Keep cross-algorithm partition primitives in `lib/partition`.
- Generic assignment, heuristics, state, setup helpers, model helpers, and
  shared partition stages belong here.
- Algorithm-specific model decisions stay in the algorithm folder.

4. Keep I/O reusable in `lib/io`.
- Streaming input follows a common pattern: read header once, then read lines
  incrementally.
- Output writers stay generic (`partition_output_writer`,
  `flatbuffer_log_writer`).

## Current Repository Structure

- `app/`
- Entry points only: `streampartition.cpp`, `streamedgepartition.cpp`

- `lib/cli/`
- CLI parsing and configuration population

- `lib/core/`
- `context/`: shared run context
- `factory/`: algorithm registry
- `interfaces/`: setup/scheduler/evaluator/reporter/algorithm contracts
- `pipeline/`: engine and pipeline execution
- `setup/`: framework-level pre-run checks
- `timing/`: scoped timing utilities

- `lib/io/`
- `stream_reader.*`: generic header + line-by-line stream reading
- `node_stream_reader.*`, `edge_stream_reader.*`: stream adapters used by
  current algorithms
- `output/`: partition/flatbuffer output writers

- `lib/algorithms/`
- `bootstrap/`: registration of built-in algorithms
- `node_partitioning/`
- `edge_partitioning/`
- Each algorithm folder has:
- `setup/`, `mode/`, `scheduler/`, `execution/`, `model/`, `evaluation/`,
  `reporting/`, `factory/`, and `<algorithm>_algorithm.*`

- `lib/partition/`
- Shared partitioning infrastructure:
- `assignment/`, `heuristics/`, `metrics/`, `model/`, `setup/`, `stages/`,
  `state/`, plus KaHIP multilevel components used by batch solving

- `lib/data_structure/`, `lib/tools/`, `lib/graph/`
- Reusable low-level data structures and utilities

## What to Implement for a New Algorithm

1. Add a new stream mode.
- Edit `lib/definitions.h` and add `StreamMode::<YourMode>`.

2. Create a new algorithm folder.
- Create `lib/algorithms/<your_algorithm>/` with:
- `<your_algorithm>_algorithm.h/.cpp`
- `setup/`
- `mode/`
- `scheduler/`
- `execution/`
- `model/`
- `evaluation/`
- `reporting/`
- `factory/`

3. Implement framework interfaces.
- `core::interfaces::IAlgorithmSetup`
- `core::interfaces::IAlgorithmScheduler`
- `core::interfaces::IAlgorithmEvaluator`
- `core::interfaces::IAlgorithmReporter`

4. Assemble and expose factory.
- Provide:
- `std::unique_ptr<core::interfaces::IStreamAlgorithm> create_algorithm();`
- Wire all interface implementations inside `factory/<...>_algorithm_factory.cpp`.

5. Register in bootstrap.
- Edit `lib/algorithms/bootstrap/register_builtin_algorithms.cpp`.
- Register `StreamMode::<YourMode>` with your factory.

6. Add build wiring.
- Add your new `.cpp` files to `LIBKAFFPA_SOURCE_FILES` in `CMakeLists.txt`.

7. Add CLI mode wiring.
- Update CLI parsing to map the user-facing mode choice to
  `StreamMode::<YourMode>`.

## Component Responsibilities

- `setup/`
- Run-scoped defaults and per-pass initialization after header parsing.

- `mode/`
- Resolve execution variant from config. Keep it pure when possible.

- `scheduler/`
- Select and execute the concrete execution path for each pass/batch.

- `execution/` 
- Main stream loop.

- `model/`
- Build algorithm-specific in-memory graph/model representation if valid.

- `evaluation/`
- Compute quality metrics when evaluation is enabled.

- `reporting/`
- Print run configuration and run summary.
- Emit optional output artifacts (partition files, flatbuffer logs).

## Dependency Rules

- `lib/core` cannot depend on algorithm folders.
- `lib/io` should remain algorithm-agnostic.
- Algorithms may depend on `core`, `io`, `partition`, `tools`,
  `data_structure`.
