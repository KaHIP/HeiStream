# HeiStream & HeiStreamE & BuffCut
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/315a384a931d4136a89a4b2846ae0475)](https://app.codacy.com/gh/KaHIP/HeiStream/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

HeiStream is a buffered streaming algorithm to heuristically solve the graph partitioning problem: dividing the nodes of a graph into k disjoint blocks of roughly the same size while minimizing the number of edges running between blocks.
HeiStream is a first attempt to close a gap observed in the space of available partitioning algorithms. 
On the one hand, there are streaming algorithms that have been adopted to partition massive graph data on small machines. 
In the streaming model, vertices arrive one at a time including their neighborhood and then have to be assigned directly to a block. 
These algorithms can partition huge graphs quickly with little memory, but they produce partitions with low solution quality. 
On the other hand, there are offline (shared-memory) multilevel algorithms that produce partitions with high quality but also need a machine with enough
memory to partition huge networks. 
HeiStream uses a buffered streaming computational model and a multilevel algorithm, which allows it to compute significantly improved partitions of huge graphs using a single machine with little memory.

<p align="center">
<img src="./img/MultilevelBufferedStreamGP.png"
  alt="Overall Structure of HeiStream"
  width="601" >
</p>

The image above illustrates the overall structure of HeiStream. 
It slides through the streamed graph G by repeating the following successive operations until all the nodes of G are assigned to blocks. 
First, it loads a batch containing the desired number of nodes alongside with their adjacency lists. 
Second, it builds a model B to be partitined. 
This model represents the already partitioned vertices as well as the nodes of the current batch. 
Third, it partitions B with a multilevel partitioning algorithm to optimize for the Fennel objective function. 
Finally, it permanently assigns the nodes from the current batch to blocks. 

HeiStream also integrates BuffCut directly into this node-streaming pipeline.
When priority buffering is enabled, the stream first inserts candidate nodes into a bounded priority buffer `Q`; high-priority nodes are then selected from `Q` to form each batch `B` before model construction and multilevel partitioning. 
This improves stream-order robustness while reusing the same core HeiStream model-building and partitioning framework.

### HeiStreamE

HeiStreamE extends HeiStream to solve the edge partitioning problem, dividing the edges of a graph into k disjoint blocks while minimizing vertex replicas.

<p align="center">
  <img src="./img/HeiStreamEOverall.png" alt="Overall Structure of HeiStreamE" width="601">
</p>

We slide through the input graph G by iteratively performing the following series of operations until all the edges of G are assigned to blocks.
First, we load a batch composed of δ vertices and their associated adjacency lists, thereby obtaining a subgraph G_b contained within the graph G.
This operation yields edges connecting vertices within the current batch, and edges connecting vertices in the current batch to vertices streamed in previous batches.
Second, we build a model β_b corresponding to G_b, where the edges of G_b are transformed into vertices. Additionally, we incorporate a representation of block assignments
from previous batches into β_b.
Third, we partition β_b using a multilevel vertex partitioning algorithm.
We conclude by permanently assigning the edges in G that correspond to vertices in our model β_b to their respective blocks.

Our experiments demonstrate that HeiStreamE produces the best solution quality (replication factor) among all known (re)streaming edge partitioners and is highly memory efficient.

## Installation Notes

### Requirements

- C++17 compiler (e.g., `g++`/`clang++`)
- CMake
- MPI development package
- OpenMP support in the compiler toolchain
- Git (needed by CMake `FetchContent` dependencies)

### Building HeiStream & HeiStreamE

To build the software, run
```shell
./compile.sh
```

Profiling/timing instrumentation can be toggled at build time:
```shell
./compile.sh OFF   # default: disable fine-grained scoped timers
./compile.sh ON    # enable fine-grained scoped timers
```

Alternatively, you can use the standard CMake build process.

The resulting binary is located in `deploy/heistream` and `deploy/heistream_edge`.

## Node Partitioning Modes and Restreaming

HeiStream node partitioning has three execution behaviors:

- No priority buffer:
  - `--buffer_size=0` (default)
  - Streamed nodes are processed in batches of `--batch_size`.
- Priority buffer (sequential):
  - `--buffer_size > 1`
  - BuffCut integration in HeiStream is enabled.
- Priority buffer (parallel):
  - `--buffer_size > 1 --run-parallel`
  - Uses the parallel priority-buffer pipeline.

Restreaming (`--num_streams_passes > 1`) is supported in all node partitioning modes:

- Pass 1 follows the selected mode above.
- Additional passes re-read the stream and refine assignments.

### Priority-Buffer Algorithm Details

When `--buffer_size > 1`, HeiStream uses the integrated BuffCut priority-buffer logic on top of batch-wise multilevel partitioning.

- It builds:
  - A bounded priority buffer `Q` with capacity `buffer_size`.
  - A batch `B` with target size `batch_size`.
- Stream Assignment Order:
  - Very high-degree nodes (controlled by threshold `--d_max`) are assigned immediately with Fennel scoring.
  - Remaining nodes are inserted into `Q` with a buffer score (`--b_score`, default `haa`).
- Eviction and batching:
  - When `Q` fills, highest-priority nodes are extracted into `B`.
  - Once `|B|` reaches `batch_size`, HeiStream builds a batch model graph and runs multilevel partitioning.
  - Assignments are generalized back to global stream state, and processing continues.
- Restreaming:
  - Additional passes (`--num_streams_passes > 1`) re-read the stream and refine assignments.
- Parallel mode:
  - `--run-parallel` enables the parallel priority-buffer pipeline for node partitioning.
  - The same high-level logic is preserved, but I/O, buffering, and partition tasks are overlapped.

## Timing and Reporting

- `Total time` is always reported in the run summary.
- Fine-grained timing components (`IO/Model/Postprocess/Buffer maintenance/Partition`) are reported only when compiled with profiling enabled (`./compile.sh ON`).
- `--evaluate` controls quality evaluation and summary printing:
  - default: `--evaluate=true`
  - `--evaluate=false` skips evaluator work and summary reporting.
- `--write_log` controls writing the flatbuffer `.bin` log file.

## Running HeiStream

To partition a graph in METIS format using baseline HeiStream (default behavior recommended when stream order is not random/adversarial), run:

```shell
./deploy/heistream <graph filename> --k=<number of blocks> [--batch_size=<nodes per batch>]
```

Priority buffering is enabled when `--buffer_size > 1`.

```shell
./deploy/heistream <graph filename> --k=<number of blocks> --batch_size=<nodes per batch> --buffer_size=<pq size>
```

To run parallel priority-buffer mode, add `--run-parallel`:

```shell
./deploy/heistream <graph filename> --k=<number of blocks> --batch_size=<nodes per batch> --buffer_size=<pq size> --run-parallel
```

For a complete list of parameters alongside with descriptions, run:

```shell
./deploy/heistream --help
```

For a description of the graph format, please have a look at the [KaHiP manual](https://github.com/KaHIP/KaHIP/raw/master/manual/kahip.pdf).

## Edge Balancing

HeiStream was not designed to balance edges, but we have implemented a temporary solution to allow this. 
If you want to balance edges instead of nodes, you can enable the --balance_edges flag within your command for executing HeiStream.

## Running HeiStreamE

To partition a graph in the METIS format using HeiStreamE, run

```shell
./deploy/heistream_edge <graph filename> --k=<number of blocks> [--batch_size=<delta vertices per buffer>]
```
By default, this command creates and writes a vector of partition IDs for all edges of size _m_ to a text file.

To avoid the creation of this vector, use the flag `--stream_output_progress`. This writes the partition IDs assigned to edges within each batch on the fly. This avoids the _O(m)_ overhead. 

To only benchmark the algorithm for runtime and memory consumption, use the flag `--benchmark`. With this command, partition IDs are not written to a file.

You may also use a light version of the in-built evaluator for partition results, which is more efficient when the graph is very large and cannot easily fit into memory using the flag `--light_evaluator`.

Additionally, you have the option to specify the name of the file to which partition IDs will be written with the flag `--output_filename=<name>`. By default, the output file is called `tmp_output.txt`.

For example, to test the provided graph(s) in `examples/`, you can use the following commands:
```shell
./deploy/heistream_edge examples/com-amazon.graph --k=32 --batch_size=32768
```
```shell
./deploy/heistream_edge examples/com-amazon.graph --k=4 --batch_size=16384 --stream_output_progress
```
```shell
./deploy/heistream_edge examples/com-amazon.graph --k=16 --batch_size=131072 --stream_output_progress --benchmark --output_filename=com-amazon_16_131072_part.txt
```

For a complete list of parameters alongside with descriptions, run:

```shell
./deploy/heistream_edge --help
```

### Edge Partitioner Restreaming

Edge partitioner mode also supports multiple passes via `--num_streams_passes`, with one important restriction:

- Restreaming in edge mode is currently supported for `minimal_mode` only.
- If `--num_streams_passes > 1` is used with non-minimal edge mode, execution aborts with an explicit message.

The default edge partitioner setup in this repository uses:

- `minimal_mode=true`
- batch alpha for Fennel alpha updates

## Licensing

HeiStream and HeiStreamE are free software provided under the MIT License.
If you publish results using our algorithms, please acknowledge our work by citing the following papers:

```
@article{heistream,
	author = {Marcelo Fonseca Faraj and Christian Schulz},
	title = {Buffered Streaming Graph Partitioning},
	year = {2022},
	publisher = {Association for Computing Machinery},
	address = {New York, NY, USA},
	issn = {1084-6654},
	url = {https://doi.org/10.1145/3546911},
	doi = {10.1145/3546911},
	journal = {ACM Journal of Experimental Algorithmics},
	month = {jun},
	keywords = {streaming algorithms, multilevel algorithms, graph partitioning}
}
```

```
@InProceedings{HeiStreamE,
  author =	{Chhabra, Adil and Fonseca Faraj, Marcelo and Schulz, Christian and Seemaier, Daniel},
  title =	{{Buffered Streaming Edge Partitioning}},
  booktitle =	{22nd International Symposium on Experimental Algorithms (SEA 2024)},
  pages =	{5:1--5:21},
  series =	{Leibniz International Proceedings in Informatics (LIPIcs)},
  ISBN =	{978-3-95977-325-6},
  ISSN =	{1868-8969},
  year =	{2024},
  volume =	{301},
  editor =	{Liberti, Leo},
  publisher =	{Schloss Dagstuhl -- Leibniz-Zentrum f{\"u}r Informatik},
  address =	{Dagstuhl, Germany},
  URL =		{https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SEA.2024.5},
  URN =		{urn:nbn:de:0030-drops-203701},
  doi =		{10.4230/LIPIcs.SEA.2024.5},
  annote =	{Keywords: graph partitioning, edge partitioning, streaming, online, buffered partitioning}
}
```

Note:
- The program writes a [flatbuffer](https://github.com/google/flatbuffers) `.bin` log only if you pass `--write_log`.
  Node streaming uses `graph_k_batch.bin` in baseline mode and `graph_k_batch_buffer.bin` when `buffer_size > 1`.
- `--evaluate=true` (default) prints run summaries and quality metrics.
  Use `--evaluate=false` to skip evaluator work.
- For a description of the METIS graph format, please have a look at the [KaHiP manual](https://github.com/KaHIP/KaHIP/raw/master/manual/kahip.pdf).
- The binary equivalents of the METIS graphs were obtained by using the converter provided in
  [KaGen](https://github.com/KarlsruheGraphGeneration/KaGen) with `output-format=parhip`.
