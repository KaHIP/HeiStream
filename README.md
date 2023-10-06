# HeiStream 

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

## Installation Notes

### Requirements

* C++-14 ready compiler 
* CMake 
* Scons (http://www.scons.org/)
* Argtable (http://argtable.sourceforge.net/)

### Building HeiStream

To build the software, run
```shell
./compile.sh
```

Alternatively, you can use the standard CMake build process.

The resulting binary is located in `deploy/heistream`.

## Running HeiStream

To partition a graph in METIS format using the basic model of HeiStream, run

```shell
./heistream <graph filename> --k=<number of blocks> --stream_buffer=<nodes per bufer>
```

To partition a graph in METIS format using the extended model of HeiStream, run

```shell
./heistream <graph filename> --k=<number of blocks> --stream_buffer=<nodes per bufer> --stream_allow_ghostnodes
```

For a complete list of parameters alongside with descriptions, run:

```shell
./heistream --help
```

For a description of the graph format, please have a look at the [KaHiP manual](https://github.com/KaHIP/KaHIP/raw/master/manual/kahip.pdf).

## Edge Balancing

HeiStream was not designed to balance edges, but we have implemented a temporary solution to allow this. 
If you want to balance edges instead of nodes, you can enable the --balance_edges flag within your command for executing HeiStream.

## Licensing

HeiStream is free software provided under the MIT License.
If you publish results using our algorithms, please acknowledge our work by citing the following paper ([pdf](https://dl.acm.org/doi/10.1145/3546911)):

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

