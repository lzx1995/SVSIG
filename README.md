
##  2. Getting Started

### 2.1 Core Organization

The `core/graphBolt/` folder contains the [SVSIG Engine](#3-SVSIG-engine), and Graphbolt's [Stream Ingestor](#5-stream-ingestor) module. The application/benchmark codes (e.g., PageRank) can be found in the `apps/` directory. Useful helper files for generating the stream of changes (`tools/generators/streamGenerator.C`), creating the graph inputs in the correct format (`tools/converters/SNAPtoAdjConverter.C` - from Ligra's codebase), and comparing the output of the algorithms (`tools/output_comparators/`) are also provided.

### 2.2 Requirements
- g++ >= 5.3.0 with support for Cilk Plus.
- [Mimalloc](https://github.com/microsoft/mimalloc) - A fast general purpose memory allocator from Microsoft (version >= 1.6).
    - Use the helper script `install_mimalloc.sh` to install mimalloc.
    - Update the LD_PRELOAD enviroment variable as specified by install_mimalloc.sh script.

**Important: SVSIG requires mimalloc to function correctly and efficiently.**

Note: gcc-5 and gcc-7 come with cilk support by default. You can easily maintain multiple versions of gcc using `update-alternatives` tool. If you currently have gcc-9, you can easily install gcc-5 and switch to it as follows:
```bash
$   # Install gcc-5
$   sudo apt install gcc-5
$   # Set the path for all gcc versions
$   sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 50
$   # gcc-9 version
$   sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 60
$   # Configure gcc to use gcc-5
$   sudo update-alternatives --config gcc
$   # Verify gcc version
$   gcc --version
```

### 2.3 Compiling and Running the Application

Compilation is done from within `apps` directory. To compile, run
Compilation is done from within apps directory. For Tiny batch mutation(<1000), to compile, run
```
cd apps
make clean
make MULTIAGGREGATIONDEVIATION=1 TINYMODE=1  -j
```
Fro Large batch mutation, to compile, run
```
cd apps
make clean
make MULTIAGGREGATIONDEVIATION=1 -j
```
 The executable takes the following command-line parameters:
 - `-s` : Optional parameter to indicate a symmetric (undirected) graph is used. 
 - `-streamPath` : Path to the input stream file or pipe (More information on the input format can be found in [Section 2.4](#24-graph-input-and-stream-input-format)).
 - `-numberOfUpdateBatches` : Optional parameter to specify the number of edge updates to be made. Default is 1.
 - `-nEdges` : Number of edge operations to be processed in a given update batch.
 - `-outputFile` : Optional parameter to print the output of a given algorithms.
 - Input graph file path (More information on the input format can be found in [Section 2.4](#24-graph-input-and-stream-input-format)).

For example,
```bash
$   # Ensure that LD_PRELOAD is set as specified by the install_mimalloc.sh
$   ./PageRank -numberOfUpdateBatches 2 -nEdges 1000 -streamPath ../inputs/sample_edge_operations.txt -outputFile /tmp/output/pr_output ../inputs/sample_graph.adj
$   ./LabelPropagation -numberOfUpdateBatches 3 -nEdges 2000 -streamPath ../inputs/sample_edge_operations.txt -seedsFile ../inputs/sample_seeds_file -outputFile /tmp/output/lp_output ../inputs/sample_graph.adj
$   ./COEM -s -numberOfUpdateBatches 3 -nEdges 2000 -streamPath ../inputs/sample_edge_operations.txt -seedsFile ../inputs/sample_seeds_file -partitionsFile ../inputs/sample_partitions_file -outputFile /tmp/output/coem_output ../inputs/sample_graph.adj
```
Other additional parameters may be required depending on the algorithm. Refer to the `Compute()` function in the application code (`apps/PageRank.C`,etc.) for the supported arguments. Additional configurations for the graph ingestor and the graph can be found in [Section 5](#4-stream-ingestor).

### 2.4 Graph Input and Stream Input Format

The initial input graph should be in the [adjacency graph format](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html). 
For example, the SNAP format (edgelist) and the adjacency graph format for a sample graph are shown below.

SNAP format:
```txt
0 1
0 2
2 0
2 1
```
 Adjacency Graph format:
```txt
AdjacencyGraph
3
4
0
2
2
1
2
0
1
```
You can use `tools/converters/SNAPtoAdjConverter` to convert an input graph in Edgelist format (SNAP format) to the adjacency graph format, as follows:
```bash
$   ./SNAPtoAdjConverter inputGraph.snap inputGraph.adj
$   # for undirected (symmetric) graphs, use the -s flag
$   ./SNAPtoAdjConverter -s inputGraph.snap inputGraphUndirected.adj 
```
The streaming input file should have the edge operation (addition/deletion) on a separate line. The edge operation should be of the format, `[d/a] source destination` where `d` indicates edge deletion and `a` indicates edge addition. Example streaming input file:
```bash
a 1 2
d 2 3
a 4 5
...
```

Edge operations can be streamed through a pipe using `tools/generators/streamGenerator.C`. It takes in the following command-line parameters:
- `-edgeOperationsFile` : Input file containing the edge operations in the format mentioned above.
- `-outputPipe` : Path of the output pipe where the edges are streamed to.

```bash
$   cd tools/generators
$   make streamGenerator
$   ./streamGenerator -edgeOperationsFile ../inputs/sample_edge_operations.txt -outputPipe ../inputs/sample_edge_operations.pipe
```
More details regarding the ingestor can be found in [Section 4](#4-stream-ingestor).

## 4. Stream Ingestor

The stream ingestor FIFO is specified by `-streamPath`. Edge operations can be written to this FIFO. `-nEdges` specifies the maximum number of edge operations that can be passed to the GraphBolt engine in a single batch. The GraphBolt engine will continue to receive batches of edges from the stream ingestor until either the stream is closed (when there are no more writers to the FIFO) or when `-numberOfBatches` has been exceeded. If the writing end of the FIFO is not opened, the GraphBolt engine (which is the reading end) will block and wait until it is opened. 

There are a few optional flags that can affect the behaviour and determine the validity of the edge operations  passed to the command line parameter `-streamPath`:

- `-fixedBatchSize`: Optional flag to ensure that the batch size is strictly adhered to. If the FIFO does not contain enough edges, the ingestor will block until it has received enough edges specified by `-nEdges` or until the stream is closed. 
- `-enforceEdgeValidity`: Optional flag to ensure that all edge operations in the batch are valid. For example, an edge deletion operation is valid only if the edge to be deleted is present in the graph. In the case of a `simple graph` (explained below), an edge addition operation is valid only if that edge does not currently exist in the graph. Invalid edges are discarded and are not included while counting the number of edges in a batch.
- `-simple`: Optional flag used to ensure that the input graph remains a simple graph (ie. no duplicate edges). The input graph is checked to remove all duplicate edges. Duplicate edges are not allowed within a batch and edge additions are checked to ensure that the edge to be added does not yet exist within the graph.


## 7. Acknowledgements
Some utility functions from [GraphBolt](https://github.com/pdclab/graphbolt), [Ligra](https://github.com/jshun/ligra) and [Problem Based Benchmark Suite](http://www.cs.cmu.edu/~pbbs/index.html) are used as part of this project. We are thankful to them for releasing their source code.
