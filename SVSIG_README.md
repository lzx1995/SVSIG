# Artifact Evaluation

1. Overview 
2. Setting up GraphBolt 
    - 2.1 Requirements 
    - 2.2 Compiling Applications 
    - 2.3 Dataset Details 
        - 2.3.1 Preparing Streaming Datasets for any Graph (Advice to Reviewer: Skip this time-consuming step by directly using provided Wiki-Vote dataset)
3. Running Applications 
4. Experiments in Paper  
    - 4.1 Evaluating the System
        - 4.1.1 Understanding the Output 
        - 4.1.2 Varying Batch Sizes 
        - 4.1.4 Number of Edges Processed 
        - 4.1.5 Sensitivity Experiments 


## 1. Overview

This readme provides instructions for reproducing the experiments in our paper, `Source Vertex Suppression Incremental Graph Processing for Streaming Graph`. 

The source vertex suppression incremental processing technique is implemented in the DZIG runtime in order to retain efficiency in presence of sparse computations, thereby pushing the boundary of dependency-driven processing of streaming graphs.

The DZIG runtime offers several different processing capabilities (e.g., different modes to ingest streaming graphs), details of which are not relevant for artifact evaluation. Hence, this readme only focuses on the necessary parts to help evaluate the artifact.

##  2. Setting up SVSIG

### 2.1 Requirements
- g++ >= 5.3.0 with support for Cilk Plus (Note: gcc-5 and gcc-7 come with cilk support by default.)
- cmake
- [Mimalloc](https://github.com/microsoft/mimalloc) - A fast general purpose memory allocator from Microsoft (version 1.6).
    - Use the helper script `install_mimalloc.sh` to install mimalloc.
      ```bash
      sh install_mimalloc.sh
      ```
    - Update the LD_PRELOAD enviroment variable as specified by install_mimalloc.sh script.

**Important: SVSIG requires mimalloc to function correctly and efficiently.**

### 2.2 Compiling Applications

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

### 2.3 Dataset Details
  The input graphs used in evaluation are: [UK](http://konect.cc/networks/dimacs10-uk-2007-05/), [TW](http://konect.cc/networks/twitter/), [TT](http://konect.cc/networks/twitter_mpi/) and [FT](http://konect.cc/networks/friendster/). For detailed instructions on converting these datasets to the `graph.snap` format, please refer to [DATASET_CONV.md](DATASET_CONV.md).
  
[Section 2.3.1](#231-preparing-streaming-datasets-for-any-graph) gives the steps to prepare correct streaming inputs from these datasets. 

  All the above input graphs are very large in size (billions of edges), and converting them to the appropriate format and creating input streams is a time consuming process. To simplify the artifact evaluation process, we have provided the input streams for a smaller graph called [Wiki-Vote](https://snap.stanford.edu/data/wiki-Vote.html) in the `inputs/wiki_vote/`. Hence, you can skip section 2.3.1 and directly start with [section 3](#3-running-applications-1-compute-second-x-5-applications--5-compute-seconds) with the Wiki-Vote graph. 

#### 2.3.1 Preparing Streaming Datasets for any Graph

The graph files are obtained in the SNAP format (edge list) format (say, `graph.snap`). In order to create the datasets for our evaluation, do the following:
1. Distribute the lines (excluding any comments) in the graph file to 2 different files as follows:
  ```bash
  # Remove any comments in the graph.snap file
  # Obtain the line numbers in the snap file and use it to divide it into 2 files
  wc -l graph.snap 
  head -n 1000 graph.snap > initial_graph.snap
  tail -n 1000+1 graph.snap > additions.snap
  ```
2. Convert the snap file to adjacency list format as follows
  ```bash
  cd tools/converters
  make SNAPtoAdjConverter
  ./SNAPtoAdjConverter initial_graph.snap initial_graph.adj
  # for undirected (symmetric) graphs, use the -s flag
  ./SNAPtoAdjConverter -s initial_graph.snap initial_graph.adj.un
  ```
  More details in section 2.4 of [README.md](README.md#24-graph-input-and-stream-input-format).

3. Creating the streaming input
  * The `additions.snap` file is used to create edge additions stream file, `additions.stream` as follows:
    ```bash
    sed -e 's/^/a\t/' additions.snap > additions.stream
    ```
  * The `initial_graph.snap` file is used to create edge deletions stream file, `deletions.stream`. This is one of the steps to ensure that all the deletion operations actually result in edge deletions.
    ```bash
    sed -e 's/^/d\t/' initial_graph.snap > deletions.stream
    ```
  * Then to combine the additions and deletions file into the same file, the following command can be used:
    ```bash
    paste -d "\n" additions.stream deletions.stream > update.stream
    ```
  * Shuffle stream file:
    ```bash
      shuf update.stream -o update_shuffle.stream
    ```
4. COEM applications work on bipartite input graphs. In order to mock bipartite graphs using non-bipartite graphs, we create a partitions file which contains all the vertices belonging to a single partition. We do this by randomly assigning each vertex in the graph to one of the partitions.
    ```bash
    cd tools/generators
    # For an input graph with 8298 vertices, do the following.
    python CreatePartitions.py --n=8298 --outputFile="Partition1_wiki"
    ```
5. Label Propagation and COEM use an initial set of seed vertices as the frontier. We randomly select 5% of the vertices as the seeds for each input graph. The seeds file can be generated as follows:
    ```bash
    cd tools/generators
    # For an input graph with 8298 vertices, do the following.
    python CreateSeeds.py --n=8298 --outputFile="LabelSeeds_wiki" --seedsPercent=0.05
    ```


### 3. Running Applications 

The command for running each application with the `wiki_vote` graph is provided below:
```bash
# Running PageRank
./PageRank -maxIters 10 -fixedBatchSize -enforceEdgeValidity -nEdges 1000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj

# Running Label Propagation
./LabelPropagation -maxIters 10 -fixedBatchSize -enforceEdgeValidity -seedsFile ../inputs/wiki_vote/LabelSeeds_wiki -features 2 -nEdges 1000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj

# Running COEM
./COEM -s -maxIters 10 -fixedBatchSize -enforceEdgeValidity -seedsFile ../inputs/wiki_vote/LabelSeeds_wiki -partitionsFile ../inputs/wiki_vote/Partition1_wiki -nEdges 1000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj.un
```

**IMPORTANT: Ensure that `-fixedBatchSize` flag is present.** This ensures that the specified number of edge operations `-nEdges` are correctly performed. In addition to this, `-enforceEdgeValidity` flag validates all edge operations before counting them towards the batch size requirement (for example, deleting edges not present in the graph is considered as an invalid edge operation).

##  4. Experiments in Paper

### 4.1 Evaluating the System

#### 4.1.1 Understanding the Output
Below is a sample output of running PageRank with `wiki_vote` graph: 
```bash
Graph created
Initializing engine ....
Number of batches: 2
Creating dependency structure ....
Initializing dependency structure ....
Finished initializing engine

   Iteration,        Time
average cache miss: 0.000000
average delta time: 0.000000
Initial graph processing : 0.006726
Number of iterations : 10

Opening Stream: Waiting for writer to open...
Stream opened
Current_batch: 1
Batch Size: 1000
Reading Time : 0.013265
Edge deletion time : 0.001108
Edge addition time : 0.000289
Edge Additions in batch: 500
Edge Deletions in batch: 500
vertex compute time: 0.001193
edge computing time: 0.004532
other time: 0.002526
Finished batch : 0.008295
```
The `Finished batch` line gives the execution time of SVSIG in seconds to process each epoch of edge updates. In the above output sample, SVSIG took 0.008295 seconds. The execution times of `SVSIG` throughout the paper (including TABLE II, Fig. 7 and Fig. 8) are obtained from the `Finished batch` line.

The time taken in seconds to apply the Edge computation, Vertex computation and Other is printed as `edge computing time`, `vertex compute time` and `other time`. In the above output sample, it took 0.004532 seconds to apply Edge computation, 0.001193 seconds to apply Vertex computation and 0.002526 seconds to apply Other.

#### 4.1.2 Varying Batch Sizes 

Throughout our evaluation, we presented results for different applications running with different input batch sizes. This is achieved by varying the `-numberOfUpdateBatches` command-line parameter.

#### 4.1.3 Number of Edges Activated

Fig. 4, Fig. 5 and Fig. 7 show the amount of work done in terms of the number of edges activated by SVSIG. This is obtained by performing the following steps.

Compile the application with `COUNTINGEDGE=1` as follows:
```bash
make clean 
export LD_PRELOAD=../lib/mimalloc/out/release/libmimalloc.so
make MULTIAGGREGATIONDEVIATION=1 COUNTINGEDGE=1 PageRank
```

And then, run the application:
```bash
./PageRank -maxIters 10 -nEdges 10000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 10 ../inputs/wiki_vote/wiki_vote_initial.adj
```

The output will include the following content:

```bash
1       40477
2       45455
3       28807
4       28639
5       28407
6       28049
7       27343
8       25481
9       22713
10      19566
all activate edges: 294937
```


In this output, numbers 1 to 10 represent the iteration times, and the numbers on the right side indicate the edges activated in this iteration. For example, `1  40477` means that 40477 edges were activated in the first iteration. "All activated edges: 294937" represents the total number of edges activated in this epoch.
