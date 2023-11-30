// Modifications Copyright (c) 2023 Zhixiong Li.
//
// Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef GRAPHBOLT_ENGINE_H
#define GRAPHBOLT_ENGINE_H

#include "../common/utils.h"
#include "AdaptiveExecutor.h"
#include "ingestor.h"
#include <vector>

enum UpdateType { edge_addition_enum, edge_deletion_enum };

#ifdef EDGEDATA
#else
struct EmptyEdgeData {};
typedef EmptyEdgeData EdgeData;

EdgeData emptyEdgeData;
#endif

// ======================================================================
// AGGREGATEVALUE AND VERTEXVALUE INITIALIZATION
// ======================================================================
// Set the initial aggregation value of the vertex. Usually, 0 for aggregations
// like sum and 1 for product.
template <class AggregationValueType, class GlobalInfoType>
inline void
initializeAggregationValue(const uintV &v,
                           AggregationValueType &v_aggregation_value,
                           const GlobalInfoType &global_info);

// Set the initial value of the vertex. For example, initial value of vertex is
// set as 1 for PageRank
template <class VertexValueType, class GlobalInfoType>
inline void initializeVertexValue(const uintV &v,
                                  VertexValueType &v_vertex_value,
                                  const GlobalInfoType &global_info);
template <class AggregationValueType>
inline void initializeSourceDeviationValue(AggregationValueType &v_source_deviation);

// Return the identity value for the aggregation value. For sum, the identity
// value is 0. For product, it is 1.
template <class AggregationValueType>
// Return the identity value for the vertex value.
inline AggregationValueType &aggregationValueIdentity();
template <class VertexValueType> inline VertexValueType &vertexValueIdentity();
template <class AggregationValueType>
inline AggregationValueType &sourceDeviationValueIdentity();

// ======================================================================
// ACTIVATE VERTEX/COMPUTE VERTEX FOR A GIVEN ITERATION
// ======================================================================
// Used to determine whether a vertex should be forced to be active in a given
// iteration. For example in PageRank, all the vertices should be active for the
// first iteration. In CF (ALS), in first iteration, all partition2 vertices
// should be active and in second iteration, all partition1 vertices should be
// active
template <class GlobalInfoType>
inline bool forceActivateVertexForIteration(const uintV &v, int iter,
                                            const GlobalInfoType &global_info);

// By default, all the vertices which receive some change from its inNeighbor
// will computes its value. For some algorithms, we need to compute the value of
// the vertex irrespective of whether the vertex receives any change for that
// iteration. This function is used to force the computation of a vertex at a
// given iteration.
template <class GlobalInfoType>
inline bool forceComputeVertexForIteration(const uintV &v, int iter,
                                           const GlobalInfoType &global_info);

// Usually, for the first iteration, the vertices propagate the entire value to
// its outNeighbor. For example, this function should return false for iter=1
// and return true for other iterations. In the case of CF, we have 2
// partitions. So, the function should return true for iter=1 or iter=2. Use
// this function to customize when delta based incremental computation should be
// performed.
inline bool shouldUseDelta(int iter);

// ======================================================================
// ADD TO OR REMOVE FROM AGGREGATION VALUES
// ======================================================================
// Function to add a value to an aggregate_value. For aggregations like sum,
// simply add the incoming_value to the aggregate_value. No need for locks or
// CAS within this function.
template <class AggregationValueType, class GlobalInfoType>
inline void addToAggregation(const AggregationValueType &incoming_value,
                             AggregationValueType &aggregate_value,
                             GlobalInfoType &global_info);

// Function to add a value to an aggregate_value. For aggregations like sum,
// simply add the incoming_value to the aggregate_value. Use CAS for atomically
// updating the aggregate_value.
template <class AggregationValueType, class GlobalInfoType>
inline void addToAggregationAtomic(const AggregationValueType &incoming_value,
                                   AggregationValueType &aggregate_value,
                                   GlobalInfoType &global_info);

template <class AggregationValueType>
inline void addToDeviation(const AggregationValueType &incoming_value,
                            AggregationValueType &v_deviation);
template <class AggregationValueType>
inline void removeFromDeviation(const AggregationValueType &incoming_value,
                            AggregationValueType &v_deviation);

template <class AggregationValueType>
inline bool DeviationNotEqualZero(const AggregationValueType &v_deviation);

// Function to remove a value from an aggregate_value. For aggregations like
// sum, simply subtract the incoming_value from the aggregate_value. No need for
// locks or CAS within this function.
template <class AggregationValueType, class GlobalInfoType>
inline void removeFromAggregation(const AggregationValueType &incoming_value,
                                  AggregationValueType &aggregate_value,
                                  GlobalInfoType &global_info);

// Function to remove a value from an aggregate_value. For aggregations like
// sum, simply subtract the incoming_value from the aggregate_value. Use CAS for
// atomically updating the aggregate_value.
template <class AggregationValueType, class GlobalInfoType>
inline void
removeFromAggregationAtomic(const AggregationValueType &incoming_value,
                            AggregationValueType &aggregate_value,
                            GlobalInfoType &global_info);

// ======================================================================
// VERTEX COMPUTE FUNCTION AND DETERMINE END OF COMPUTATION
// ======================================================================
// Function to compute the value of a vertex 'v', provided its current
// aggregation_value and previous vertex_value. For PageRank, we can directly
// compute vertex_value_curr from its aggregation_value. For Label Propagation,
// the computation is based on the aggregation_value as well as its previous
// vertex_value_prev.
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline void computeFunction(const uintV &v,
                            const AggregationValueType &aggregation_value,
                            const VertexValueType &vertex_value_prev,
                            VertexValueType &vertex_value_curr,
                            GlobalInfoType &global_info);

// Function which determines whether a vertex value has converged. For PageRank,
// if the |vertex_value_curr - vertex_value_prev| < 0.01, the vertex value does
// not have significant change to push to its outNeighbors in the next
// iteration. In this case, return false. Otherwise, return true. Note: Even
// when the vertex does not have significant change in its value at iteration
// 'i', the change could become significant in a future iteration as it keeps
// accumulating these changes at every iteration.
template <class VertexValueType, class GlobalInfoType>
inline bool notDelZero(const VertexValueType &vertex_value_prev,
                      const VertexValueType &vertex_value_curr,
                      GlobalInfoType &global_info);

// ======================================================================
// EDGE FUNCTIONS
// ======================================================================
// For an edge (u,v), the edge function is performed in 3 steps: (1) The change
// in contribution for the source vertex, which handles the computation specific
// to that source vertex, (2) The edgeFunction which incorporates the
// computations specific to that edge, and (3) The final change_in_contribution
// for the edge (u, v) is aggregated into the aggregate_value of v using
// addToAggregationAtomic.

// For a vertex u, given its value in current iteraion and previous iteraion,
// compute the change_in_contribution. For a vertex v in PageRank, the same
// value (PR[v]/out_degree[v]) has to be pushed to all its outNeighbors. In the
// case of delta-based computation, it pushes out (PR_curr[v] -
// PR_prev[v])/out_degree[v]. NOTE: change_in_contribution can be directly
// updated in place. No LOCKS/CAS required for updating u_change_in_contribution
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline void sourceChangeInContribution(
    const uintV &v, AggregationValueType &v_change_in_contribution,
    const VertexValueType &v_value_prev, const VertexValueType &v_value_curr,
    GlobalInfoType &static_info);

template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline void sourceChangeInContribution(
    const uintV &v, AggregationValueType &v_change_in_contribution,
    const VertexValueType &v_value_prev_old, const VertexValueType &v_value_curr_old,
    const VertexValueType &v_value_prev_new, const VertexValueType &v_value_curr_new,
    GlobalInfoType &static_info_old, GlobalInfoType &static_info_new,
    bool &new_change, bool &old_change);
// For a vertex u, given its value in current iteration and its computed
// change_in_contribution, update change_in_contribution for a given edge (u,
// v). For example, in Label Propagation, the edge_data of (u, v) is
// multiplied to the change_in_contribution of each factor . Return false if the
// value from u should not be included in the aggregation value of v. Return
// true otherwise. NOTE: The changes to change_in_contribution should be made in
// place. No LOCKS/CAS required for updating u_change_in_contribution
template <class AggregationValueType, class VertexValueType, class EdgeDataType,
          class GlobalInfoType>
inline bool edgeFunction(const uintV &u, const uintV &v,
                         const EdgeDataType &edge_data,
                         const VertexValueType &u_value,
                         AggregationValueType &u_change_in_contribution,
                         GlobalInfoType &global_info);

// ======================================================================
// INCREMENTAL COMPUTING / DETERMINING FRONTIER
// ======================================================================
// For a given edge addition (u, v), define how it affects the source vertex in
// the first iteration.
// activateInCurrentIteration defines if the source vertex will be
// active for the first iteration. For example, in PageRank, if the out_degree
// of u has changed, then the value of u pushed to all its outNeighbors is
// changed too. So, if out_degree(u) in the updated graph is not equal to its
// out_degree in the previous version of the graph, return true. Otherwise,
// return false. NOTE: Ensure that you store the relevant information required
// (for algorithm/application) for a given graph version in the global_info
// object. Example: out_degree of all vertices has to be stored for a given
// graph.
// forceComputeInCurrentIteration defines if the source vertex will have to be
// recomputed in the first iteration. For example, in COEM, if the sum of
// inWeights of a vertex changes, then computeFuntion() should be
// called for that vertex in the first iteration.
template <class GlobalInfoType>
inline void hasSourceChangedByUpdate(const uintV &v, UpdateType update_type,
                                     bool &activateInCurrentIteration,
                                     bool &forceComputeInCurrentIteration,
                                     GlobalInfoType &global_info,
                                     GlobalInfoType &global_info_old);

// For a given edge addition (u, v), define if the destination vertex has
// changed.
// activateInCurrentIteration defines if the destination vertex will be active
// for the first iteration. NOTE: Ensure that you store the relevant information
// required (for algorithm/application) for a given graph version in the
// global_info object. Example: out_degree of all vertices has to be stored for
// a given graph.
// forceComputeInCurrentIteration defines if the destination vertex will have to
// be recomputed in the first iteration.
template <class GlobalInfoType>
inline void hasDestinationChangedByUpdate(const uintV &v,
                                          UpdateType update_type,
                                          bool &activateInCurrentIteration,
                                          bool &forceComputeInCurrentIteration,
                                          GlobalInfoType &global_info,
                                          GlobalInfoType &global_info_old);

// ======================================================================
// HELPER FUNCTIONS
// ======================================================================
// Helper function for printing additional data while printing to output file
template <class GlobalInfoType>
void printAdditionalData(ofstream &output_file, const uintV &v,
                         GlobalInfoType &info);

// Helper function for printing the dependency data - Useful for debugging
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
void printHistory(const uintV &v, AggregationValueType **agg_values,
                  VertexValueType **vertex_values, GlobalInfoType &info,
                  int history_iterations);

// ======================================================================
// GRAPHBOLT ENGINE
// ======================================================================
template <class vertex, class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
class GraphBoltEngine {

public:
  graph<vertex> &my_graph;
  commandLine config;

  // TODO: Currently, history_iterations = max_iterations
  int max_iterations;
  int history_iterations;
  int converged_iteration;
  bool use_lock;

  double average_cache_miss;
  double average_delta_time;
  double average_delta_compute_time;

  RWLock *vertex_locks;

  // Dependency information
  AggregationValueType **aggregation_values;
  VertexValueType **vertex_values;

  // Current graph
  long n;
  GlobalInfoType &global_info;
  AggregationValueType *delta;
  bool use_source_contribution;
  AggregationValueType *source_change_in_contribution;
  #ifdef MULTI
  AggregationValueType **source_deviation;
  AggregationValueType *temp_source_deviation;
  AggregationValueType *del_source_deviation;
  bool *is_modify_source;
  #endif
  #ifdef ONE
   AggregationValueType *temp_source_deviation;
  #endif
  bool first_delta_compute;

  // Previous graph
  long n_old;
  GlobalInfoType global_info_old;

  // temporary structures
  VertexValueType *vertex_value_old_next;
  VertexValueType *vertex_value_old_curr;
  VertexValueType *vertex_value_old_prev;

  // TODO : Replace with more bitmaps
  bool *all;
  bool *frontier_curr;
  bool *frontier_next;
  bool *changed;
  #ifdef TINY
  bool *first_changed;
  #endif
  bool *retract;
  bool *propagate;

  // Stream Ingestor
  Ingestor<vertex> ingestor;
  int current_batch;

  // For adaptive switching
  AdaptiveExecutor adaptive_executor;
  bool ae_enabled;

  long total_work_done;

  unsigned long activate_edge;
  unsigned long all_activate_edge;
  double edge_map_time;
  double vertex_compute_time;
  double other_time;
  double all_edge_map_time;
  double all_vertex_compute_time;
  double all_other_time;

  // ======================================================================
  // CONSTRUCTOR / INIT
  // ======================================================================
  GraphBoltEngine(graph<vertex> &_my_graph, int _max_iter,
                  GlobalInfoType &_global_info, bool _use_lock,
                  commandLine _config)
      : my_graph(_my_graph), max_iterations(_max_iter),
        history_iterations(_max_iter), converged_iteration(0),
        global_info(_global_info), use_lock(_use_lock), global_info_old(),
        config(_config), ingestor(_my_graph, _config), current_batch(0),
        adaptive_executor(history_iterations), total_work_done(0) {
    n = my_graph.n;
    n_old = 0;
    first_delta_compute = 1;
    all_edge_map_time = 0; all_other_time = 0; all_vertex_compute_time = 0;
    if (use_lock) {
      cout << "Using locks for edge operations\n";
      createLocks();
      initLocks();
    }
    ae_enabled = config.getOptionValue("-ae");
  }

  void init() {
    cout << "Creating dependency structure ....\n";
    createDependencyData();
    createTemporaryStructures();
    createVertexSubsets();
    cout << "Initializing dependency structure ....\n";
    initVertexSubsets();
    initTemporaryStructures();
    initDependencyData();
    activate_edge = 0;
    all_activate_edge = 0;
  }

  ~GraphBoltEngine() {
    freeDependencyData();
    freeVertexSubsets();
    if (use_lock) {
      destroyLocks(n);
      deleteA(vertex_locks);
    }
    global_info.cleanup();
  }

  // ======================================================================
  // CREATE / DESTROY LOCKS
  // ======================================================================
  void createLocks() { vertex_locks = newA(RWLock, n); }
  void resizeLocks() {
    vertex_locks = renewA(RWLock, vertex_locks, n);
    initLocks(n_old, n);
  }
  void initLocks() { initLocks(0, n); }
  void initLocks(long start_index, long end_index) {
    parallel_for(long i = start_index; i < end_index; i++) {
      vertex_locks[i].init();
    }
  }
  void destroyLocks(long array_size) {
    parallel_for(long i = 0; i < array_size; i++) { vertex_locks[i].destroy(); }
  }

  // ======================================================================
  // DEPENDENCY DATA STORAGE
  // ======================================================================
  void createDependencyData() {
    aggregation_values = newA(AggregationValueType *, history_iterations);
    vertex_values = newA(VertexValueType *, history_iterations);
#ifdef MULTI
    // source_deviation = newA(AggregationValueType *, history_iterations+1);
    source_deviation = newA(AggregationValueType *, history_iterations);
#endif
    for (int i = 0; i < history_iterations; i++) {
      aggregation_values[i] = newA(AggregationValueType, n);
      vertex_values[i] = newA(VertexValueType, n);
#ifdef MULTI
      source_deviation[i] = newA(AggregationValueType, n);
#endif
    }
  }
  void resizeDependencyData() {
    for (int i = 0; i < history_iterations; i++) {
      aggregation_values[i] =
          renewA(AggregationValueType, aggregation_values[i], n);
      vertex_values[i] = renewA(VertexValueType, vertex_values[i], n);
#ifdef MULTI
      source_deviation[i] = renewA(AggregationValueType, source_deviation[i], n);
#endif
    }
    initDependencyData(n_old, n);
  }
  void freeDependencyData() {
    for (int i = 0; i < history_iterations; i++) {
      deleteA(aggregation_values[i]);
      deleteA(vertex_values[i]);
#ifdef MULTI
      deleteA(source_deviation[i]);
#endif
    }
// #ifdef MULTI
//     deleteA(source_deviation[history_iterations]);
// #endif
    deleteA(aggregation_values);
    deleteA(vertex_values);
#ifdef MULTI
    deleteA(source_deviation);
#endif
  }
  void initDependencyData() { initDependencyData(0, n); }
  void initDependencyData(long start_index, long end_index) {
    for (int iter = 0; iter < history_iterations; iter++) {
      parallel_for(long v = start_index; v < end_index; v++) {
        initializeAggregationValue<AggregationValueType, GlobalInfoType>(
            v, aggregation_values[iter][v], global_info);
        initializeVertexValue<VertexValueType, GlobalInfoType>(
            v, vertex_values[iter][v], global_info);
#ifdef MULTI
        initializeSourceDeviationValue<AggregationValueType>(
            source_deviation[iter][v]);
#endif
      }
    }
// #ifdef MULTI
//     parallel_for(long v = start_index; v < end_index; v++) {
//       initializeSourceDeviationValue<AggregationValueType>(
//           source_deviation[history_iterations][v]);
//     }
// #endif
  }

  // ======================================================================
  // TEMPORARY STRUCTURES USED BY THE BSP ENGINE
  // ======================================================================
  virtual void createTemporaryStructures() {
    vertex_value_old_next = newA(VertexValueType, n);
    vertex_value_old_curr = newA(VertexValueType, n);
    vertex_value_old_prev = newA(VertexValueType, n);
    delta = newA(AggregationValueType, n);
#ifdef MULTI
    temp_source_deviation = newA(AggregationValueType, n);
    del_source_deviation = newA(AggregationValueType, n);
    is_modify_source = newA(bool, n);
#endif 
#ifdef ONE
    temp_source_deviation = newA(AggregationValueType, n);
#endif
    if (use_source_contribution)
      source_change_in_contribution = newA(AggregationValueType, n);
  }
  virtual void resizeTemporaryStructures() {
    vertex_value_old_next = renewA(VertexValueType, vertex_value_old_next, n);
    vertex_value_old_curr = renewA(VertexValueType, vertex_value_old_curr, n);
    vertex_value_old_prev = renewA(VertexValueType, vertex_value_old_prev, n);
    delta = renewA(AggregationValueType, delta, n);
#ifdef MULTI
    temp_source_deviation = renewA(AggregationValueType, temp_source_deviation, n);
    del_source_deviation = renewA(AggregationValueType, del_source_deviation, n);
    is_modify_source = renewA(bool, is_modify_source, n);
#endif
#ifdef ONE
    temp_source_deviation = renewA(AggregationValueType, temp_source_deviation, n);
#endif
    if (use_source_contribution)
      source_change_in_contribution =
          renewA(AggregationValueType, source_change_in_contribution, n);
  }
  virtual void freeTemporaryStructures() {
    deleteA(vertex_value_old_next);
    deleteA(vertex_value_old_curr);
    deleteA(vertex_value_old_prev);
    deleteA(delta);
#ifdef MULTI
    deleteA(temp_source_deviation);
    deleteA(del_source_deviation);
    deleteA(is_modify_source);
#endif
#ifdef ONE
    deleteA(temp_source_deviation);
#endif
    if (use_source_contribution)
      deleteA(source_change_in_contribution);
  }
  virtual void initTemporaryStructures() { initTemporaryStructures(0, n); }
  virtual void initTemporaryStructures(long start_index, long end_index) {
    parallel_for(long v = start_index; v < end_index; v++) {
      vertex_value_old_next[v] = vertexValueIdentity<VertexValueType>();
      vertex_value_old_curr[v] = vertexValueIdentity<VertexValueType>();
      vertex_value_old_prev[v] = vertexValueIdentity<VertexValueType>();
      delta[v] = aggregationValueIdentity<AggregationValueType>();
#ifdef MULTI
      temp_source_deviation[v] = sourceDeviationValueIdentity<AggregationValueType>();
      del_source_deviation[v] = sourceDeviationValueIdentity<AggregationValueType>();
#endif
#ifdef ONE 
      temp_source_deviation[v] = sourceDeviationValueIdentity<AggregationValueType>();
#endif
      if (use_source_contribution)
        source_change_in_contribution[v] =
            aggregationValueIdentity<AggregationValueType>();
    }
  }

  // ======================================================================
  // VERTEX SUBSETS USED BY THE BSP ENGINE
  // ======================================================================
  void createVertexSubsets() {
    all = newA(bool, n);
    frontier_curr = newA(bool, n);
    frontier_next = newA(bool, n);
    changed = newA(bool, n);
#ifdef TINY
    first_changed = newA(bool, n);
#endif
    retract = newA(bool, n);
    propagate = newA(bool, n);
  }
  void resizeVertexSubsets() {
    all = renewA(bool, all, n);
    frontier_curr = renewA(bool, frontier_curr, n);
    frontier_next = renewA(bool, frontier_next, n);
    changed = renewA(bool, changed, n);
#ifdef TINY
    first_changed = renewA(bool, first_changed, n);
#endif
    retract = renewA(bool, retract, n);
    propagate = renewA(bool, propagate, n);
    initVertexSubsets(n_old, n);
  }
  void freeVertexSubsets() {
    deleteA(all);
    deleteA(frontier_curr);
    deleteA(frontier_next);
    deleteA(changed);
#ifdef TINY
    deleteA(first_changed);
#endif
    deleteA(retract);
    deleteA(propagate);
  }
  void initVertexSubsets() { initVertexSubsets(0, n); }
  void initVertexSubsets(long start_index, long end_index) {
    parallel_for(long j = start_index; j < end_index; j++) {
      all[j] = 1;
      frontier_curr[j] = 0;
      frontier_next[j] = 0;
      changed[j] = 0;
#ifdef TINY
      first_changed[j] = 0;
#endif
      retract[j] = 0;
      propagate[j] = 0;
    }
  }

  // ======================================================================
  // PROCESS VERTEX ADDITION
  // ======================================================================
  void processVertexAddition(uintV maxVertex) {
    n_old = n;
    n = maxVertex + 1;
    resizeDependencyData();
    resizeTemporaryStructures();
    resizeVertexSubsets();
    if (use_lock) {
      cout << "Resizing locks\n";
      resizeLocks();
    }
  }

  // ======================================================================
  // HELPER FUNCTIONS TO PRINT OUTPUT
  // ======================================================================
  void testPrint() {
    // cout << setprecision(VAL_PRECISION);
    // for (auto curr : debug_vertices) {
    //   cout << "Vertex " << curr << "\n";
    //   cout << "Indegree " << my_graph.V[curr].getInDegree() << "\n";
    //   cout << "Outdegree " << my_graph.V[curr].getOutDegree() << "\n";
    //   printHistory(curr, aggregation_values, vertex_values, global_info,
    //                converged_iteration);
    // }
    printHistory(0, aggregation_values, vertex_values, global_info,
                   converged_iteration);
  }

  void printOutput0() {}
  void printOutput() {
    string output_file_path = config.getOptionValue("-outputFile", "/tmp/");
    bool should_print = true;
    if (output_file_path.compare("/tmp/") == 0) {
      should_print = false;
    }
    if (should_print) {
      string curr_output_file_path =
          output_file_path + to_string(current_batch);
      std::cout << "Printing to file : " << curr_output_file_path << "\n";
      ofstream output_file;
      output_file.open(curr_output_file_path, ios::out);
      output_file << fixed;
      output_file << setprecision(VAL_PRECISION2);
      
      for (uintV v = 0; v < n; v++) {
        //output_file << v << " " << my_graph.V[v].getInDegree() << " "
          //          << my_graph.V[v].getOutDegree() << " ";
        // printAdditionalData(output_file, v, global_info);
        output_file << vertex_values[converged_iteration][v] << "\n";
      }
      
    }
    cout << "\n";
  }
  // ======================================================================
  // RUN AND INITIAL COMPUTE
  // ======================================================================
  void run() {
    // ofstream file("neighbor.txt");
    // for(uintV v = 0; v < n; v++) {
    //   intE deg = my_graph.V[v].getOutDegree();
    //   for(intE e = 0; e < deg; e++) {
    //     file << my_graph.V[v].getOutNeighbor(e) << endl;
    //   }
    // }
    // file.close();
    initialCompute();

    // ======================================================================
    // Incremental Compute - Get the next update batch from ingestor
    // ======================================================================
    ingestor.validateAndOpenFifo();
    average_delta_compute_time = 0;
    while (ingestor.processNextBatch()) {
      current_batch++;
      // if (current_batch == ingestor.number_of_batches) {
      //   ofstream file("neighbor.txt");
      //   for(uintV v = 0; v < n; v++) {
      //     intE deg = my_graph.V[v].getOutDegree();
      //     if (deg > 1) {
      //       std::sort(my_graph.V[v].getOutNeighbors(), my_graph.V[v].getOutNeighbors()+deg);
      //     }
          
      //     for(intE e = 0; e < deg; e++) {
      //       file << my_graph.V[v].getOutNeighbor(e) << endl;
      //     }
      //   }
      //   file.close();
      // }
      total_work_done = 0;
      edgeArray &edge_additions = ingestor.getEdgeAdditions();
      edgeArray &edge_deletions = ingestor.getEdgeDeletions();
      // ingestor.edge_additions and ingestor.edge_deletions have been added
      // to the graph datastructure. Now, refine using it.
#ifdef MULTI
      if(ingestor.edge_count >= EDGE_MAX_COUNT) {
        if (edge_additions.maxVertex >= n) {
          processVertexAddition(edge_additions.maxVertex);
        }
        global_info.processUpdates(edge_additions, edge_deletions);
        // global_info_old.reProcess(global_info);
        // parallel_for(uintV u = 0; u < n; u++) {
        //   initializeAggregationValue<AggregationValueType, GlobalInfoType>(
        //       u, aggregation_value[0][u], global_info);
        //   initializeVertexValue<VertexValueType, GlobalInfoType>(
        //       u, vertex_values[0][u], global_info);
        // }
        // for(int i = 0; i < max_iterations; i++) {
        //   parallel_for(uintV u = 0; u < n; u++) {
        //     source_deviation[i][u] = sourceDeviationValueIdentity<AggregationValueType>();
        //   }
        // }
        parallel_for(uintV v = 0; v < n; v++) {
          frontier_curr[v] = forceActivateVertexForIteration(v, 1, global_info);
          frontier_next[v] = 0;
        }
        traditionalIncrementalComputation(1);
        ingestor.edge_count = 0;
      } else {
        deltaCompute(edge_additions, edge_deletions);
      }
#else
      deltaCompute(edge_additions, edge_deletions);
      // if(current_batch % 10 == 0) {
      //   printOutput();
      // }
#endif
#ifdef TENTIME
      if(current_batch %10 == 0 && current_batch >= 500) {
      // if(current_batch %10 == 0) {
        printOutput();
      }
#endif
    }
    cout << "average deltaCompute time: " << average_delta_compute_time / ingestor.number_of_batches << endl;
    cout << "average reading time: " << ingestor.reading_time / ingestor.number_of_batches << endl;
    cout << "edge average processing time: " << ingestor.edge_process_time / ingestor.number_of_batches << endl;
    cout << "average IO time: " << (ingestor.reading_time+ingestor.edge_process_time) / ingestor.number_of_batches << endl;
    cout << "average processing time: " << (average_delta_compute_time+ingestor.reading_time+ingestor.edge_process_time) / ingestor.number_of_batches << endl;
#ifdef SUBTIME
    cout << "average vertex computing time: " << all_vertex_compute_time / ingestor.number_of_batches << endl;
    cout << "average edge computing time: " << all_edge_map_time / ingestor.number_of_batches << endl;
    cout << "average other time: " << all_other_time / ingestor.number_of_batches << endl;
#endif
    // printOutput();
    
    // ofstream file("neighbor.txt");
    // for(uintV v = 0; v < n; v++) {
    //   intE deg = my_graph.V[v].getOutDegree();
    //   if (deg > 1) {
    //     std::sort(my_graph.V[v].getOutNeighbors(), my_graph.V[v].getOutNeighbors()+deg);
    //   }
      
    //   for(intE e = 0; e < deg; e++) {
    //     file << my_graph.V[v].getOutNeighbor(e) << endl;
    //   }
    // }
    // file.close();
    freeTemporaryStructures();
  }

  void initialCompute() {
    timer full_timer, t1;
    full_timer.start();
    t1.start();

    global_info.init();

    // Initilaize frontier
    // The other values are already initialized with the default values during
    // initialization of the GraphBoltEngine
    parallel_for(uintV v = 0; v < n; v++) {
      frontier_next[v] = 0;
      frontier_curr[v] = 0;
      frontier_curr[v] = forceActivateVertexForIteration(v, 1, global_info);
    }

    MY_TIMER_LOGS([&] {
      cout << "\n" << setw(PRINT_WIDTH + 1) << "Iteration,";
#ifdef EDGEWORK
      cout << setw(PRINT_WIDTH + 1) << "T_Edges," << setw(PRINT_WIDTH + 1)
           << "D_Edges,";
#endif
      cout << setw(PRINT_WIDTH + 1) << "Time\n";
    });

    int iters = traditionalIncrementalComputation(1);
    //
    // full_timer.next();

    // for (int i = 0; i < 10; i++) {
    //   initVertexSubsets();
    //   initTemporaryStructures();
    //   initDependencyData();
    //   parallel_for(uintV v = 0; v < n; v++) {
    //     frontier_next[v] = 0;
    //     frontier_curr[v] = 0;
    //     frontier_curr[v] = forceActivateVertexForIteration(v, 1, global_info);
    //   }
    //   traditionalIncrementalComputation(1);

    // }
    cout << "average cache miss: " << average_cache_miss << endl;
    cout << "average delta time: " << average_delta_time << endl;

    cout << "Initial graph processing : " << full_timer.stop() << "\n";
    cout << "Number of iterations : " << iters << "\n";
    printOutput();
    // testPrint();
  }

  void reset() {
  #ifdef MULTI
    // for(int i = 0; i < converged_iteration+1; i++) {
    //   parallel_for(uintV u = 0; u < n; u++) {
    //     source_deviation[u] = sourceDeviationValueIdentity<AggregationValueType>();
    //   }
    // }
    for(int i = 0; i < converged_iteration; i++) {
      parallel_for(uintV u = 0; u < n; u++) {
        source_deviation[u] = sourceDeviationValueIdentity<AggregationValueType>();
      }
    }
  #endif
  }

  virtual int traditionalIncrementalComputation(int start_iteration) = 0;
  virtual void deltaCompute(edgeArray &edge_additions,
                            edgeArray &edge_deletions) = 0;

  // ======================================================================
  // ADAPTIVE SWITCHING TO TRADITIONAL INCREMENTAL COMPUTATION
  // ======================================================================
  bool shouldSwitch(int iter, double dz_inc_iter_time) {
    if (iter > 0) {
      if (adaptive_executor.approximateTimeForCurrIter() < dz_inc_iter_time) {
        return true;
      }
    }
    // Update approximate_time_for_prev_iteration for next iteration
    parallel_for(uintV v = 0; v < n; v++) {
      if ((iter > 0 && notDelZero(vertex_values[iter][v],
                                 vertex_values[iter - 1][v], global_info)) ||
          forceActivateVertexForIteration(v, iter + 1, global_info)) {
        // We will use frontier_next as it is not being used by
        // deltaCompute
        frontier_next[v] = 1;
      } else {
        frontier_next[v] = 0;
      }
    }
    long active_edges =
        sequence::plusReduceDegree(my_graph.V, frontier_next, (long)n);
    adaptive_executor.updateApproximateTimeForEdges(active_edges);

    return false;
  }

  int performSwitch(int iter) {
    // If called at beginning of iteration, use iter-1 and iter-2 to decide
    // whether a vertex is active
    parallel_for(uintV v = 0; v < n; v++) {
      if (notDelZero(vertex_values[iter - 1][v], vertex_values[iter - 2][v],
                    global_info) ||
          forceActivateVertexForIteration(v, iter, global_info)) {
        // Update frontier_curr with active vertices for traditional
        // incremental processing
        frontier_curr[v] = 1;
        // Update source change in contribution
        sourceChangeInContribution<AggregationValueType, VertexValueType,
                                   GlobalInfoType>(
            v, source_change_in_contribution[v], vertex_values[iter - 2][v],
            vertex_values[iter - 1][v], global_info);

      } else {
        frontier_curr[v] = 0;
      }
      frontier_next[v] = 0;
    }
    cout << "iter: " << iter;
    cout << "*\n";
    // cout << "Switching to traditional incremental computation at iteration: "
    //      << iter << "\n";
    return traditionalIncrementalComputation(iter);
  }
};

#endif
