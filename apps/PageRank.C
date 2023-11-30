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

#include "../core/common/utils.h"
#include "../core/graphBolt/GraphBoltEngine_simple.h"
#include "../core/main.h"
#include <math.h>
#include "../core/graph/graph.h"

// ======================================================================
// PAGERANKINFO
// ======================================================================
template <class vertex>
class PageRankInfo {
public:
  // Should I just use vectors for this?
  graph<vertex> *my_graph;
  uintV n;
  double epsilon;
  double deviation;
  double damping;
  long *out_degrees;

  PageRankInfo() : my_graph(nullptr), n(0), epsilon(0), damping(0), out_degrees(nullptr) {}

  PageRankInfo(graph<vertex> *_my_graph, uintV _n, double _epsilon, double _damping, double _deviation)
      : my_graph(_my_graph), n(_n), epsilon(_epsilon), damping(_damping), deviation(_deviation) {
    if (n > 0) {
      out_degrees = newA(long, n);
      parallel_for(uintV i = 0; i < n; i++) { out_degrees[i] = 0; }
    }
  }

  void init(){
    parallel_for(uintV i = 0; i < n; i++) {
      out_degrees[i] = my_graph->V[i].getOutDegree();
    }    
  }

  void copy(const PageRankInfo &object) {
    if (object.n > n) {
      if (n == 0) {
        n = object.n;
        out_degrees = newA(long, n);
      } else {
        // realloc
        n = object.n;
        out_degrees = renewA(long, out_degrees, n);
      }
    }
    long min_n = std::min(object.n, n);
    parallel_for(uintV i = 0; i < min_n; i++) {
      out_degrees[i] = object.out_degrees[i];
    }
    epsilon = object.epsilon;
    damping = object.damping;
  }

  void reProcess(const PageRankInfo &object) {

  }

  ~PageRankInfo() {
    if (n > 0)
      deleteA(out_degrees);
  }

  void processUpdates(edgeArray &edge_additions, edgeArray &edge_deletions) {
    // Increase out_degrees array size
    if (edge_additions.maxVertex >= n) {
      uintV n_old = n;
      n = edge_additions.maxVertex + 1;
      out_degrees = renewA(long, out_degrees, n);
      parallel_for(uintV i = n_old; i < n; i++) { out_degrees[i] = 0; }
    }

    parallel_for(long i = 0; i < edge_additions.size; i++) {
      uintV source = edge_additions.E[i].source;
      uintV destination = edge_additions.E[i].destination;
      writeAdd(&out_degrees[source], (long)1);
    }
    parallel_for(long i = 0; i < edge_deletions.size; i++) {
      uintV source = edge_deletions.E[i].source;
      uintV destination = edge_deletions.E[i].destination;
      writeAdd(&out_degrees[source], (long)-1);
    }
  }

  void cleanup() {}

  inline long getOutDegree(const uintV &v) {
    return (v < n) ? out_degrees[v] : 0;
  }
};

// ======================================================================
// AGGREGATEVALUE AND VERTEXVALUE INITIALIZATION
// ======================================================================
double initial_source_deviation_value = 0;
double initial_aggregation_value = 0;
double initial_vertex_value = 1;
// double initial_vertex_value = 0.15;
double aggregation_value_identity = 0;
double source_deviation_value_identity = 0;
double vertex_value_identity = 0;
template <class AggregationValueType, class GlobalInfoType>
inline void
initializeAggregationValue(const uintV &v,
                           AggregationValueType &v_aggregation_value,
                           const GlobalInfoType &global_info) {
  v_aggregation_value = initial_aggregation_value;
}

template <class VertexValueType, class GlobalInfoType>
inline void initializeVertexValue(const uintV &v,
                                  VertexValueType &v_vertex_value,
                                  const GlobalInfoType &global_info) {
  v_vertex_value = initial_vertex_value;
}

template <class AggregationValueType>
inline void initializeSourceDeviationValue(AggregationValueType &v_source_deviation) {
  v_source_deviation = initial_source_deviation_value;
}

template <class AggregationValueType>
inline AggregationValueType &aggregationValueIdentity() {
  return aggregation_value_identity;
}

template <class AggregationValueType>
inline AggregationValueType &sourceDeviationValueIdentity() {
  return source_deviation_value_identity;
}

template <class VertexValueType> inline VertexValueType &vertexValueIdentity() {
  return vertex_value_identity;
}
// ======================================================================
// ACTIVATE VERTEX/COMPUTE VERTEX FOR A GIVEN ITERATION
// ======================================================================
template <class GlobalInfoType>
inline bool forceActivateVertexForIteration(const uintV &v, int iter,
                                            const GlobalInfoType &global_info) {
  if (iter == 1) {
    return true;
  } else {
    return false;
  }
}

template <class GlobalInfoType>
inline bool forceComputeVertexForIteration(const uintV &v, int iter,
                                           const GlobalInfoType &global_info) {
  if (iter == 1) {
    return true;
  } else {
    return false;
  }
}
inline bool shouldUseDelta(int iter) {
  if (iter == 1) {
    return false;
  } else {
    return true;
  }
}

// ======================================================================
// ADD TO OR REMOVE FROM AGGREGATION VALUES
// ======================================================================
template <class AggregationValueType, class GlobalInfoType>
inline void addToAggregation(const AggregationValueType &incoming_value,
                             AggregationValueType &aggregate_value,
                             GlobalInfoType &global_info) {
  aggregate_value += incoming_value;
}

template <class AggregationValueType, class GlobalInfoType>
inline void addToAggregationAtomic(const AggregationValueType &incoming_value,
                                   AggregationValueType &aggregate_value,
                                   GlobalInfoType &global_info) {
  writeAdd(&aggregate_value, incoming_value);
}
template <class AggregationValueType>
inline void addToDeviation(const AggregationValueType &incoming_value,
                            AggregationValueType &v_deviation) {
  v_deviation += incoming_value;
}

template <class AggregationValueType>
inline void removeFromDeviation(const AggregationValueType &incoming_value,
                            AggregationValueType &v_deviation) {
  v_deviation -= incoming_value;
}

template <class AggregationValueType>
inline bool DeviationNotEqualZero(const AggregationValueType &v_deviation) {
  return fabs(v_deviation) > 1e-15;
}

template <class AggregationValueType, class GlobalInfoType>
inline void removeFromAggregation(const AggregationValueType &incoming_value,
                                  AggregationValueType &aggregate_value,
                                  GlobalInfoType &global_info) {
  aggregate_value -= incoming_value;
}

template <class AggregationValueType, class GlobalInfoType>
inline void
removeFromAggregationAtomic(const AggregationValueType &incoming_value,
                            AggregationValueType &aggregate_value,
                            GlobalInfoType &global_info) {
  writeAdd(&aggregate_value, -incoming_value);
}

// ======================================================================
// VERTEX COMPUTE FUNCTION AND DETERMINE END OF COMPUTATION
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline void computeFunction(const uintV &v,
                            const AggregationValueType &aggregation_value,
                            const VertexValueType &vertex_value_curr,
                            VertexValueType &vertex_value_next,
                            GlobalInfoType &global_info) {
  vertex_value_next =
      (1 - global_info.damping) + (global_info.damping * aggregation_value);
  // vertex_value_next = aggregation_value;
}

template <class VertexValueType, class GlobalInfoType>
inline bool notDelZero(const VertexValueType &value_curr,
                      const VertexValueType &value_next,
                      GlobalInfoType &global_info) {
  return (fabs(value_next - value_curr) > global_info.epsilon);
}

template <class AggregationValueType, class GlobalInfoType>
inline bool
sourceNotDelZero(const uintV &v,
                  const AggregationValueType &default_source,
                  const AggregationValueType &source_change,
                  GlobalInfoType &global_info) {
  if (global_info.getOutDegree(v) == 0) return false;
  // return source_change != default_source;
  // return fabs(default_source - source_change) > global_info.epsilon*0.0001;
  return fabs((default_source - source_change)*global_info.getOutDegree(v)) > global_info.epsilon*global_info.deviation;
  return fabs((default_source - source_change)*global_info.getOutDegree(v)) > global_info.epsilon*0.0001;
  return (fabs((default_source - source_change)*global_info.getOutDegree(v)) > global_info.epsilon);
}

// ======================================================================
// EDGE FUNCTIONS
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class StaticInfoType>
inline void sourceChangeInContribution(
    const uintV &v, AggregationValueType &v_change_in_contribution,
    const VertexValueType &v_value_prev, const VertexValueType &v_value_curr,
    StaticInfoType &global_info) {
  v_change_in_contribution =
      (global_info.getOutDegree(v) != 0)
          ? (v_value_curr - v_value_prev) / global_info.getOutDegree(v)
          : 0;
}

template <class AggregationValueType>
inline void removeAndAdd(const AggregationValueType &remove_value,
                          const AggregationValueType &add_value,
                          AggregationValueType &v_change_in_contribution) {
  v_change_in_contribution = add_value - remove_value;
}

template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline void sourceChangeInContribution(
    const uintV &v, AggregationValueType &v_change_in_contribution,
    const VertexValueType &v_value_prev_old, const VertexValueType &v_value_curr_old,
    const VertexValueType &v_value_prev_new, const VertexValueType &v_value_curr_new,
    GlobalInfoType &static_info_old, GlobalInfoType &static_info_new,
    bool &new_change, bool &old_change) {
  AggregationValueType value_old_change = 
    (static_info_old.getOutDegree(v) != 0)
        ? (v_value_curr_old - v_value_prev_old) / static_info_old.getOutDegree(v)
        : 0;
  AggregationValueType value_new_change = 
    (static_info_new.getOutDegree(v) != 0)
        ? (v_value_curr_new - v_value_prev_new) / static_info_new.getOutDegree(v)
        : 0;
  removeAndAdd(value_old_change, value_new_change, v_change_in_contribution);
}

template <class AggregationValueType, class VertexValueType, class EdgeDataType,
          class GlobalInfoType>
inline bool edgeFunction(const uintV &u, const uintV &v,
                         const EdgeDataType &edge_weight,
                         const VertexValueType &u_value,
                         AggregationValueType &u_change_in_contribution,
                         GlobalInfoType &global_info) {
  return true;
}

// ======================================================================
// INCREMENTAL COMPUTING / DETERMINING FRONTIER
// =====================================================================
template <class GlobalInfoType>
inline void hasSourceChangedByUpdate(const uintV &v, UpdateType update_type,
                                     bool &activateInCurrentIteration,
                                     bool &forceComputeInCurrentIteration,
                                     GlobalInfoType &global_info,
                                     GlobalInfoType &global_info_old) {
  if (global_info.getOutDegree(v) != global_info_old.getOutDegree(v))
    activateInCurrentIteration = true;
}

template <class GlobalInfoType>
inline void hasDestinationChangedByUpdate(const uintV &v,
                                          UpdateType update_type,
                                          bool &activateInCurrentIteration,
                                          bool &forceComputeInCurrentIteration,
                                          GlobalInfoType &global_info,
                                          GlobalInfoType &global_info_old) {}

// ======================================================================
// HELPER FUNCTIONS
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
void printHistory(const uintV &v, AggregationValueType **agg_values,
                  VertexValueType **vertex_values, GlobalInfoType &info,
                  int history_iterations) {
  // for (int iter = 0; iter < history_iterations; iter++) {
  //   cout << iter << "," << agg_values[iter][v] << "," << vertex_values[iter][v]
  //        << "\n";
  // }
  for (int iter = 0; iter < history_iterations; iter++) {
    ofstream output_file;
    string curr_output_file_path = "../output/" + to_string(iter);
    output_file.open(curr_output_file_path, ios::out);
    output_file << fixed;
    output_file << setprecision(VAL_PRECISION2);
    for (uintV v = 0; v < info.n; v++) {
      output_file << agg_values[iter][v] << "\n";
      // output_file << agg_values[iter][info.my_graph->old_ids[v]] << "\n";
    }
    output_file.close();
  }
}

template <class GlobalInfoType>
void printAdditionalData(ofstream &output_file, const uintV &v,
                         GlobalInfoType &info) {}

// ======================================================================
// COMPUTE FUNCTION
// ======================================================================
template <class vertex> void compute(graph<vertex> &G, commandLine config) {
  uintV n = G.n;
  // int max_iters = config.getOptionLongValue("-maxIters", 10);
  int max_iters = config.getOptionLongValue("-maxIters", 10);
  double epsilon = config.getOptionDoubleValue("-epsilon", 0.01d);
  double deviation = config.getOptionDoubleValue("-deviation", 0.0001d);
  max_iters += 1;
  double damping = 0.85;

  PageRankInfo<vertex> global_info(&G, n, epsilon, damping, deviation);

  cout << "Initializing engine ....\n";
  GraphBoltEngineSimple<vertex, double, double, PageRankInfo<vertex>> engine(
      G, max_iters, global_info, false, config);
  engine.init();
  cout << "Finished initializing engine\n";

  engine.run();
}
