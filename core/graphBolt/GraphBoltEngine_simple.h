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

#ifndef GRAPHBOLT_ENGINE_SIMPLE_H
#define GRAPHBOLT_ENGINE_SIMPLE_H

#include "GraphBoltEngine.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>
#include <sys/ioctl.h>

// ======================================================================
// GRAPHBOLTENGINESIMPLE
// ======================================================================
template <class vertex, class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
class GraphBoltEngineSimple
    : public GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                             GlobalInfoType> {
public:
  GraphBoltEngineSimple(graph<vertex> &_my_graph, int _max_iter,
                        GlobalInfoType &_static_data, bool _use_lock,
                        commandLine _config)
      : GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>(_my_graph, _max_iter, _static_data,
                                        _use_lock, _config) {
    use_source_contribution = true;
  }

  // ======================================================================
  // TEMPORARY STRUCTURES USED BY THE SIMPLE ENGINE
  // ======================================================================
  void createTemporaryStructures() {
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::createTemporaryStructures();
  }
  void resizeTemporaryStructures() {
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::resizeTemporaryStructures();
    initTemporaryStructures(n_old, n);
  }
  void freeTemporaryStructures() {
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::freeTemporaryStructures();
  }
  void initTemporaryStructures() { initTemporaryStructures(0, n); }
  void initTemporaryStructures(long start_index, long end_index) {
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::initTemporaryStructures(start_index,
                                                             end_index);
  }

  static long perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                            int cpu, int group_fd, unsigned long flags) {
    return syscall(__NR_perf_event_open, hw_event, pid, cpu,
                   group_fd, flags);
  }
  // ======================================================================
  // TRADITIONAL INCREMENTAL COMPUTATION
  // ======================================================================
  // TODO : Currently, max_iterations = history_iterations.
  // Need to implement computation without history.
  int traditionalIncrementalComputation(int start_iteration) {
    timer iteration_timer, phase_timer;
    double misc_time, copy_time, phase_time, iteration_time;
    // timer my_time;
    double source_change_time, delta_time, all_delta_time = 0;
    unsigned long long all_count = 0;

    vertexSubset frontier_curr_vs(n, frontier_curr);
    bool use_delta = true;
    int iter = start_iteration;
    // unsigned long long activate_edge = 0;
    // unsigned long long all_activate_edge = 0;
    unsigned long long activate_vertex = 0;
    unsigned long long all_activate_vertex = 0;
    // ofstream file("inDegree.txt");
    // for(uintV v = 0; v < n; v++) {
    //   intE inDegree = my_graph.V[v].getInDegree();
    //   intE outDegree = my_graph.V[v].getOutDegree();
    //   file << inDegree << "\n";
    //   file << outDegree << "\n";
    // }
    // file.close();
    // cout << "n is: " << n << endl;
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
    // ofstream file("degree1.txt");
    // for(uintV i = 0; i < n; i++) {
    //   file << my_graph.V[i].getOutDegree() << endl;
    // }
    // file.close();
    // ofstream file("neighbor0.txt");
    // for(uintV v = 0; v < n; v++) {
    //   intE deg = my_graph.V[v].getOutDegree();
    //   // file << deg << endl;
    //   // continue;
    //   if (deg > 1) {
    //     std::sort(my_graph.V[v].getOutNeighbors(), my_graph.V[v].getOutNeighbors()+deg);
    //   }
      
    //   for(intE e = 0; e < deg; e++) {
    //     // cout << v << " " << my_graph.V[v].getOutNeighbor(e) << endl;
    //     file << my_graph.V[v].getOutNeighbor(e) << endl;
    //   }
    // }
    // file.close();
#ifdef EDGEWORK
    long em_work = 0;
#endif
    // ofstream file("neighbor.txt");
    // for(uintV v = 0; v < n; v++) {
    //   file << my_graph.V[v].getOutDegree() << endl; 
    // }
    // for(uintV i = n-5; i < n; i++) {
    //   intE deg = my_graph.V[i].getOutDegree();
    //   cout << "vertex: " << i << " degree: " << deg << " ngh: ";
    //   for(intE d = 0; d < deg; d++) {
    //     cout << my_graph.V[i].getOutNeighbor(d) << " ";
    //   }
    //   cout << endl;
    // }
    if (frontier_curr_vs.numNonzeros() == 0) {
      converged_iteration = start_iteration;

    } else {
      // struct perf_event_attr pe;
      // int fd;

      // memset(&pe, 0, sizeof(struct perf_event_attr));
      // pe.type = PERF_TYPE_HARDWARE;
      // pe.config = PERF_COUNT_HW_CACHE_MISSES; // Cache misses
      // pe.size = sizeof(struct perf_event_attr);
      // pe.disabled = 1; // Disabled initially
      // pe.exclude_kernel = 1;
      // pe.exclude_hv = 1;
      

      // fd = perf_event_open(&pe, 0, -1, -1, 0);
      // if (fd == -1) {
      //     perror("perf_event_open");
      //     exit(EXIT_FAILURE);
      // }
      for (iter = start_iteration; iter < max_iterations; iter++) {
        // initialize timers
        {
          phase_timer.start();
          misc_time = 0;
          copy_time = 0;
          source_change_time = 0;
          delta_time = 0;
#ifdef EDGEWORK
          em_work = 0;
#endif
        }
        MY_TIMER_LOGS([&] {
          iteration_timer.start();
          // cout << setw(PRINT_WIDTH) << iter << ",";
        });

        // ========== COPY - Prepare curr iteration ==========
        if (iter > 0) {
          // Copy the aggregate and actual value from iter-1 to iter
          parallel_for(uintV v = 0; v < n; v++) {
            vertex_values[iter][v] = vertex_values[iter - 1][v];
            aggregation_values[iter][v] = aggregation_values[iter - 1][v];
            delta[v] = aggregationValueIdentity<AggregationValueType>();
#ifdef MULTI
            source_deviation[iter][v] = sourceDeviationValueIdentity<AggregationValueType>();
#endif
          }
        }
        use_delta = shouldUseDelta(iter);
                // int activate_count = 0;
        // for(int i = 0; i < n; i++) {
        //   if(frontier_curr[i])
        //     activate_count++;
        // }
        // cout << "iter: " << iter << " activate_count:" << activate_count << endl;
        // ========== MISC - count active edges for AE ==========
        phase_time = phase_timer.next();
        adaptive_executor.updateCopyTime(iter, phase_time);
        adaptive_executor.updateEdgesProcessed(iter, my_graph,
                                               frontier_curr_vs);
        misc_time = phase_timer.next();
        adaptive_executor.updateMiscTime(iter, phase_timer.next());

        // my_time.start();

        // ========== EDGE COMPUTATION ==========
        if ((use_source_contribution) && (iter == 1)) {
          // Compute source contribution for first iteration
          parallel_for(uintV u = 0; u < n; u++) {
            if (frontier_curr[u]) {
              // compute source change in contribution
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  u, source_change_in_contribution[u],
                  vertexValueIdentity<VertexValueType>(),
                  vertex_values[iter - 1][u], global_info);
            }
          }
        }
        // source_change_time = my_time.next();
        // cout << "source_change_time: " << source_change_time << endl;

        // ioctl(fd, PERF_EVENT_IOC_RESET, 0);
        // ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
        parallel_for(uintV u = 0; u < n; u++) {
          if (frontier_curr[u]) {
            // writeAdd(&activate_vertex, (unsigned long long)1);
#ifdef EDGEWORK
            long curr_work = 0;
#endif
            // check for propagate and retract for the vertices.
            intE outDegree = my_graph.V[u].getOutDegree();
            // writeAdd(&activate_edge, (unsigned long)outDegree);
            granular_for(j, 0, outDegree, (outDegree > 1024), {
              uintV v = my_graph.V[u].getOutNeighbor(j);
              AggregationValueType contrib_change =
                  use_source_contribution
                      ? source_change_in_contribution[u]
                      : aggregationValueIdentity<AggregationValueType>();
#ifdef EDGEDATA
              EdgeData *edge_data = my_graph.V[u].getOutEdgeData(j);
#else
              EdgeData *edge_data = &emptyEdgeData;
#endif
              bool ret =
                  edgeFunction(u, v, *edge_data, vertex_values[iter - 1][u],
                               contrib_change, global_info);
              if (ret) {
#ifdef EDGEWORK
                writeAdd(&curr_work, (long)1);
#endif
                if (use_lock) {
                  vertex_locks[v].writeLock();
                  addToAggregation(contrib_change, delta[v], global_info);
                  vertex_locks[v].unlock();
                } else {
                  addToAggregationAtomic(contrib_change, delta[v], global_info);
                }
                if (!frontier_next[v])
                  frontier_next[v] = 1;
              }
            });
#ifdef EDGEWORK
            writeAdd(&em_work, curr_work);
#endif
          }
        }
        // ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
        // unsigned long long count;
        // ssize_t read_result = read(fd, &count, sizeof(count));
        // if (read_result == -1) {
        //     perror("read");
        //     exit(EXIT_FAILURE);
        // }

        // printf("Cache misses: %llu\n", count);
        // all_count += count;
        // delta_time = my_time.next();
        // cout << "delta_time: " << delta_time << endl;
        // all_delta_time += delta_time;

        phase_time = phase_timer.next();
        adaptive_executor.updateEdgeMapTime(iter, phase_time);
#ifdef EDGEWORK
        total_work_done += em_work;
#endif
#ifdef EDGEWORK
          cout << setw(PRINT_WIDTH) << em_work << ",";
          cout << setw(PRINT_WIDTH) << "0"
               << ",";
#endif

        // ========== VERTEX COMPUTATION ==========
        // for(uintV v = 0; v < n; v++) {
        parallel_for(uintV v = 0; v < n; v++) {
          // Reset frontier for next iteration
          frontier_curr[v] = 0;
          // Process all vertices affected by EdgeMap
          if (frontier_next[v] ||
              forceComputeVertexForIteration(v, iter, global_info)) {

            frontier_next[v] = 0;
            // Update aggregation value and reset change received[v] (i.e.
            // delta[v])
            addToAggregation(delta[v], aggregation_values[iter][v],
                             global_info);
            delta[v] = aggregationValueIdentity<AggregationValueType>();

            // Calculate new_value based on the updated aggregation value
            VertexValueType new_value;
            computeFunction(v, aggregation_values[iter][v],
                            vertex_values[iter - 1][v], new_value, global_info);
            // Check if change is significant
            if (notDelZero(new_value, vertex_values[iter - 1][v], global_info)) {
              // change is significant. Update vertex_values
              vertex_values[iter][v] = new_value;
              // Set active for next iteration.
              frontier_curr[v] = 1;
              

            } else {
              // change is not significant. Copy vertex_values[iter-1]
              // vertex_values[iter][v] = vertex_values[iter - 1][v];
            }
          }
          frontier_curr[v] =
              frontier_curr[v] ||
              forceActivateVertexForIteration(v, iter + 1, global_info);
          if (frontier_curr[v]) {
            if (use_source_contribution) {
              // update source_contrib for next iteration
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, source_change_in_contribution[v],
                  vertex_values[iter - 1][v], vertex_values[iter][v],
                  global_info);
            } else {
              source_change_in_contribution[v] =
                  aggregationValueIdentity<AggregationValueType>();
            }
          }
        }
        phase_time = phase_timer.stop();
        adaptive_executor.updateVertexMapTime(iter, phase_time);

        vertexSubset temp_vs(n, frontier_curr);
        frontier_curr_vs = temp_vs;
        misc_time += phase_timer.next();
        iteration_time = iteration_timer.stop();
        // MY_TIMER_LOGS([&] {
        //   cout << setw(PRINT_WIDTH) << iteration_time << "\n";
        // });
        if (ae_enabled && iter == 1) {
          adaptive_executor.setApproximateTimeForCurrIter(iteration_time);
        }

        converged_iteration = iter;
        
        // cout << iter << " " << activate_edge << endl;
        if (frontier_curr_vs.isEmpty()) {
          break;
        }
        all_activate_edge += activate_edge;
        activate_edge = 0;
        all_activate_vertex += activate_vertex;
        activate_vertex = 0;
      }
    }
    
    if (ae_enabled) {
      adaptive_executor.updateEquation(converged_iteration);
    }
    // cout << "all cache miss: " << all_count << endl;
    // cout << "all delta time: " << all_delta_time << endl;
    // static bool flag = false;
    // if (flag == false) {
    //   flag = true;
    // } else {
    //   average_cache_miss += all_count;
    //   average_delta_time += all_delta_time;
    // }

    // cout << "all activate edges: " << all_activate_edge << endl;
    // cout << "all activate vertexs: " << all_activate_vertex << endl;
    // ofstream file("neighbor.txt");
    // for(uintV v = 0; v < n; v++) {
    //   file << vertex_values[converged_iteration][v] << endl; 
    // }
    // file.close();
    
    return converged_iteration;
  }

  // ======================================================================
  // DELTACOMPUTE
  // ======================================================================
  void deltaCompute(edgeArray &edge_additions, edgeArray &edge_deletions) {
    timer iteration_timer, phase_timer, full_timer, pre_compute_timer;
    double misc_time, copy_time, phase_time, iteration_time, pre_compute_time;
    iteration_time = 0;
    full_timer.start();

    // double delta_time = 0, compute_time = 0;
    double compute_time = 0, edge_time = 0, all_copy_time;
    edge_map_time = 0, vertex_compute_time = 0, other_time = 0;
    // timer delta_timer;
#ifdef EDGEWORK
    long em_work = 0;
    long de_work = 0;
#endif
    phase_timer.start();
    // TODO : Realloc addition of new vertices
    n_old = n;
    if (edge_additions.maxVertex >= n) {
      processVertexAddition(edge_additions.maxVertex);
    }

    // Reset values before incremental computation
    parallel_for(uintV v = 0; v < n; v++) {
      frontier_curr[v] = 0;
      frontier_next[v] = 0;
      changed[v] = 0;
#ifdef TINY
      first_changed[v] = 0;
#endif
#ifdef ONE
      temp_source_deviation[v] = sourceDeviationValueIdentity<AggregationValueType>();
#endif
#ifdef MULTI
      // temp_source_deviation[v] = sourceDeviationValueIdentity<AggregationValueType>();
      is_modify_source[v] = 0;
#endif

      initializeVertexValue<VertexValueType>(v, vertex_value_old_next[v],
                                             global_info);
      delta[v] = aggregationValueIdentity<AggregationValueType>();
#ifdef TINY

#else
      // vertex_value_old_prev[v] = vertexValueIdentity<VertexValueType>();
      // vertex_value_old_curr[v] = vertexValueIdentity<VertexValueType>();
      if (use_source_contribution) {
        source_change_in_contribution[v] =
            aggregationValueIdentity<AggregationValueType>();
      }
#endif
    }

    // ==================== UPDATE GLOBALINFO ===============================
    // deltaCompute/initCompute Save a copy of global_info before we lose any
    // relevant information of the old graph For example, In PageRank, we need
    // to save the outDegree for all vertices corresponding to the old graph
    global_info_old.copy(global_info);

    // Update global_info based on edge additions or deletions. This is
    // application specific. For example, for pagerank, the the outDegree of
    // vertices with edgeAddition will increase and those with edgeDeletions
    // will decrease
    global_info.processUpdates(edge_additions, edge_deletions);
    // for(uintV v = 0; v < n; v++) {
    //   cout << "v: " << v << " activate: " << frontier_curr[v] << " changed: " << changed[v] << endl;
    // }
    global_info_old.reProcess(global_info);
    // parallel_for(uintV v = 0; v < n; v++) {
    //   frontier_next[v] = 0;
    //   frontier_curr[v] = 0;
    //   frontier_curr[v] = forceActivateVertexForIteration(v, 1, global_info);
    // }
    // traditionalIncrementalComputation(1);
    // return;
    // ========== EDGE COMPUTATION - DIRECT CHANGES - for first iter ==========
    int fal_act = 0;
    pre_compute_timer.start();
    parallel_for(long i = 0; i < edge_additions.size; i++) {
      uintV source = edge_additions.E[i].source;
      uintV destination = edge_additions.E[i].destination;

      // Update frontier and changed values
      hasSourceChangedByUpdate(source, edge_addition_enum,
                               frontier_curr[source], changed[source],
                               global_info, global_info_old);
      hasSourceChangedByUpdate(destination, edge_addition_enum,
                               frontier_curr[destination], changed[destination],
                               global_info, global_info_old);
      // cout << "has source change: " << source << " isChanged: " << changed[source] << " dest: " <<  destination << " isChanged: " << changed[destination] << endl;
#ifdef MULTI
      if(forceActivateVertexForIteration(source, 1, global_info_old) ||
          DeviationNotEqualZero(source_deviation[0][source])) {
        // if ((!forceActivateVertexForIteration(source, 1, global_info_old)) && DeviationNotEqualZero(source_deviation[0][source]))
        //   writeAdd(&fal_act, 1);
#else
      if (forceActivateVertexForIteration(source, 1, global_info_old)) {
#endif
        if (frontier_curr[source]) {
          changed[source] = true;
        }
        if (frontier_curr[destination]) {
          changed[source] = true;
        }

        AggregationValueType contrib_change;
        if (use_source_contribution) {
          sourceChangeInContribution<AggregationValueType, VertexValueType,
                                     GlobalInfoType>(
              source, contrib_change, vertexValueIdentity<VertexValueType>(),
              vertex_values[0][source], global_info_old);
#ifdef MULTI
          if (DeviationNotEqualZero(source_deviation[0][source])) {
            removeFromDeviation(source_deviation[0][source], contrib_change);
          }
#endif
        }

// Do repropagate for edge source->destination.
#ifdef EDGEDATA
        EdgeData *edge_data = edge_additions.E[i].edgeData;
#else
        EdgeData *edge_data = &emptyEdgeData;
#endif
        bool ret =
            edgeFunction(source, destination, *edge_data,
                         vertex_values[0][source], contrib_change, global_info);
        if (ret) {
#ifdef EDGEWORK
          writeAdd(&de_work, (long)1);
#endif
          if (use_lock) {
            vertex_locks[destination].writeLock();
            addToAggregation(contrib_change, delta[destination],
                             global_info_old);
            vertex_locks[destination].unlock();
          } else {
            addToAggregationAtomic(contrib_change, delta[destination],
                                   global_info_old);
          }
          if (!changed[destination]) {
            changed[destination] = true;
            // cout << "...changed: " << destination << endl;
          }
          

        }
      }
#ifdef MULTI
      is_modify_source[source] = 1;
#endif
    }

    parallel_for(long i = 0; i < edge_deletions.size; i++) {
      uintV source = edge_deletions.E[i].source;
      uintV destination = edge_deletions.E[i].destination;

      hasSourceChangedByUpdate(source, edge_deletion_enum,
                               frontier_curr[source], changed[source],
                               global_info, global_info_old);
      hasSourceChangedByUpdate(destination, edge_deletion_enum,
                               frontier_curr[destination], changed[destination],
                               global_info, global_info_old);
#ifdef MULTI
      if(forceActivateVertexForIteration(source, 1, global_info_old) ||
          DeviationNotEqualZero(source_deviation[0][source])) {
        // if ((!forceActivateVertexForIteration(source, 1, global_info_old)) && DeviationNotEqualZero(source_deviation[0][source]))
        //   writeAdd(&fal_act, 1);
#else
      if (forceActivateVertexForIteration(source, 1, global_info_old)) {
#endif
        // Update frontier and changed values
        if (frontier_curr[source]) {
          changed[source] = true;
        }
        if (frontier_curr[destination]) {
          changed[source] = true;
        }

        AggregationValueType contrib_change;
        if (use_source_contribution) {
          sourceChangeInContribution<AggregationValueType, VertexValueType,
                                     GlobalInfoType>(
              source, contrib_change, vertexValueIdentity<VertexValueType>(),
              vertex_values[0][source], global_info_old);
#ifdef MULTI
          if (DeviationNotEqualZero(source_deviation[0][source])) {
            removeFromDeviation(source_deviation[0][source], contrib_change);
          }
#endif
        }

// Do retract for edge source->destination
#ifdef EDGEDATA
        EdgeData *edge_data = edge_deletions.E[i].edgeData;
#else
        EdgeData *edge_data = &emptyEdgeData;
#endif
        bool ret = edgeFunction(source, destination, *edge_data,
                                vertex_values[0][source], contrib_change,
                                global_info_old);
        if (ret) {
#ifdef EDGEWORK
          writeAdd(&de_work, (long)1);
#endif
          if (use_lock) {
            vertex_locks[destination].writeLock();
            removeFromAggregation(contrib_change, delta[destination],
                                  global_info_old);
            vertex_locks[destination].unlock();
          } else {
            removeFromAggregationAtomic(contrib_change, delta[destination],
                                        global_info_old);
          }
          if (!changed[destination])
            changed[destination] = true;

        }
      }
#ifdef MULTI
      is_modify_source[source] = 1;
#endif
    }
    pre_compute_time = pre_compute_timer.stop();

    // =============== INCREMENTAL COMPUTE - REFINEMENT START ================
//     MY_TIMER_LOGS([&] {
//       cout << "\n"
//            << setw(PRINT_WIDTH + 1) << "Iteration,";
// #ifdef EDGEWORK
//       cout << setw(PRINT_WIDTH + 1) << "T_Edges," << setw(PRINT_WIDTH + 1)
//            << "D_Edges,";
// #endif
//       cout << setw(PRINT_WIDTH + 1) << "vertices,";
//       cout << setw(PRINT_WIDTH + 1) << "edges,";
//       cout << setw(PRINT_WIDTH + 1) << "Time\n";

//     });

    vertexSubset frontier_curr_vs(n, frontier_curr);
    bool should_switch_now = false;
    bool use_delta = true;

    if (ae_enabled && shouldSwitch(0, 0)) {
      should_switch_now = true;
    }
    int delta_activate_vertex_count = 0;
    int delta_activate_edge_count = 0;
    // unsigned long all_activate_edges = 0;
    // unsigned long iter_activate_edges = 0;
    unsigned long long activate_vertexs = 0;
    unsigned long long all_activate_vertex = 0;
    all_activate_edge = 0;
    activate_edge = 0;

    // parallel_for(uintV u = 0; u < n; u++) {
    //   delta[u] = aggregationValueIdentity<AggregationValueType>();
    //   source_change_in_contribution[u] = aggregationValueIdentity<AggregationValueType>();
    //   frontier_curr[u] = forceActivateVertexForIteration(u, 1, global_info);;
    //   frontier_next[u] = 0;
    // }
    // traditionalIncrementalComputation(1);
    // return ;

    // delta_timer.start();
    // cout << "-------------------" << endl;
    // for(uintV v = 0; v < n; v++) {
    //   // cout << "v: " << v << " activate: " << frontier_curr[v] << " changed: " << changed[v] << endl;
    // }
    // cout << "false activate vertex count: " << fal_act << endl;
    other_time = phase_timer.stop();
    for (int iter = 1; iter < max_iterations; iter++) {
      activate_edge = 0;
      activate_vertexs = 0;
      // Perform switch if needed
      if (should_switch_now) {
        converged_iteration = performSwitch(iter);
        break;
      }

      // initialize timers
      {
        iteration_timer.start();
        phase_timer.start();
        iteration_time = 0;
        misc_time = 0;
        copy_time = 0;
#ifdef EDGEWORK
        em_work = 0;
#endif

      }
      use_delta = shouldUseDelta(iter);

      // ================ COPY - PREPARE CURRENT ITERATION ================
      {
#ifdef TINY
        if(iter == 1) {
          VertexValueType *temp1 = vertex_value_old_prev;
          vertex_value_old_prev = vertex_value_old_curr;
          vertex_value_old_curr = vertex_value_old_next;
          vertex_value_old_next = temp1;
        }
#else
        VertexValueType *temp1 = vertex_value_old_prev;
        vertex_value_old_prev = vertex_value_old_curr;
        vertex_value_old_curr = vertex_value_old_next;
        vertex_value_old_next = temp1;
#endif
        if (iter <= converged_iteration) {
          parallel_for(uintV v = 0; v < n; v++) {
#ifdef TINY
#else
            vertex_value_old_next[v] = vertex_values[iter][v];
#endif
            if(is_modify_source[v]) {
              del_source_deviation[v] = source_deviation[iter][v];
            }
          }
        } else {
          converged_iteration = performSwitch(iter);
          break;
        }
      }
      // MY_TIMER_LOGS([&] { cout << setw(PRINT_WIDTH) << iter << ","; });
      // copy_time += phase_timer.next();
      other_time += phase_timer.next();
      // all_copy_time += copy_time;
      // ========== EDGE COMPUTATION - TRANSITIVE CHANGES ==========

      if ((use_source_contribution) && (iter == 1)) {
        // Compute source contribution for first iteration
        parallel_for(uintV u = 0; u < n; u++) {
          if (frontier_curr[u]) {
            // compute source change in contribution
#ifdef MULTI
            bool old_change, new_change;
#ifdef TINY
            source_change_in_contribution[u] =
                aggregationValueIdentity<AggregationValueType>();
#endif
            old_change = notDelZero(vertex_value_old_curr[u], 
                                vertexValueIdentity<VertexValueType>(), global_info_old);
            new_change = notDelZero(vertex_values[iter-1][u],
                                vertexValueIdentity<VertexValueType>(), global_info);
            AggregationValueType contrib_change =
                aggregationValueIdentity<AggregationValueType>();
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                        GlobalInfoType>(u, contrib_change, 
                  vertexValueIdentity<VertexValueType>(), vertex_value_old_curr[u],
                  vertexValueIdentity<VertexValueType>(), vertex_values[iter-1][u],
                  global_info_old, global_info, new_change, old_change);
            // sourceChangeInContribution<AggregationValueType, VertexValueType,
            //                             GlobalInfoType>(u, contrib_change, 
            //       vertex_value_old_curr[u], vertexValueIdentity<VertexValueType>(),
            //       vertex_values[iter-1][u], vertexValueIdentity<VertexValueType>(),
            //       global_info_old, global_info, new_change, old_change);
            addToDeviation(contrib_change, source_deviation[iter-1][u]);
            if(sourceNotDelZero(u, source_deviation[iter-1][u],
                                aggregationValueIdentity<AggregationValueType>(),
                                global_info)) {
              source_change_in_contribution[u] = source_deviation[iter-1][u];
              source_deviation[iter-1][u] = sourceDeviationValueIdentity<AggregationValueType>();
            } else {
              frontier_curr[u] = 0;
            }
#endif
            /*
            // AggregationValueType contrib_change =
            //     aggregationValueIdentity<AggregationValueType>();
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                u, contrib_change, vertexValueIdentity<VertexValueType>(),
                vertex_values[iter - 1][u], global_info);
            addToAggregation(contrib_change, source_change_in_contribution[u],
                             global_info);
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                u, contrib_change, vertexValueIdentity<VertexValueType>(),
                vertex_value_old_curr[u], global_info_old);
            removeFromAggregation(
                contrib_change, source_change_in_contribution[u], global_info);
            */
          }
        }
      }
      // delta_timer.next();
      parallel_for(uintV u = 0; u < n; u++) {
        if (frontier_curr[u]) {
          // writeAdd(&activate_vertexs, (unsigned long long)1);
          // check for propagate and retract for the vertices.
          intE outDegree = my_graph.V[u].getOutDegree();
#ifdef COUNT
          writeAdd(&activate_edge, (unsigned long)outDegree);
#endif
#ifdef EDGEWORK
          long curr_work = 0;
#endif
          granular_for(i, 0, outDegree, (outDegree > 1024), {
            uintV v = my_graph.V[u].getOutNeighbor(i);
            bool ret = false;
            AggregationValueType contrib_change =
                use_source_contribution
                    ? source_change_in_contribution[u]
                    : aggregationValueIdentity<AggregationValueType>();

#ifdef EDGEDATA
            EdgeData *edge_data = my_graph.V[u].getOutEdgeData(i);
#else
            EdgeData *edge_data = &emptyEdgeData;
#endif
            ret = edgeFunction(u, v, *edge_data, vertex_values[iter - 1][u],
                               contrib_change, global_info);

            if (ret) {
#ifdef EDGEWORK
              writeAdd(&curr_work, (long)1);
#endif
              if (use_lock) {
                vertex_locks[v].writeLock();
                if (ret) {
                  addToAggregation(contrib_change, delta[v], global_info);
                }
                vertex_locks[v].unlock();

              } else {
                if (ret) {
                  addToAggregationAtomic(contrib_change, delta[v], global_info);
                }
              }

              if (!changed[v])
                changed[v] = 1;
            }
          });
#ifdef EDGEWORK
          writeAdd(&em_work, curr_work);
#endif
        }
      }
      // delta_time += delta_timer.next();
      phase_time = phase_timer.next();
      edge_map_time += phase_time;
      // edge_time += phase_time;
#ifdef EDGEWORK
        total_work_done += em_work;
        total_work_done += de_work;
        cout << setw(PRINT_WIDTH) << em_work << ",";
        cout << setw(PRINT_WIDTH) << de_work << ",";
        de_work = 0;
#endif

      // ========== VERTEX COMPUTATION  ==========
      bool use_delta_next_iteration = shouldUseDelta(iter + 1);
      parallel_for(uintV v = 0; v < n; v++) {
        // changed vertices need to be processed
        frontier_curr[v] = 0;
        if ((v >= n_old) && (changed[v] == false)) {
          changed[v] = forceComputeVertexForIteration(v, iter, global_info);
        }

        if (changed[v]) {
          // cout << "changedv: " << v << endl;
          frontier_curr[v] = 0;

          // delta has the current cumulative change for the vertex.
          // Update the aggregation value in history
          addToAggregation(delta[v], aggregation_values[iter][v], global_info);

          VertexValueType new_value;
          computeFunction(v, aggregation_values[iter][v],
                          vertex_values[iter - 1][v], new_value, global_info);

          if (forceActivateVertexForIteration(v, iter + 1, global_info)) {
            frontier_curr[v] = 1;
          }
          AggregationValueType contrib_change =
              aggregationValueIdentity<AggregationValueType>();
          source_change_in_contribution[v] =
              aggregationValueIdentity<AggregationValueType>();
          bool old_change = 0, new_change = 0;
#ifdef TINY
          if(!first_changed[v]) {
            vertex_value_old_curr[v] = vertex_values[iter-1][v];
            first_changed[v] = 1;
          } else {
            vertex_value_old_curr[v] = vertex_value_old_next[v];
          }
          
          vertex_value_old_next[v] = vertex_values[iter][v];
#endif
          if (notDelZero(new_value, vertex_values[iter - 1][v], global_info)) {
            // change is significant. Update vertex_values
            vertex_values[iter][v] = new_value;
            new_change = 1;
            frontier_curr[v] = 1;
#ifdef MULTI
            // cout << "ignore0...." << endl;
            //--------------------------
            if (use_delta_next_iteration) {
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, contrib_change, vertex_values[iter - 1][v],
                  vertex_values[iter][v], global_info);
            } else {
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, contrib_change, vertexValueIdentity<VertexValueType>(),
                  vertex_values[iter][v], global_info);
            }
            // addToAggregation(contrib_change, source_change_in_contribution[v],
            //                  global_info);
            addToDeviation(contrib_change, source_deviation[iter][v]);

            //--------------------------------------
#endif
          } else {
            // change is not significant. Copy vertex_values[iter-1]
            vertex_values[iter][v] = vertex_values[iter - 1][v];
          }
          

          if (notDelZero(vertex_value_old_next[v], vertex_value_old_curr[v],
                        global_info_old)) {
            // change is significant. Update v_change
            frontier_curr[v] = 1;
            old_change = 1;
#ifdef MULTI
            //===============================
            // cout << "ignore1...." << endl;
            if (use_delta_next_iteration) {
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, contrib_change, vertex_value_old_curr[v],
                  vertex_value_old_next[v], global_info_old);
            } else {

              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, contrib_change, vertexValueIdentity<VertexValueType>(),
                  vertex_value_old_next[v], global_info_old);
            }
            // removeFromAggregation(contrib_change,
            //                       source_change_in_contribution[v],
            //                       global_info_old);
            removeFromDeviation(contrib_change, source_deviation[iter][v]);
            //-----------------------------------
#endif
          }
#ifdef ONE
          // cout << "one aggregate deviation....." << endl;
          if (frontier_curr[v]) {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(v, contrib_change,
                  vertex_value_old_curr[v], vertex_value_old_next[v],
                  vertex_values[iter-1][v], vertex_values[iter][v],
                  global_info_old, global_info, new_change, old_change);
            addToDeviation(contrib_change, temp_source_deviation[v]);
            if(sourceNotDelZero(v, temp_source_deviation[v],
                                aggregationValueIdentity<AggregationValueType>(),
                                global_info)) {
              source_change_in_contribution[v] = temp_source_deviation[v];
              temp_source_deviation[v] = sourceDeviationValueIdentity<AggregationValueType>();
            } else {
              frontier_curr[v] = 0;
            }
          }
#endif

#ifdef MULTI
          if (frontier_curr[v]) {
            // sourceChangeInContribution<AggregationValueType, VertexValueType,
            //                              GlobalInfoType>(v, contrib_change,
            //       vertex_value_old_curr[v], vertex_value_old_next[v],
            //       vertex_values[iter-1][v], vertex_values[iter][v],
            //       global_info_old, global_info, new_change, old_change);
            // addToDeviation(contrib_change, source_deviation[iter][v]);
            if(sourceNotDelZero(v, source_deviation[iter][v],
                                aggregationValueIdentity<AggregationValueType>(),
                                global_info)) {
              source_change_in_contribution[v] = source_deviation[iter][v];
              source_deviation[iter][v] = sourceDeviationValueIdentity<AggregationValueType>();
            } else {
              frontier_curr[v] = 0;
            }
          }
#endif
          //-----------------------------------------------------------------
#ifdef IGNORE
          // cout << "ignore0...." << endl;
          if(frontier_curr[v] &&!sourceNotDelZero(v, source_change_in_contribution[v],
                                  aggregationValueIdentity<AggregationValueType>(),
                                  global_info)) {
              frontier_curr[v] = 0;
          }
#endif
          //--------------------------------------------------------------------
        }
      }
      // compute_time += delta_timer.next();
      phase_time = phase_timer.next();
      vertex_compute_time += phase_time;
      // compute_time += phase_time;
      // ========== EDGE COMPUTATION - DIRECT CHANGES - for next iter ==========

      bool has_direct_changes = false;
      parallel_for(long i = 0; i < edge_additions.size; i++) {
        uintV source = edge_additions.E[i].source;
        uintV destination = edge_additions.E[i].destination;
        AggregationValueType contrib_change;
#ifdef MULTI
  #ifdef TINY
        VertexValueType old_val = changed[source] ? vertex_value_old_curr[source] : vertex_values[iter-1][source];
        VertexValueType new_val = changed[source] ? vertex_value_old_next[source] : vertex_values[iter][source];
        if (notDelZero(old_val, new_val, global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old)) ||
            DeviationNotEqualZero(del_source_deviation[source])) {
  #else
        if (notDelZero(vertex_value_old_curr[source],
                      vertex_value_old_next[source], global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old)) ||
            DeviationNotEqualZero(del_source_deviation[source])) {
  #endif
#else
        if (notDelZero(vertex_value_old_curr[source],
                      vertex_value_old_next[source], global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old))) {
#endif
          if (use_delta_next_iteration) {
#ifdef TINY
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, old_val,
                new_val, global_info_old);
#else
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, vertex_value_old_curr[source],
                vertex_value_old_next[source], global_info_old);
#endif
#ifdef MULTI
            if (DeviationNotEqualZero(del_source_deviation[source])) {
              removeFromDeviation(del_source_deviation[source], contrib_change);
            }
#endif
          } else {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, vertexValueIdentity<VertexValueType>(),
                vertex_value_old_next[source], global_info_old);
          }
// Do repropagate for edge source->destination.
#ifdef EDGEDATA
          EdgeData *edge_data = edge_additions.E[i].edgeData;
#else
          EdgeData *edge_data = &emptyEdgeData;
#endif
          bool ret = edgeFunction(source, destination, *edge_data,
                                  vertex_values[0][source], contrib_change,
                                  global_info);

          if (ret) {
#ifdef EDGEWORK
            writeAdd(&de_work, (long)1);
#endif
            if (use_lock) {
              vertex_locks[destination].writeLock();
              addToAggregation(contrib_change, delta[destination],
                               global_info_old);
              vertex_locks[destination].unlock();
            } else {
              addToAggregationAtomic(contrib_change, delta[destination],
                                     global_info_old);
            }
            if (!changed[destination])
              changed[destination] = 1;
            if (!has_direct_changes)
              has_direct_changes = true;
          }
        }
      }

      parallel_for(long i = 0; i < edge_deletions.size; i++) {
        uintV source = edge_deletions.E[i].source;
        uintV destination = edge_deletions.E[i].destination;
        AggregationValueType contrib_change;
#ifdef MULTI
  #ifdef TINY
        VertexValueType old_val = changed[source] ? vertex_value_old_curr[source] : vertex_values[iter-1][source];
        VertexValueType new_val = changed[source] ? vertex_value_old_next[source] : vertex_values[iter][source];
        if (notDelZero(old_val, new_val, global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old)) ||
            DeviationNotEqualZero(del_source_deviation[source])) {
  #else
        if (notDelZero(vertex_value_old_curr[source],
                      vertex_value_old_next[source], global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old)) ||
            DeviationNotEqualZero(del_source_deviation[source])) {
  #endif
#else
        if (notDelZero(vertex_value_old_curr[source],
                      vertex_value_old_next[source], global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old))) {
#endif
          // Do repropagate for edge source->destination.
          if (use_delta_next_iteration) {
#ifdef TINY
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, old_val,
                new_val, global_info_old);
#else
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, vertex_value_old_curr[source],
                vertex_value_old_next[source], global_info_old);
#endif
#ifdef MULTI
            if(DeviationNotEqualZero(del_source_deviation[source])) {
              removeFromDeviation(del_source_deviation[source], contrib_change);
            }
#endif
          } else {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change, vertexValueIdentity<VertexValueType>(),
                vertex_value_old_next[source], global_info_old);
          }
#ifdef EDGEDATA
          EdgeData *edge_data = edge_deletions.E[i].edgeData;
#else
          EdgeData *edge_data = &emptyEdgeData;
#endif
          bool ret = edgeFunction(source, destination, *edge_data,
                                  vertex_values[0][source], contrib_change,
                                  global_info);

          if (ret) {
#ifdef EDGEWORK
            writeAdd(&de_work, (long)1);
#endif
            if (use_lock) {
              vertex_locks[destination].writeLock();
              removeFromAggregation(contrib_change, delta[destination],
                                    global_info_old);
              vertex_locks[destination].unlock();

            } else {
              removeFromAggregationAtomic(contrib_change, delta[destination],
                                          global_info_old);
            }
            if (!changed[destination])
              changed[destination] = 1;
            if (!has_direct_changes)
              has_direct_changes = true;
          }
        }
      }
      phase_time = phase_timer.next();
      other_time += phase_time;
#ifdef COUNT
      cout << iter << "\t" << activate_edge << endl;
#endif

      misc_time += phase_timer.next();
      iteration_time = iteration_timer.next();

      // Convergence check
      
      if(iter == converged_iteration) {
        vertexSubset temp_vs(n, frontier_curr);
        frontier_curr_vs = temp_vs;
        if (!has_direct_changes && frontier_curr_vs.isEmpty()) {
          break;
        }
      }
      /*
      if (!has_direct_changes && frontier_curr_vs.isEmpty()) {
        // There are no more active vertices
        if (iter == converged_iteration) {
          // MY_TIMER_LOGS([&] {
          //   cout << setw(PRINT_WIDTH) << iteration_time << "\n";
          // });
          break;
        } else if (iter > converged_iteration) {
          assert(("Missed switching to Traditional incremental computing when "
                  "iter == converged_iter",
                  false));
        } else {
          // Values stable for the changed vertices at this iteration.
          // But, the changed vertices might receive new changes. So,
          // continue loop until iter == converged_iteration vertices may
          // still not have converged. So, keep continuing until
          // converged_iteration is reached.
        }
      }*/
      if (iter == 1) {
        iteration_time += pre_compute_time;
      }
#ifdef COUNT
      all_activate_edge += activate_edge; 
      // all_activate_vertex += activate_vertexs;
#endif
      if (ae_enabled && shouldSwitch(iter, iteration_time)) {
        should_switch_now = true;
      }
      misc_time += phase_timer.stop();
      iteration_time += iteration_timer.stop();
      // MY_TIMER_LOGS([&] {
      //   cout << setw(PRINT_WIDTH) << iteration_time << "\n";
      // });
    }
    // cout << "compute time: " << compute_time << endl;
    // cout << "all delta time: " << delta_time << endl;

    // cout << "vertex copy time: " << all_copy_time << endl;
    // cout << "vertex compute time: " << compute_time << endl;
#ifdef SUBTIME
    cout << "vertex compute time: " << vertex_compute_time << endl;
    cout << "edge computing time: " << edge_map_time << endl;
    cout << "other time: " << other_time << endl;
    all_edge_map_time += edge_map_time;
    all_other_time += other_time;
    all_vertex_compute_time += vertex_compute_time;
#endif
#ifdef COUNT
    cout << "all activate edges: " << all_activate_edge << endl;
#endif
    // cout << "all activate vertexs: " << all_activate_vertex << endl;
    double full_time = full_timer.stop();
    average_delta_compute_time += full_time;
    cout << "Finished batch : " << full_time << "\n";
    // cout << "Finished batch : " << full_timer.stop() << "\n";
#ifdef EDGEWORK
    cout << "Edges processed : " << total_work_done << "\n";
#endif
    cout << "Number of iterations : " << converged_iteration << "\n";
    first_delta_compute = false;
    // testPrint();
    // printOutput();
  }

  // Refactor this in a better way
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::my_graph;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::config;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::max_iterations;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::history_iterations;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::converged_iteration;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::use_lock;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_locks;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::aggregation_values;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_values;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::n;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::global_info;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::delta;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::use_source_contribution;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::source_change_in_contribution;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::n_old;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::global_info_old;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_value_old_next;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_value_old_curr;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_value_old_prev;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::all;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::frontier_curr;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::frontier_next;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::changed;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::ingestor;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::current_batch;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::adaptive_executor;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::ae_enabled;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::total_work_done;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::testPrint;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::printOutput;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::shouldSwitch;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::performSwitch;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::processVertexAddition;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::average_cache_miss;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::average_delta_time;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::average_delta_compute_time;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::all_activate_edge;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::activate_edge;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::first_delta_compute;
#ifdef MULTI
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::source_deviation;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::temp_source_deviation;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::is_modify_source;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::del_source_deviation;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::edge_map_time;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_compute_time;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::other_time;
#endif
#ifdef ONE
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::temp_source_deviation;
#endif
#ifdef TINY
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::first_changed;
#endif
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::all_edge_map_time;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::all_other_time;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::all_vertex_compute_time;
};
#endif
