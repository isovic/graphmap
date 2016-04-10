/*
 * vertices.h
 *
 *  Created on: Feb 13, 2015
 *      Author: isovic
 */

#ifndef VERTICES_H_
#define VERTICES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include "log_system/log_system.h"



// Quite an ugly data structure, but very cache friendly.
class Vertices {
 public:
  int64_t *timestamps;
  int64_t *reference_starts;
  int64_t *reference_ends;
  int64_t *query_starts;
  int64_t *query_ends;
  int64_t *num_kmers;
  int64_t *covered_bases_queries;
  int64_t *covered_bases_references;
  int64_t *registry_numbers;

  int64_t num_vertices;
  int64_t container_capacity;

  Vertices();
  ~Vertices();
  void Clear();

  inline int Init(int64_t dest_vertex_idx, int64_t timestamp, int64_t reference_start,
           int64_t query_start, int64_t kmer_length, int64_t registry_number) {
    return Set(dest_vertex_idx, timestamp, reference_start, reference_start, query_start, query_start, 1, kmer_length, kmer_length, registry_number);
  }

  inline int Set(int64_t dest_vertex_idx, int64_t timestamp, int64_t reference_start,
          int64_t reference_end, int64_t query_start, int64_t query_end,
           int64_t num_kmer, int64_t covered_bases_query,
           int64_t covered_bases_reference, int64_t registry_number) {
    if (dest_vertex_idx >= num_vertices || dest_vertex_idx < 0) {
      return 1;
    }

    timestamps[dest_vertex_idx] = timestamp;
    reference_starts[dest_vertex_idx] = reference_start;
    reference_ends[dest_vertex_idx] = reference_end;
    query_starts[dest_vertex_idx] = query_start;
    query_ends[dest_vertex_idx] = query_end;
    num_kmers[dest_vertex_idx] = num_kmer;
    covered_bases_queries[dest_vertex_idx] = covered_bases_query;
    covered_bases_references[dest_vertex_idx] = covered_bases_reference;
    registry_numbers[dest_vertex_idx] = registry_number;

    return 0;
  }

  inline int Add(int64_t timestamp, int64_t reference_start,
           int64_t reference_end, int64_t query_start, int64_t query_end,
           int64_t num_kmer, int64_t covered_bases_query,
           int64_t covered_bases_reference, int64_t registry_number) {
    if (num_vertices >= container_capacity) {
      Reserve(container_capacity + capacity_increment_size_);
    }

    num_vertices += 1;
    Set((num_vertices - 1), timestamp, reference_start, reference_end, query_start, query_end, num_kmer, covered_bases_query, covered_bases_reference, registry_number);

    return 0;
  }

  inline int Add(const Vertices &src_vertices, int64_t src_vertex_idx) {
    return Add(src_vertices.timestamps[src_vertex_idx],
                src_vertices.reference_starts[src_vertex_idx],
                src_vertices.reference_ends[src_vertex_idx],
                src_vertices.query_starts[src_vertex_idx],
                src_vertices.query_ends[src_vertex_idx],
                src_vertices.num_kmers[src_vertex_idx],
                src_vertices.covered_bases_queries[src_vertex_idx],
                src_vertices.covered_bases_references[src_vertex_idx],
                src_vertices.registry_numbers[src_vertex_idx]);
  }

  void Reserve(int64_t size);
  void Resize(int64_t size);

  inline int CopyValuesWithin(int64_t source_idx, int64_t dest_idx) {
    if (source_idx >= num_vertices || dest_idx >= num_vertices || source_idx < 0 || dest_idx < 0) {
      LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "When CopyValuesWithin is called. source_idx = %ld, dest_idx = %ld, num_vertices = %ld\n", source_idx, dest_idx, num_vertices));
      return 1;
    }

    timestamps[dest_idx] = timestamps[source_idx];
    reference_starts[dest_idx] = reference_starts[source_idx];
    reference_ends[dest_idx] = reference_ends[source_idx];
    query_starts[dest_idx] = query_starts[source_idx];
    query_ends[dest_idx] = query_ends[source_idx];
    num_kmers[dest_idx] = num_kmers[source_idx];
    covered_bases_queries[dest_idx] = covered_bases_queries[source_idx];
    covered_bases_references[dest_idx] = covered_bases_references[source_idx];
    registry_numbers[dest_idx] = registry_numbers[source_idx];

    return 0;
  }

  inline int CopyValuesFromOut(Vertices &src_vertices, int64_t src_vertex_idx, int64_t dest_idx) {
    return Set(dest_idx,
               src_vertices.timestamps[src_vertex_idx],
               src_vertices.reference_starts[src_vertex_idx],
               src_vertices.reference_ends[src_vertex_idx],
               src_vertices.query_starts[src_vertex_idx],
               src_vertices.query_ends[src_vertex_idx],
               src_vertices.num_kmers[src_vertex_idx],
               src_vertices.covered_bases_queries[src_vertex_idx],
               src_vertices.covered_bases_references[src_vertex_idx],
               src_vertices.registry_numbers[src_vertex_idx]);
  }

  inline void EraseValues() {
    if (num_vertices <= 0)
      return;

    memset(timestamps, -1, num_vertices);
    memset(reference_starts, 0, num_vertices);
    memset(reference_ends, 0, num_vertices);
    memset(query_starts, 0, num_vertices);
    memset(query_ends, 0, num_vertices);
    memset(num_kmers, 0, num_vertices);
    memset(covered_bases_queries, 0, num_vertices);
    memset(covered_bases_references, 0, num_vertices);
    memset(registry_numbers, -1, num_vertices);
  }

  inline float CalculateRatio(int64_t vertex_idx) {
    float ratio = 0.0f;
    int64_t query_start = query_starts[vertex_idx];
    int64_t query_end = query_ends[vertex_idx];
    int64_t reference_start = reference_starts[vertex_idx];
    int64_t reference_end = reference_ends[vertex_idx];

    int64_t query_distance = (query_end >= query_start) ? (query_end - query_start) : (query_start - query_end);
    int64_t ref_distance = (reference_end >= reference_start) ? (reference_end - reference_start) : (reference_start - reference_end);

    if (query_distance != 0)
      ratio = ((float) std::min(query_distance, ref_distance)) / ((float) std::max(query_distance, ref_distance));
    else
      ratio = 1.0f;

    return ratio;
  }

  inline float CalculateSuppress(int64_t vertex_idx) {
    float ratio = 0.0f, ratio_suppress = 0.0f;

    ratio = CalculateRatio(vertex_idx);

    ratio_suppress = (ratio < 1.0f) ? (1.0f - ratio) : (ratio - 1.0f);

    return ratio_suppress;
  }

  inline std::string VerboseToString(int64_t vertex_idx) const {
    std::stringstream ret;

    if (vertex_idx < 0 || vertex_idx >= num_vertices) {
      ret << "Error with vertex_idx! vertex_idx = " << vertex_idx << ", containter_capacity = " << container_capacity << ", num_vertices = " << num_vertices;
      return ret.str();
    }

    ret << "timestamp = " << timestamps[vertex_idx];
    ret <<  "; q[" << query_starts[vertex_idx] << ", " << query_ends[vertex_idx] << "]; r[" << reference_starts[vertex_idx]<< ", " <<
        reference_ends[vertex_idx] <<
            "]; d[" << (query_ends[vertex_idx] - query_starts[vertex_idx]) << ", " << (reference_ends[vertex_idx] - reference_starts[vertex_idx]) <<
            "]; length = " << num_kmers[vertex_idx] <<
            "; dist_ratio = " << ((double) std::min((reference_ends[vertex_idx] - reference_starts[vertex_idx]), (query_ends[vertex_idx] - query_starts[vertex_idx]))) / ((double) std::max((reference_ends[vertex_idx] - reference_starts[vertex_idx]), (query_ends[vertex_idx] - query_starts[vertex_idx]))) <<
            "; cov_bases_query = " << covered_bases_queries[vertex_idx] << "; cov_bases_ref = " << covered_bases_references[vertex_idx] << "; registry_num = " << registry_numbers[vertex_idx];

    return ret.str();
  }

 private:
  inline int ReallocArray_(int64_t **array_ptr, int64_t size);
  int64_t capacity_increment_size_;
};

#endif /* VERTICES_H_ */
