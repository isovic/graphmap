/*
 * filter_anchors.h
 *
 *  Created on: Mar 22, 2016
 *      Author: isovic
 */

#ifndef SRC_GRAPHMAP_FILTER_ANCHORS_H_
#define SRC_GRAPHMAP_FILTER_ANCHORS_H_

#include <string>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include "sequences/single_sequence.h"
#include "sequences/sequence_file.h"
#include "containers/vertices.h"
#include "program_parameters.h"

#include "containers/score_registry.h"
#include "utility/utility_general.h"
#include "containers/region.h"
#include "containers/mapping_data.h"
#include "containers/vertices.h"

/// These are some constants used for filtering shady anchors.
/// TODO: This can be omitted if dynamic programming was used to penalize the anchor distances.
/// int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * MIN_CLUSTER_COVERAGE_FACTOR;
#define MIN_CLUSTER_COVERAGE_FACTOR 0.05f
/// int64_t min_cluster_length = read->get_sequence_length() * MIN_CLUSTER_LENGTH_FACTOR;
#define MIN_CLUSTER_LENGTH_FACTOR 0.03f

using int128_t = __int128;
using uint128_t = unsigned __int128;

#define get128_qid(x)   ((int32_t) (x & 0x0FFFFFFFF))
#define get128_rpos(x)  ((int32_t) ((x >> 32) & 0x0FFFFFFFF))
#define get128_qpos(x)  ((int32_t) ((x >> 64) & 0x0FFFFFFFF))
#define get128_rid(x)   ((int32_t) ((x >> 96) & 0x0FFFFFFFF))
///      d                  c                 b             a
/// ref_id << 96 | query_start << 64 | ref_start << 32 | query_id
#define pack128(qstart,rstart,qid,rid)         ((((uint128_t) rid) << 96) | (((uint128_t) qstart) << 64) | (((uint128_t) rstart) << 32) | ((uint128_t) qid))

struct ClusterAndIndices {
  Range query;
  Range ref;
  int32_t num_anchors = 0;
  int32_t coverage = 0;
  std::vector<int> lcskpp_indices;
};

int64_t CalcScore(int32_t qpos, int32_t rpos, int32_t next_qpos, int32_t next_rpos, double indel_bandwidth_margin, int32_t fwd_length, int32_t dist_aab, int32_t dist_dbm, double *score_gap, double *score_dist);

void GetPositionsFromRegistry2(const Vertices& registry_entries, int64_t vertex_id, int32_t *qpos_start, int32_t *rpos_start, int32_t *qpos_end, int32_t *rpos_end);
void GetPositionsFromRegistry(const Vertices& registry_entries, const std::vector<int> &lcskpp_indices, int64_t lcskpp_id, int32_t *qpos_start, int32_t *rpos_start, int32_t *qpos_end, int32_t *rpos_end);
void GetPositionsFrom128bit(const std::vector<uint128_t> &hits, const std::vector<int> &lcskpp_indices, int64_t lcskpp_id, int32_t seed_len, int32_t *qpos_start, int32_t *rpos_start, int32_t *qpos_end, int32_t *rpos_end);

int FilterAnchorsByDiff(const SingleSequence* read, ScoreRegistry* local_score, const ProgramParameters *parameters,
                  const std::vector<int> &lcskpp_indices, std::vector<int> &ret_filtered_lcskpp_indices);

int FilterAnchorsByChaining(const SingleSequence* seq, ScoreRegistry* local_score, const ProgramParameters *parameters,
                  const std::vector<int> &lcskpp_indices, double indel_bandwidth_margin, int32_t max_dist, int32_t lookahead_dist_factor, int64_t min_covered_bases, int32_t cluster_size_cutoff,
                  std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids);

int GenerateClusters(int64_t min_num_anchors_in_cluster, int64_t min_cluster_length, int64_t min_cluster_covered_bases, float min_cluster_coverage, std::vector<int> &lcskpp_indices,
                     ScoreRegistry* local_score, MappingData* mapping_data,
                     const SingleSequence* read, const ProgramParameters* parameters, std::vector<ClusterAndIndices *> &ret_clusters,
                     std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids);
int GenerateClustersDummy(int64_t min_cluster_length, float min_cluster_coverage, std::vector<int> &lcskpp_indices,
                     ScoreRegistry* local_score, MappingData* mapping_data,
                     const SingleSequence* read, const ProgramParameters* parameters, std::vector<ClusterAndIndices *> &ret_clusters,
                     std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids);

int VerboseClustersToFile_(std::string out_file, const ScoreRegistry* local_score, const MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters, const std::vector<ClusterAndIndices *> &clusters);

#endif /* SRC_GRAPHMAP_FILTER_ANCHORS_H_ */
