/*
 * lcskpp.h
 *
 *  Created on: Nov 10, 2015
 *      Author: isovic
 */

#ifndef LCSKPP_H_
#define LCSKPP_H_

#include <stdlib.h>
#include <stdint.h>
#include <algorithm>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cassert>
#include <cmath>
#include <cstdlib>

#include "algorithm/fenwick.h"
#include "sequences/single_sequence.h"
#include "program_parameters.h"
#include "utility/utility_general.h"
#include "log_system/log_system.h"
#include "containers/clusters.h"

using int128_t = __int128;
using uint128_t = unsigned __int128;

#define get_lcsk128_qpos(x)  ((int32_t) (((x) >> 0) & 0x0FFFFFFFF))
#define get_lcsk128_qid(x)   ((int32_t) (((x) >> 32) & 0x0FFFFFFFF))
#define get_lcsk128_rpos(x)  ((int32_t) (((x) >> 64) & 0x0FFFFFFFF))
#define get_lcsk128_rid(x)   ((int32_t) (((x) >> 96) & 0x0FFFFFFFF))
///      d                  c                 b             a
/// ref_id << 96 | query_start << 64 | ref_start << 32 | query_id
#define pack_lcsk128(qstart,rstart,qid,rid)         ((((uint128_t) rid) << 96) | (((uint128_t) rstart) << 64) | (((uint128_t) qid) << 32) | ((uint128_t) qstart))

/// packed_hits needs to be sorted.
int LCSkPacked(const uint128_t *packed_hits, int64_t n_hits, int32_t k, int64_t *ret_lcskpp_length, std::vector<int32_t> *ret_lcskpp_indices, bool print_debug=false);

int FilterAnchorsByChainingPacked(const uint128_t *packed_hits, int64_t n_hits, int32_t k, const SingleSequence* read, const ProgramParameters *parameters,
                  const std::vector<int> &lcskpp_indices, double indel_bandwidth_margin, int32_t max_dist, int32_t lookahead_dist_factor, int64_t min_covered_bases, int32_t cluster_size_cutoff,
                  std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids);

int GenerateClustersPacked(const uint128_t *packed_hits, int64_t n_hits, int32_t k, int64_t min_num_anchors_in_cluster, int64_t min_cluster_length, int64_t min_cluster_covered_bases,
                     float min_cluster_coverage, std::vector<int> &lcskpp_indices,
                     const SingleSequence* read, const ProgramParameters* parameters, std::vector<ClusterAndIndices *> &ret_clusters,
                     std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids);

int WriteLCSkDebug(std::string out_file, std::string qname, int64_t qlen, std::string rname, int64_t rlen, int32_t k, const uint128_t *packed_hits, int64_t n_hits, std::vector<int32_t> *lcskpp_indices, std::vector<int32_t> *cluster_ids);

void TestLCSk();
void TestLCSk2();
void TestLCSk3();

#endif /* LCSKPP_H_ */
