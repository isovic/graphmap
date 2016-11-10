/*
 * owler_result.cc
 *
 *  Created on: Jun 6, 2016
 *      Author: isovic
 */

#include <owler2/overlap_factory.h>
#include "owler2/lcskpp_packed.h"

typedef unsigned __int128 uint128_t;

OverlapFactory::OverlapFactory() {

}

OverlapFactory::~OverlapFactory() {

}

void OverlapFactory::GetOverlap(
    const uint128_t* packed_hits, int64_t n_hits,
    const std::vector<ClusterAndIndices*>& clusters,
    const std::vector<int>& lcskpp_indices,
    const std::vector<int32_t>& cluster_ids,
    int64_t read_id, const std::string &read_name, int64_t read_len,
    int64_t ref_id, const std::string &ref_name, int64_t ref_len, bool is_ref_rev, OverlapLine &ret) {

  if (lcskpp_indices.size() == 0) {
    return;
  }
//  int32_t Aid = get_lcsk128_qid(packed_hits[lcskpp_indices.front()]);
//  int32_t Bid = get_lcsk128_rid(packed_hits[lcskpp_indices.front()]);
  int32_t Astart = get_lcsk128_qpos(packed_hits[lcskpp_indices.front()]) & 0x0FFFFFFFF;
  int32_t Bstart = get_lcsk128_rpos(packed_hits[lcskpp_indices.front()]) & 0x0FFFFFFFF;
  int32_t Aend = get_lcsk128_qpos(packed_hits[lcskpp_indices.back()]) & 0x0FFFFFFFF;
  int32_t Bend = get_lcsk128_rpos(packed_hits[lcskpp_indices.back()]) & 0x0FFFFFFFF;

  int64_t shared_minmers = lcskpp_indices.size();
  double perc_err = std::min(((double) shared_minmers * 12) / ((double) (Aend - Astart)), 1.0);

  if (is_ref_rev) {
    int64_t temp = Bstart;
    Bstart = ref_len - Bend - 1;
    Bend = ref_len - temp - 1;
  }

  ret = OverlapLine(read_id, ref_id, read_name, ref_name, perc_err, shared_minmers,
                                                   0, Astart, Aend, read_len,
                                                   (is_ref_rev == false) ? 0 : 1, Bstart, Bend, ref_len);

//  return ret;
}
