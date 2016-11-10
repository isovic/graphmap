/*
 * owler_result.h
 *
 *  Created on: Jun 6, 2016
 *      Author: isovic
 */

#ifndef SRC_OWLER2_OVERLAP_FACTORY_H_
#define SRC_OWLER2_OVERLAP_FACTORY_H_

#include "overlap.h"
#include "containers/clusters.h"
#include <stdint.h>
#include <vector>
#include <memory>

class OverlapFactory {
 public:
  static void GetOverlap(const unsigned __int128 *packed_hits, int64_t n_hits, const std::vector<ClusterAndIndices *> &clusters,
              const std::vector<int> &lcskpp_indices, const std::vector<int32_t> &cluster_ids,
              int64_t read_id, const std::string &read_name, int64_t read_len,
              int64_t ref_id, const std::string &ref_name, int64_t ref_len, bool is_ref_rev, OverlapLine &ret);
 private:
  OverlapFactory();
  ~OverlapFactory();
};

#endif /* SRC_OWLER2_OVERLAP_FACTORY_H_ */
