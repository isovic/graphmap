/*
 * clusters.h
 *
 *  Created on: Jul 25, 2016
 *      Author: isovic
 */

#ifndef SRC_CONTAINERS_CLUSTERS_H_
#define SRC_CONTAINERS_CLUSTERS_H_

#include <stdint.h>
#include <vector>
#include "containers/range.h"

struct ClusterAndIndices {
  Range query;
  Range ref;
  int32_t num_anchors = 0;
  int32_t coverage = 0;
  std::vector<int> lcskpp_indices;
};

#endif /* SRC_CONTAINERS_CLUSTERS_H_ */
