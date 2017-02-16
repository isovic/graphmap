/*
 * owler_data.h
 *
 *  Created on: Jul 2, 2015
 *      Author: isovic
 */

#ifndef OWLER_DATA_H_
#define OWLER_DATA_H_

#include <memory>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <sstream>
#include <vector>

#include "minimizer_index/minimizer_index.h"
#include "containers/range.h"

class PairwiseOverlap {
 public:
  PairwiseOverlap() : qid(0), tid(0), num_seeds(0), cov_bases(0) { }

  Range query, target;
  int64_t qid, tid;
  int64_t num_seeds;
  int64_t cov_bases;
};

class OwlerData {
 public:
  OwlerData() { };
  ~OwlerData() { };

  std::vector<uint128_t> hits;
  std::vector<PairwiseOverlap> overlaps;
  std::string unmapped_reason;
  std::vector<std::string> out_lines;
};

#endif /* OWLER_DATA_H_ */
