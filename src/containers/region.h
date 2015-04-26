/*
 * region.h
 *
 *  Created on: Dec 21, 2014
 *      Author: ivan
 */

#ifndef REGION_H_
#define REGION_H_

#include <stdlib.h>
#include <string>
#include <sstream>
#include "index/index.h"

struct Region {
  int64_t start = 0;
  int64_t end = 0;
  int64_t reference_id = -1;
  std::string rname;
  int64_t region_index = -1;
  int64_t region_votes = 0;
  bool is_split = false;
  int64_t split_start = 0;
  int64_t split_end = 0;
};

// If the region is split in two parts, that is if the genome is circular, this function copies both parts in a new data array.
// It is users responsibility to free the allocated space using delete[].
int ConcatenateSplitRegion(const Index *index_reference, const Region &region, int8_t **ret_concatenated_data, int64_t *ret_data_length, int64_t *ret_start_offset, int64_t *ret_position_of_ref_end);
// Simply verbose region's details.
std::string VerboseRegionAsString(Region &region);

#endif /* REGION_H_ */
