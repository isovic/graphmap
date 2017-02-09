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
// #include "index/index.h"
#include "minimizer_index/minimizer_index.h"

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

//// Creates a copy of the region data from the Index.
//int CopyLinearRegion(const MinimizerIndex *index_reference, const Region *region, int8_t **ret_concatenated_data, int64_t *ret_data_length, int64_t *ret_start_offset);

// If the region is split in two parts, that is if the genome is circular, this function copies both parts in a new data array.
// It is users responsibility to free the allocated space using delete[].
int ConcatenateSplitRegion(std::shared_ptr<is::MinimizerIndex> index_reference, const Region *region, int8_t **ret_concatenated_data, int64_t *ret_data_length, int64_t *ret_start_offset, int64_t *ret_position_of_ref_end);

// Checks if the region is linear or split. If the region is linear, it returns the pointer to the existing part of the Index data and is_cleanup_required is set to false.
// Otherwise, a new data array is allocated and the data copied from the split parts of the Index.
// If the is_cleanup_required parameter is true, region_data needs to be freed by the user using free().
int GetRegionData(std::shared_ptr<is::MinimizerIndex> index, const Region *region,
                  int8_t **region_data, int64_t *data_len, int64_t *index_pos, int64_t *index_pos_of_ref_end, bool *is_cleanup_required);

//// Checks if the region is linear or split. It copies the data to a new array, and returns the pointer to the region data.
//// region_data needs to be freed by the user using free().
//int GetRegionDataCopy(const MinimizerIndex *index, const Region *region,
//                  int8_t **region_data, int64_t *data_len, int64_t *index_pos, int64_t *index_pos_of_ref_end);

// Simply verbose region's details.
std::string VerboseRegionAsString(Region &region);

#endif /* REGION_H_ */
