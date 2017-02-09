/*
 * region.cc
 *
 *  Created on: Dec 26, 2014
 *      Author: isovic
 */

#include "containers/region.h"

//int CopyLinearRegion(const Index *index_reference, const Region *region, int8_t **ret_concatenated_data, int64_t *ret_data_length, int64_t *ret_start_offset) {
//  if (region->is_split == true)
//    return 1;
//
//  int8_t *data_copy = new int8_t[(region->end - region->start + 1) + 1];
//  if (data_copy == NULL) {
//    return 3;
//  }
//
//  memmove((data_copy), &(index_reference->get_data()[region->start]), (region->end - region->start + 1));
//
//  data_copy[(region->end - region->start + 1)] = '\0';
//
//  *ret_concatenated_data = data_copy;
//  *ret_data_length = (region->end - region->start + 1);
//  *ret_start_offset = region->start;
//
//  return 0;
//}

int ConcatenateSplitRegion(std::shared_ptr<is::MinimizerIndex> index_reference, const Region *region, int8_t **ret_concatenated_data, int64_t *ret_data_length, int64_t *ret_start_offset, int64_t *ret_position_of_ref_end) {
  if (region->is_split == false)
    return 1;

  int64_t region_length_first = (region->end - region->start + 1);
  int64_t region_length_second = (region->split_end - region->split_start + 1);
  int64_t region_length_joined = region_length_first + region_length_second;
  if (region_length_first <= 0 || region_length_second <= 0 || region_length_joined <= 0)
    return 2;

  int8_t *data_copy = new int8_t[region_length_joined + 1];
  if (data_copy == NULL) {
    return 3;
  }

  int64_t start_offset = 0;
  int64_t position_of_ref_end = 0;

  // If the main region is at the beginning of the reference. The region is then expanded towards left and right, but on the left it zips back
  // to the end of the circular reference.
  if (region->start < region->split_start) {
    memmove(data_copy, &(index_reference->get_data()[region->split_start]), region_length_second);
    memmove((data_copy + region_length_second), &(index_reference->get_data()[region->start]), region_length_first);
    position_of_ref_end = region->split_end - region->split_start; // + 1;
    start_offset = region->split_start;

    // If the main region is at the end of the reference. The region is then expanded towards left and right, but on the right it zips back
    // to the beginning of the circular reference.
  } else {
    memmove((data_copy), &(index_reference->get_data()[region->start]), region_length_first);
    memmove((data_copy + region_length_first), &(index_reference->get_data()[region->split_start]), region_length_second);
    position_of_ref_end = region->end - region->start;
    start_offset = region->start;

  }

  data_copy[region_length_joined] = '\0';

  *ret_concatenated_data = data_copy;
  *ret_data_length = region_length_joined;
  *ret_start_offset = start_offset;
  *ret_position_of_ref_end = position_of_ref_end;

  return 0;
}

int GetRegionData(std::shared_ptr<is::MinimizerIndex> index, const Region *region,
                  int8_t **region_data, int64_t *data_len, int64_t *index_reg_start, int64_t *pos_of_ref_end, bool *is_cleanup_required) {

  if (region->is_split == false) {
    *region_data = (int8_t *) (&index->get_data()[0] + region->start);
    *data_len = (region->end - region->start);
    *index_reg_start = region->start;
    *pos_of_ref_end = -1;
    *is_cleanup_required = false;

  } else {
    ConcatenateSplitRegion(index, region, region_data, data_len, index_reg_start, pos_of_ref_end);
    *is_cleanup_required = true;

  }

  return 0;
}

//int GetRegionDataCopy(const Index *index, const Region *region,
//                  int8_t **region_data, int64_t *data_len, int64_t *index_reg_pos, int64_t *reg_pos_of_ref_end) {
//
//  if (region->is_split == false) {
//    CopyLinearRegion(index, region, region_data, data_len, index_reg_pos);
//    *reg_pos_of_ref_end = -1;
//
//  } else {
//    ConcatenateSplitRegion(index, region, region_data, data_len, index_reg_pos, reg_pos_of_ref_end);
//
//  }
//
//  return 0;
//}

std::string VerboseRegionAsString(Region &region) {
  std::stringstream ss;

  ss << "start = " << region.start;
  ss << ", end = " << region.end;
  ss << ", reference_id = " << region.reference_id;
  ss << ", region_index = " << region.region_index;
  ss << ", region_votes = " << region.region_votes;
  ss << ", is_split = " << ((int) region.is_split);
  ss << ", split_start = " << region.split_start;
  ss << ", split_end = " << region.split_end;

  return ss.str();
}
