/*
 * mapping_data.cc
 *
 *  Created on: Mar 19, 2015
 *      Author: isovic
 */

#include "mapping_data.h"

MappingData::MappingData() {
  bins.clear();
  intermediate_mappings.clear();
  final_mapping_ptrs.clear();

  bin_size = -1;

  num_seeds_over_limit = 0;
  num_seeds_with_no_hits = 0;
  num_seeds_errors = 0;

  num_similar_mappings = 0;
  num_same_mappings = 0;
  avg_covered_bases_of_all_mappings = 0;
  std_covered_bases_of_all_mappings = 0;
  median_covered_bases_of_all_mappings = 0;

  iteration = 0;

  unmapped_reason = std::string("");

  num_region_iterations = 0;
  mapping_quality = 0;
  metagen_alignment_score = 0;

  time_region_selection = 0.0;
  time_mapping = 0.0;
  time_alignment = 0.0;
  time_region_seed_lookup = 0.0;
  time_region_hitsort = 0.0;
  time_region_conversion = 0.0;
  time_region_alloc = 0.0;
  time_region_counting = 0.0;
}

MappingData::~MappingData() {
  clear();
}

bool MappingData::IsMapped() {
  for (int32_t i=0; i<final_mapping_ptrs.size(); i++) {
    if (final_mapping_ptrs[i]->IsMapped() == true) { return true; };
  }
  return false;
}

bool MappingData::IsAligned() {
  for (int32_t i=0; i<final_mapping_ptrs.size(); i++) {
    if (final_mapping_ptrs[i]->IsAligned() == true) { return true; };
  }
  return false;
}

std::string MappingData::VerboseMappingDataToString_(const std::vector<PathGraphEntry *> *mapping_data, std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read) const {
  std::stringstream ss;

  int64_t reference_length = index->get_data().size();
  int64_t read_length = read->get_data_length();

  ss << "-----------------------\n";
  ss << "--- num_entries = " << mapping_data->size() << "\n";
  ss << "--- read id = " << read->get_sequence_absolute_id() << "\n";
  ss << "--- read name = " << read->get_header() << "\n";
  ss << "--- read_length = " << read_length << "\n";
  ss << "--- reference_length = " << reference_length << "\n";

  for (int64_t i = (mapping_data->size() - 1); i >= 0; i--) {
//    ss << "--- [" << i << "] ";
    ss << "[" << i << "/" << mapping_data->size() << "] ";
    int64_t start_location = 0, start_location_raw = 0;

    ss << "local_score_id = " << mapping_data->at(i)->get_mapping_data().local_score_id;
    ss << "\n      ° " << mapping_data->at(i)->VerboseToString();
    ss << "\n      ° r_id = " << mapping_data->at(i)->get_region_data().reference_id << ", fwd_r_id = " << (mapping_data->at(i)->get_region_data().reference_id % index->get_num_sequences_forward()) << ", region_index = " << mapping_data->at(i)->get_region_data().region_index;
    ss << "\n        ° \"" << index->get_headers()[mapping_data->at(i)->get_region_data().reference_id % index->get_num_sequences_forward()] << "\"";
    ss << "\n      ° Unmapped reason: \"" << unmapped_reason << "\"";
    int64_t relative_position = 0;
    int64_t absolute_position = 0;
    SeqOrientation orientation = kForward;

    ///// TODO: 06.02.2017.
//     This chunk below was removed due to the incompatibilities with the new index.
//        int64_t reference_id = index->RawPositionConverter(start_location, 0, &absolute_position, &relative_position, &orientation);
    int64_t reference_id = mapping_data->at(i)->get_region_data().reference_id;

        int64_t reference_start = mapping_data->at(i)->get_mapping_data().ref_coords.start;
//        index->RawPositionConverter(mapping_data->at(i)->get_mapping_data().ref_coords.start, 0, &absolute_position, &reference_start, &orientation);
        int64_t reference_end = mapping_data->at(i)->get_mapping_data().ref_coords.end;
//        index->RawPositionConverter(mapping_data->at(i)->get_mapping_data().ref_coords.end, 0, &absolute_position, &reference_end, &orientation);

        for (int64_t j = 0; j < mapping_data->at(i)->get_alignments().size(); j++) {
          ss << "\n      ° Alignment " << j << " / " << mapping_data->at(i)->get_alignments().size();
          ss << "\n        ° r_id = " << mapping_data->at(i)->get_region_data().reference_id << ", region_index = " << mapping_data->at(i)->get_region_data().region_index << ", region_votes = " << mapping_data->at(i)->get_region_data().region_votes << ", position = " << relative_position << ", r1[" << reference_start << ", " << reference_end << "], " << ((orientation == kForward) ? "forward" : "reverse");
          ss << ", sam_NM = " << mapping_data->at(i)->get_alignments()[j].edit_distance << ", sam_AS = " << mapping_data->at(i)->get_alignments()[j].alignment_score << ", sam_evalue = " << mapping_data->at(i)->get_alignments()[j].evalue << ", sam_pos = " << mapping_data->at(i)->get_alignments()[j].ref_start << ", sam_mapq = " << ((int64_t) mapping_data->at(i)->get_alignments()[j].mapping_quality) << ", relative_position = " << relative_position;
          ss << "\n        ° r_len = " << index->get_reference_lengths()[mapping_data->at(i)->get_region_data().reference_id] << ", l1_l = " << mapping_data->at(i)->get_l1_data().l1_l <<
              ", match_rate = " << ((float) mapping_data->at(i)->get_alignments()[j].num_eq_ops) / ((float) mapping_data->at(i)->get_alignments()[j].nonclipped_length) <<
              ", error_rate = " << ((float) mapping_data->at(i)->get_alignments()[j].num_x_ops + mapping_data->at(i)->get_alignments()[j].num_d_ops + mapping_data->at(i)->get_alignments()[j].num_i_ops) / ((float) mapping_data->at(i)->get_alignments()[j].nonclipped_length) <<
              " (X: = " << ((float) mapping_data->at(i)->get_alignments()[j].num_x_ops) / ((float) mapping_data->at(i)->get_alignments()[j].nonclipped_length) <<
              ", I = " << ((float) mapping_data->at(i)->get_alignments()[j].num_i_ops) / ((float) mapping_data->at(i)->get_alignments()[j].nonclipped_length) <<
              ", D: = " << ((float) mapping_data->at(i)->get_alignments()[j].num_d_ops) / ((float) mapping_data->at(i)->get_alignments()[j].nonclipped_length) << ")";

          ss << "\n        ° \"" << index->get_headers()[mapping_data->at(i)->get_region_data().reference_id % index->get_num_sequences_forward()] << "\"";
        }
    ss << "\n-----------";
    if (i == 0) {
      ss << "\n";
    }
    ss << "\n";
  }

  return ss.str();
}

std::string MappingData::VerboseFinalMappingsToString(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read) const {
  return VerboseMappingDataToString_(&final_mapping_ptrs, index, read);
}

std::string MappingData::VerboseIntermediateMappingsToString(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read) const {
  return VerboseMappingDataToString_(&intermediate_mappings, index, read);
}

void MappingData::clear() {
  vertices.Clear();
  bins.clear();
  for (int64_t i = 0; i < intermediate_mappings.size(); i++) {
    if (intermediate_mappings[i])
      delete intermediate_mappings[i];
    intermediate_mappings[i] = NULL;
  }
  intermediate_mappings.clear();
  unmapped_reason = std::string("");
  num_region_iterations = 0;
  mapping_quality = 0;
  metagen_alignment_score = 0;
}
