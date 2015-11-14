/*
 * process_read.cc
 *
 *  Created on: Mar 20, 2015
 *      Author: isovic
 */

#include <limits>
#include <algorithm>
#include "owler/owler.h"
//#include "sam/sam_entry.h"
#include "index/index_hash.h"

#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include "utility/fenwick.h"



int Owler::ProcessRead(OwlerData *owler_data, std::vector<Index *> indexes, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params) {
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("\n"), "[]");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("Entered function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024)), "ProcessRead");

  // If the read length is too short, call it unmapped.
  if (read->get_sequence_length() < 80) {
    std::stringstream ss;
    ss << "Unmapped_5__readlength_too_short" << "__readlength=" << read->get_sequence_length() << "__limit=" << 80;
    owler_data->unmapped_reason += ss.str();
    return 0;
  }

//  IndexSpacedHash test_index;
//  test_index.GenerateFromSingleSequence(*read);

  owler_data->seed_hits2.reserve(500000);
  /// Check if it's a case of self-overlap. In this case, overlap can be performed faster, because the index will already have pre-processed seeds of all reads.
  if (parameters->reads_path == parameters->reference_path)
    CollectSeedHitsExperimentalSubseededIndex(owler_data, indexes, read, parameters);
  else
    CollectSeedHitsExperimentalCalcSubseedsFast(owler_data, indexes, read, parameters);



//  CollectSeedHits(owler_data, indexes, read, parameters);
//  FilterUnlikelyOverlaps(owler_data, indexes, read, parameters);
//  ApplyLCS(owler_data, indexes, read, parameters);

  ApplyLCS2(owler_data, indexes, read, parameters);

//  std::string overlap_lines;
//  CollectMHAPLines(std::string &overlap_lines, owler_data, read, parameters)

//  int Owler::CollectAMOSLines(std::string &ret_amos_lines, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters)



//  float threshold_step = 0.10f;
//  float bin_value_threshold = mapping_data->bins.front().bin_value;
//  bin_value_threshold = (mapping_data->bins.front().bin_value * (1.0f - threshold_step));
//  bin_value_threshold = std::max(bin_value_threshold, 2.0f);
//  int64_t num_regions_within_threshold = CountBinsWithinThreshold_(mapping_data, bin_value_threshold);
//  float min_allowed_bin_value = 0.0f;
//  if (parameters->parsimonious_mode) {
//    min_allowed_bin_value = (0.50f * mapping_data->bins.front().bin_value);
//
//  } else {
//    min_allowed_bin_value = (0.75f * mapping_data->bins.front().bin_value);
////    min_allowed_bin_value = 0.0f;
//  }
//
//  int64_t max_num_regions = parameters->max_num_regions;
//  if (parameters->alignment_approach == "overlapper") {
//    max_num_regions = mapping_data->bins.size();
////    int64_t smaller_seq_length = std::min(read->get_sequence_length(), index->get_reference_lengths()[0]);
////    min_allowed_bin_value = std::max((0.10f * mapping_data->bins.front().bin_value), 0.10f * read->get_sequence_length());
////    min_allowed_bin_value = 0.10f * read->get_sequence_length();
////    min_allowed_bin_value = 0.01f * read->get_sequence_length();
//    min_allowed_bin_value = std::min(0.10f * read->get_sequence_length(), 100.0f);
//  }
//
//  mapping_data->num_region_iterations = 0;
//
//
//  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Top 10 scoring bins:\n"), "ProcessRead");
//  for (int64_t i = 0; i < mapping_data->bins.size() && i < 10; i++) {
//    Region region = CalcRegionFromBin_(i, mapping_data, read, parameters);
//    ScoreRegistry local_score(region, i);
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[i = %ld] location_start = %ld, location_end = %ld, is_reverse = %d, vote = %ld, region_index = %ld\n", i, region.start, region.end, (int) (region.start >= index_->get_data_length_forward()), region.region_votes, region.region_index), "ProcessRead");
//  }
//
//  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n\n", "ProcessRead");
//
//  // Process regions one by one.
//  for (int64_t i = 0; i < mapping_data->bins.size() && i <max_num_regions; i++) {
////    if (parameters->alignment_algorithm == "overlapper") {
////      if (index->get_headers()[mapping_data->bins[i].reference_id % index->get_num_sequences_forward()] == ((std::string) read->get_header())) {
////        continue;
////      }
////    }
//
////  for (int64_t i = 0; i < mapping_data->bins.size(); i++) {
//    // If the ret_check value is zero, then just continue as normal.
//    int ret_check = 0;
//
//    if (parameters->alignment_approach != "overlapper") {
//      ret_check = CheckRegionSearchFinished_(i, min_allowed_bin_value, threshold_step, &bin_value_threshold, mapping_data, read, parameters);
//    } else {
//      if (mapping_data->bins[i].bin_value < 1 || mapping_data->bins[i].bin_value < min_allowed_bin_value)
//        break;
//    }
//
//    // Region search needs to stop.
//    if (ret_check < 0) {
//
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("CheckRegionSearchFinished returned with value to break! ret_check = %d\n", ret_check), "ProcessRead");
//
//      break;
//
//      // Another iteration needs to be performed.
//    } else if (ret_check > 0) {
//      mapping_data->num_region_iterations += 1;
//
////      if (mapped_data.num_region_iterations > 1 && (num_similar_mappings_ - num_same_mappings_) > 10 && ret_check == 5) {
////        std::stringstream ss;
////        ss << "Unmapped_15_mapping_quality_is_zero_after_batch_completed." << "__readlength=" << read->get_sequence_length() << "__max_region_votes=" << mapped_data.bins.front().bin_value << "__num_region_iterations=" << mapped_data.num_region_iterations;
////        mapped_data.unmapped_reason += ss.str();
////        return 0;
////      }
//    }
//  // for (int64_t i = 0; i < mapping_data->bins.size(); i++) {
//
//  //   if (mapping_data->bins[i].bin_value < min_allowed_bin_value)
//  //     break;
//  //   if (i >= parameters->max_num_regions && ((i > 0 && mapping_data->bins[i].bin_value != mapping_data->bins[i-1].bin_value) || i == 0))
//  //     break;
//
//    Region region = CalcRegionFromBin_(i, mapping_data, read, parameters);
//    ScoreRegistry local_score(region, i);
//
//    bool is_reverse = (region.start >= index_->get_data_length_forward());
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[i = %ld] location_start = %ld, location_end = %ld, is_reverse = %d, vote = %ld, region_index = %ld\n", i, region.start, region.end, (int) (region.start >= index_->get_data_length_forward()), region.region_votes, region.region_index), "ProcessRead");
//
//    // Perform the GraphMap on a single region.
//    GraphMap_(&local_score, index_read, mapping_data, index_primary, index_secondary, read, parameters);
//
//    // Just verbose.
//    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Local scores (raw, before LCSk):\n"), "ProcessRead");
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("%s", local_score.VerboseToString().c_str()), "ProcessRead");
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Running PostProcessRegionWithLCS_. j = %ld / %ld, local_score.size() = %ld\n", i, mapping_data->bins.size(), local_score.get_registry_entries().num_vertices), "ProcessRead");
//    }
//
////    #ifndef RELEASE_VERSION
//    if (parameters->alignment_algorithm == "myers" || parameters->alignment_algorithm == "gotoh") {
//      int ret_value_lcs = PostProcessRegionWithLCS_(&local_score, mapping_data, index_primary, index_secondary, read, parameters);
//    } else {
//      int ret_value_lcs = ExperimentalPostProcessRegionWithLCS_(&local_score, mapping_data, index_primary, index_secondary, read, parameters);
//    }
//
////    if (parameters->alignment_algorithm == "anchor") {
////      int ret_value_lcs = ExperimentalPostProcessRegionWithLCS_(&local_score, mapping_data, index, index_secondary, read, parameters);
////    } else {
////      int ret_value_lcs = PostProcessRegionWithLCS_(&local_score, mapping_data, index, index_secondary, read, parameters);
////    }
//
////    #else
////        int ret_value_lcs = PostProcessRegionWithLCS_(&local_score, mapping_data, index, index_secondary, read, parameters);
////    #endif
//    local_score.Clear();
//
//    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("-----\n\n"), "ProcessRead");
//    }
//
////    CheckMinimumMappingConditions(mapping_data, parameters);
//  }
//
//  if (index_read)
//    delete index_read;
//  index_read = NULL;
//  mapping_data->vertices.Clear();
//
//  GenerateAlignments_(mapping_data, index_primary, read, parameters, evalue_params);

  // Just verbose.
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n"), "[]");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, ((parameters->num_threads == 1) || read->get_sequence_id() == parameters->debug_read), FormatString("Exiting function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024)), "ProcessRead");
  }

  return 0;
}


int Owler::CollectSeedHits(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  if (owler_data == NULL || read == NULL || parameters == NULL)
    return 1;
  if (indexes.size() == 0 || (indexes.size() > 0 && indexes[0] == NULL))
    return 2;

  int64_t readlength = read->get_sequence_length();

  /// Initialize the data structures to hold the results.
  owler_data->Init((SingleSequence*) read, indexes);

  for (int64_t i = 0; i < readlength; i += parameters->kmer_step) {
    int8_t *seed = (int8_t *) &(read->get_data()[i]);

    int32_t seed_shape_id = 0;

    /// Loop through all indexes to collect all different hits.
    for (int64_t j = 0; j < indexes.size(); j++) {
      seed_shape_id = j;

      Index *index = indexes[j];
      if (index == NULL)
        continue;

      int64_t k = (int64_t) ((IndexSpacedHash *) index)->get_shape_index_length();
      if ((i + k) >= readlength)
        continue;

      uint64_t hits_start = 0, num_hits = 0;
      int64_t *hits = NULL;

      int ret_search = index->FindAllRawPositionsOfSeed(seed, k, parameters->max_num_hits, &hits, &hits_start, &num_hits);

      // Check if there is too many hits (or too few).
      if (ret_search == 1) {
        owler_data->num_seeds_with_no_hits += 1;
      } else if (ret_search == 2) {
        owler_data->num_seeds_over_limit += 1;
      } else if (ret_search > 2) {
        owler_data->num_seeds_errors += 1;
      }

      /// Counting kmers in regions of bin_size on the genome
      for (int64_t j1 = hits_start; j1 < (hits_start + num_hits); j1++) {
        int64_t position = hits[j1];

        /// Find the index of the reference that was hit. This also includes the reverse sequences.
        /// Reverse sequences are considered the same as any other reference sequence.
        int64_t reference_index = index->RawPositionToReferenceIndexWithReverse(position);
        int64_t reference_length = index->get_reference_lengths()[reference_index];
        int64_t reference_start = index->get_reference_starting_pos()[reference_index];
        int64_t reference_end = reference_start + reference_length;
        int64_t position_local = position - reference_start;

        if (reference_index < 0) {
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index.\n"), "SelectRegionsWithHoughAndCircular");
          continue;
        }
        /// Don't count self hits
        if (index->get_headers()[reference_index % index->get_num_sequences_forward()] == ((std::string) read->get_header())) {
          continue;
        }
        /// Count unique hits for a pair of reads.
        if (owler_data->overlaps[reference_index].last_update < (i + 1)) {
          owler_data->overlaps[reference_index].num_unique_hits += 1;
        }

        owler_data->overlaps[reference_index].last_update = (i + 1);
        SeedHit seed_hit;
        owler_data->overlaps[reference_index].seed_hits.push_back(SeedHit((uint32_t) i, (uint32_t) position_local, seed_shape_id));
      }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)

      if (index->IsManualCleanupRequired("FindAllRawPositionsOfSeed") == 0) {
        if (hits)
          free(hits);
        hits = NULL;
      }
      hits = NULL;
    }
  }

  return 0;
}

int Owler::FilterUnlikelyOverlaps(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  std::sort(owler_data->overlaps.begin(), owler_data->overlaps.end(), overlaps_greater_than_key());
  int64_t min_num_hits = std::min(100.0f, 0.10f * read->get_sequence_length());
//  for (int64_t i = 0; i < owler_data->overlaps.size(); i++) {
//    printf ("owler_data->overlaps[%ld].num_unique_hits = %ld\n", i, owler_data->overlaps[i].num_unique_hits);
//  }

  int64_t num_hits_to_leave = owler_data->overlaps.size();
  for (int64_t i = 0; i < owler_data->overlaps.size(); i++) {
    if (owler_data->overlaps[i].num_unique_hits < min_num_hits) {
      num_hits_to_leave = i;
      break;
    }
  }
  if (num_hits_to_leave > 0 && num_hits_to_leave < owler_data->overlaps.size()) {
    owler_data->overlaps.resize(num_hits_to_leave);
  } else {
    owler_data->overlaps.clear();
  }



  /// Just verbose.
  for (int64_t i=0; i<owler_data->overlaps.size(); i++) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[%ld / %ld] %s\n", i, owler_data->overlaps.size(), owler_data->overlaps[i].VerboseToString().c_str()), "[]");
  }

#ifndef RELEASE_VERSION
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing gapped qgrams to file.\n"), "FilterUnlikelyOverlaps");
    for (int64_t i=0; i<owler_data->overlaps.size(); i++) {
      std::string out_path = FormatString("temp/overlaps/scores-%ld.csv", owler_data->overlaps[i].ref_id_);
      FILE *fp = fopen(out_path.c_str(), "w");
      if (fp != NULL) {
        fprintf (fp, "%s\t0\t%ld\t0\t%ld\t0.0\n", read->get_header(), read->get_sequence_length(), indexes[0]->get_reference_lengths()[owler_data->overlaps[i].ref_id_]);
        for (int64_t j=0; j<owler_data->overlaps[i].seed_hits.size(); j++) {
          uint32_t seed_len = owler_data->seed_types[owler_data->overlaps[i].seed_hits[j].seed_type].length();
          fprintf (fp, "%ld\t%ld\n", owler_data->overlaps[i].seed_hits[j].query_pos, owler_data->overlaps[i].seed_hits[j].ref_pos);
          fprintf (fp, "%ld\t%ld\n", owler_data->overlaps[i].seed_hits[j].query_pos + seed_len, owler_data->overlaps[i].seed_hits[j].ref_pos + seed_len);
        }
        fclose(fp);
      }
    }

//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing gapped qgrams to file.\n"), "FilterUnlikelyOverlaps");
//    for (int64_t i=0; i<owler_data->overlaps.size(); i++) {
//      std::string out_path = FormatString("temp/overlaps/LCS-%ld.csv", owler_data->overlaps[i].ref_id_);
//      FILE *fp = fopen(out_path.c_str(), "w");
//      if (fp != NULL) {
//        fprintf (fp, "%s\t0\t%ld\t0\t%ld\t0.0\n", read->get_header(), read->get_sequence_length(), indexes[0]->get_reference_lengths()[owler_data->overlaps[i].ref_id_]);
//        for (int64_t j=0; j<owler_data->overlaps[i].seed_hits.size(); j++) {
//          uint32_t seed_len = owler_data->seed_types[owler_data->overlaps[i].seed_hits[j].seed_type].length();
//          fprintf (fp, "%ld\t%ld\n", owler_data->overlaps[i].seed_hits[j].query_pos, owler_data->overlaps[i].seed_hits[j].ref_pos);
//          fprintf (fp, "%ld\t%ld\n", owler_data->overlaps[i].seed_hits[j].query_pos + seed_len, owler_data->overlaps[i].seed_hits[j].ref_pos + seed_len);
//        }
//        fclose(fp);
//      }
//    }

//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing gapped qgrams to file.\n"), "FilterUnlikelyOverlaps");
    for (int64_t i=0; i<owler_data->overlaps.size(); i++) {
      std::string out_path = FormatString("temp/overlaps/LCSL1-%ld.csv", owler_data->overlaps[i].ref_id_);
      FILE *fp = fopen(out_path.c_str(), "w");
      if (fp != NULL) {
        fprintf (fp, "%s\t0\t%ld\t0\t%ld\t0.0\n", read->get_header(), read->get_sequence_length(), indexes[0]->get_reference_lengths()[owler_data->overlaps[i].ref_id_]);
        for (int64_t j=0; j<owler_data->overlaps[i].seed_hits.size(); j++) {
          uint32_t seed_len = owler_data->seed_types[owler_data->overlaps[i].seed_hits[j].seed_type].length();
          fprintf (fp, "%ld\t%ld\n", owler_data->overlaps[i].seed_hits[j].query_pos, owler_data->overlaps[i].seed_hits[j].ref_pos);
          fprintf (fp, "%ld\t%ld\n", owler_data->overlaps[i].seed_hits[j].query_pos + seed_len, owler_data->overlaps[i].seed_hits[j].ref_pos + seed_len);
        }
        fclose(fp);
      }
    }

//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing gapped qgrams to file.\n"), "FilterUnlikelyOverlaps");
    for (int64_t i=0; i<owler_data->overlaps.size(); i++) {
      std::string out_path = FormatString("temp/overlaps/double_LCS-%ld.csv", owler_data->overlaps[i].ref_id_);
      FILE *fp = fopen(out_path.c_str(), "w");
        if (fp != NULL) {
        fprintf (fp, "%s\t0\t%ld\t0\t%ld\t0.0\n", read->get_header(), read->get_sequence_length(), indexes[0]->get_reference_lengths()[owler_data->overlaps[i].ref_id_]);
        for (int64_t j=0; j<owler_data->overlaps[i].seed_hits.size(); j++) {
          uint32_t seed_len = owler_data->seed_types[owler_data->overlaps[i].seed_hits[j].seed_type].length();
          fprintf (fp, "%ld\t%ld\n", owler_data->overlaps[i].seed_hits[j].query_pos, owler_data->overlaps[i].seed_hits[j].ref_pos);
          fprintf (fp, "%ld\t%ld\n", owler_data->overlaps[i].seed_hits[j].query_pos + seed_len, owler_data->overlaps[i].seed_hits[j].ref_pos + seed_len);
        }
        fclose(fp);
      }
    }

  }
#endif

  return 0;
}

int Owler::ApplyLCS(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  for (int64_t i = 0; i < owler_data->overlaps.size(); i++) {
    int32_t lcskpp_length;
    std::vector<int> lcskpp_indices;
    CalcLCSFromLocalScoresCacheFriendly_(owler_data, i, &lcskpp_length, &lcskpp_indices);

#ifndef RELEASE_VERSION
    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    std::vector<SeedHit> filtered_seed_hits;
    filtered_seed_hits.resize(lcskpp_indices.size());
    for (int64_t j = 0; j < lcskpp_indices.size(); j++) {
      filtered_seed_hits[j] = owler_data->overlaps[i].seed_hits[lcskpp_indices[j]];
    }
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing gapped qgrams to file.\n"), "FilterUnlikelyOverlaps");
      std::string out_path = FormatString("temp/overlaps/LCS-%ld.csv", owler_data->overlaps[i].ref_id_);
      FILE *fp = fopen(out_path.c_str(), "w");
      if (fp != NULL) {
        fprintf (fp, "%s\t0\t%ld\t0\t%ld\t0.0\n", read->get_header(), read->get_sequence_length(), indexes[0]->get_reference_lengths()[owler_data->overlaps[i].ref_id_]);
        for (int64_t j=0; j<filtered_seed_hits.size(); j++) {
          uint32_t seed_len = owler_data->seed_types[filtered_seed_hits[j].seed_type].length();
          fprintf (fp, "%ld\t%ld\n", filtered_seed_hits[j].query_pos, filtered_seed_hits[j].ref_pos);
          fprintf (fp, "%ld\t%ld\n", filtered_seed_hits[j].query_pos + seed_len, filtered_seed_hits[j].ref_pos + seed_len);
        }
        fclose(fp);
      }
    }
#endif

  }

  return 0;
}

int Owler::CollectMHAPLines(std::string &ret_overlap_lines, OwlerData *owler_data, const SingleSequence *read, const ProgramParameters *parameters) {
  return STATE_MAPPED;
}



int Owler::CollectAMOSLines(std::string &ret_amos_lines, OwlerData *owler_data, const SingleSequence *read, const ProgramParameters *parameters) {
  std::stringstream ss;

//  int64_t num_mapped_alignments = 0;
//  int64_t num_unmapped_alignments = 0;
//  for (int64_t i = 0; i < mapping_data->final_mapping_ptrs.size(); i++) {
//    if (mapping_data->final_mapping_ptrs.at(i)->get_alignment_primary().is_aligned == false) {
//      num_unmapped_alignments += 1;
//      continue;
//    }
//    if (ss.str().size() > 0)
//      ss << "\n";
//    ss << mapping_data->final_mapping_ptrs.at(i)->GenerateAMOS();
//    num_mapped_alignments += 1;
//  }
//
////  int64_t num_mapped_alignments = 0;
////  int64_t num_unmapped_alignments = 0;
////  for (int64_t i = 0; i < mapping_data->final_mapping_ptrs.size(); i++) {
////    if (mapping_data->final_mapping_ptrs.at(i)->get_mapping_data().is_mapped == false) {
////      num_unmapped_alignments += 1;
////      continue;
////    }
////    if (ss.str().size() > 0)
////      ss << "\n";
////    ss << mapping_data->final_mapping_ptrs.at(i)->GenerateAMOS();
////    num_mapped_alignments += 1;
////  }
//
//  // If all reported alignments have been declared unmapped for some reason, output one of them to be consistent
//  // and provide an alignment line for each read.
//  // Also, there will always be at leas one final_mapping_ptrs_ entry, because we have handled the size() == 0 case
//  // a little bit up.
//  if (num_unmapped_alignments == mapping_data->final_mapping_ptrs.size() || mapping_data->unmapped_reason.size() > 0) {
//    std::stringstream temp_ss;
//    temp_ss << "__num_region_iterations=" << mapping_data->num_region_iterations;
//    mapping_data->unmapped_reason += FormatString("__num_unmapped_alignments=%ld", num_unmapped_alignments) + temp_ss.str();
//  }
//
//  if (mapping_data->unmapped_reason.size() > 0) {
////    ret_sam_lines = GenerateUnmappedSamLine_(mapping_data, parameters->verbose_sam_output, read);
//    return STATE_UNMAPPED;
//  }
//
//  ret_amos_lines = ss.str();

  return STATE_MAPPED;
}



// Assumes vertices are sorted.
// If use_l1_filtering is true, then all vertices/anchors that have coordinates further than allowed_dist from the L1 line are filtered out.
// Otherwise, all vertices will be used.
// The L1 line is specified with k = 1 and l parameters (y = k*x + l).
void Owler::CalcLCSFromLocalScoresCacheFriendly_(OwlerData* owler_data, int64_t overlap_id, int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices) {
  uint32_t num_vertices = owler_data->overlaps[overlap_id].seed_hits.size();
  PairwiseOverlapData *overlaps = &(owler_data->overlaps[overlap_id]);

  if (num_vertices <= 0)
    return;

  unsigned __int128 *events = NULL;
  uint64_t *matches_starts = NULL;
  uint64_t *matches_dists_ref = NULL;
  uint64_t *matches_indices = NULL;

  uint32_t n = 0;
  int64_t num_matches = 0;

  int64_t num_events = 0;
  int64_t lcskpp_length = 0;

  uint32_t min_ref = overlaps->seed_hits[0].ref_pos, min_query = overlaps->seed_hits[0].query_pos;

  for (uint32_t i=0; i<num_vertices; i++) {
    min_ref = std::min(min_ref, overlaps->seed_hits[i].ref_pos);
    min_query = std::min(min_query, overlaps->seed_hits[i].query_pos);
  }

  events = (unsigned __int128 *) malloc(sizeof(unsigned __int128) * num_vertices * 2);
  matches_starts = (uint64_t *) malloc(sizeof(uint64_t) * num_vertices);
  matches_dists_ref = (uint64_t *) malloc(sizeof(uint64_t) * num_vertices);
  matches_indices = (uint64_t *) malloc(sizeof(uint64_t) * num_vertices);

  for (uint32_t i=0; i<num_vertices; i++) {
    uint32_t seed_len = owler_data->seed_types[overlaps->seed_hits[i].seed_type].length();
    uint32_t ref_start = overlaps->seed_hits[i].ref_pos - min_ref;
    uint32_t ref_end = ref_start + seed_len;
    uint32_t query_start = overlaps->seed_hits[i].query_pos - min_query;
    uint32_t query_end = query_start + seed_len;

    unsigned __int128 event1 = (((unsigned __int128) ref_start) << (8 * 8)) | (((unsigned __int128) query_start) << (4 * 8)) | (((unsigned __int128) (i + num_vertices)));
    events[num_events] = event1;
    num_events += 1;

    unsigned __int128 event2 = (((unsigned __int128) ref_end) << (8 * 8)) | (((unsigned __int128) query_end) << (4 * 8)) | ((((unsigned __int128) i)));
    events[num_events] = event2;
    num_events += 1;

    matches_starts[num_matches] = (uint64_t) (event1 >> (4 * 8));
    matches_dists_ref[num_matches] = seed_len;
    matches_indices[num_matches] = i;

    num_matches += 1;

    n = std::max(n, ref_end);
    n = std::max(n, query_end);
  }

  std::sort(events, (events + num_events));

  // Indexed by column, first:dp value, second:index in matches.
  FenwickMax<std::pair<int, int> > dp_col_max(n);
  std::vector<int> dp(num_matches);
  std::vector<int> recon(num_matches);
  std::vector<int> continues(num_matches, -1);

//    matches bi trebao biti sortiran na pocetku!

  for (int64_t curr = 0; curr < num_matches; curr++) {
    uint64_t G = 0;
    uint64_t G1 = (((matches_starts[curr] >> (4 * 8)) & (0x00000000FFFFFFFF)) - 1);
    uint64_t G2 = (((matches_starts[curr]) & (0x00000000FFFFFFFF)) - 1);
    G = (G1 << (4 * 8)) | G2;

    auto prev = std::lower_bound(matches_starts, (matches_starts + num_matches), G);

    if (prev != (matches_starts + num_matches) && *prev == G) {
      continues[curr] = prev - matches_starts;
    }
  }

  int best_idx = 0;
  lcskpp_length = 0;

  for (int64_t current_event=0; current_event<num_events; current_event++) {
    int64_t raw_idx = (int64_t) ((events[current_event] & (0x00000000FFFFFFFF)));

    int64_t idx = (raw_idx >= ((int64_t) num_matches)) ? (raw_idx - ((int64_t) num_matches)) : (raw_idx);
    bool is_beginning = (raw_idx >= num_matches);
    uint64_t i = (uint64_t) ((events[current_event] >> (8 * 8)) & (0x00000000FFFFFFFF));
    uint64_t j = (uint64_t) ((events[current_event] >> (4 * 8)) & (0x00000000FFFFFFFF));
    int primary_diagonal = n - 1 + i - j;

    if (is_beginning) { // begin
      std::pair<int, int> prev_dp = dp_col_max.get(j);
      uint64_t k_length = matches_dists_ref[idx];
      dp[idx] = k_length;      // k
      recon[idx] = -1;

      if (prev_dp.first > 0) {
        dp[idx] = prev_dp.first + k_length;
        recon[idx] = prev_dp.second;
      }
    } else {
      if (continues[idx] != -1) {
        if (dp[continues[idx]] + 1 > dp[idx]) {
          dp[idx] = dp[continues[idx]] + 1;
          recon[idx] = continues[idx];
        }
      }

      dp_col_max.update(j, std::make_pair(dp[idx], idx));

      if (dp[idx] > lcskpp_length) {
        lcskpp_length = dp[idx];
        best_idx = idx;
      }
    }
  }

  ret_lcskpp_indices->clear();
  ret_lcskpp_indices->reserve(num_matches);

  if (best_idx != -1 && recon.size() > 0) {
    ret_lcskpp_indices->push_back(matches_indices[best_idx]);

    for (int i1 = best_idx; i1 != -1; i1 = recon[i1]) {
      if (recon[i1] != -1) {
        ret_lcskpp_indices->push_back(matches_indices[recon[i1]]);
      }
    }
  }

  *ret_lcskpp_length = lcskpp_length;

  if (events)
    free(events);
  events = NULL;
  if (matches_starts)
    free(matches_starts);
  matches_starts = NULL;
  if (matches_dists_ref)
    free(matches_dists_ref);
  matches_dists_ref = NULL;
  if (matches_indices)
    free(matches_indices);
  matches_indices = NULL;
}










void WriteHits(std::string out_path, std::vector<SeedHit2> &seed_hits, int64_t hits_start, int64_t hits_end, int64_t ref_id, std::string read_header, int64_t read_length, std::string reference_header, int64_t reference_length, std::vector<int> *indices_to_output, std::vector<int32_t> *cluster_ids) {
  std::vector<SeedHit2> filtered_seed_hits;

  if (indices_to_output != NULL) {
    filtered_seed_hits.resize(indices_to_output->size());
    for (int64_t j = 0; j < indices_to_output->size(); j++) {
      filtered_seed_hits[j] = seed_hits[(*indices_to_output)[j] + hits_start];
    }
  } else {
    filtered_seed_hits.resize((hits_end - hits_start + 1));
    for (int64_t j = 0; j < (hits_end - hits_start + 1); j++) {
      filtered_seed_hits[j] = seed_hits[j + hits_start];
    }
  }

  FILE *fp = fopen(out_path.c_str(), "w");
  if (fp != NULL) {
    fprintf (fp, "%s\t0\t%ld\t%s\t0\t%ld\t0.0\n", read_header.c_str(), read_length, reference_header.c_str(), reference_length);
    for (int64_t j=0; j<filtered_seed_hits.size(); j++) {
      uint32_t seed_len = 12; // owler_data->seed_types[filtered_seed_hits[j].seed_type].length();
      int32_t cluster_id = 0;
      if (cluster_ids) {
        cluster_id = cluster_ids->at(j);
      }
      fprintf (fp, "%ld\t%ld\t%d\n", filtered_seed_hits[j].query_pos, filtered_seed_hits[j].ref_pos, cluster_id);
      fprintf (fp, "%ld\t%ld\t%d\n", filtered_seed_hits[j].query_pos + seed_len, filtered_seed_hits[j].ref_pos + seed_len, cluster_id);
    }
    fclose(fp);
  }
}

int Owler::CalcCoveredBases(std::vector<SeedHit2> &seed_hits, int64_t seed_length, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t *ret_cov_A, int64_t *ret_cov_B) {
  int64_t cov_bases = 0, cov_bases_A = 0, cov_bases_B = 0;

//  for (int64_t i = 0; i < lcskpp_indices.size(); i++) {
//    int64_t A_start = seed_hits[hits_start + lcskpp_indices[i]].query_pos;
//    int64_t B_start = seed_hits[hits_start + lcskpp_indices[i]].ref_pos;
//    printf ("[%ld] [A_start, B_start] = [%ld, %ld]\n", i, A_start, B_start);
//    fflush(stdout);
//  }
//  printf ("\n");
//  fflush(stdout);

  /// End coordinates will be inclusive, as well as start coordinates.
  int64_t A_end_prev = -1;
  int64_t B_end_prev = -1;
  for (int64_t i = 0; i < lcskpp_indices.size(); i++) {
//  for (int64_t i=(lcskpp_indices.size()-1); i>=0; i--) {
    int64_t A_start = seed_hits[hits_start + lcskpp_indices[i]].query_pos;
    int64_t A_end = A_start + seed_length - 1;
    int64_t B_start = seed_hits[hits_start + lcskpp_indices[i]].ref_pos;
    int64_t B_end = B_start + seed_length - 1;

    cov_bases_A += std::min((A_end - A_end_prev), seed_length);
    cov_bases_B += std::min((B_end - B_end_prev), seed_length);
//    if (A_end_prev == -1 && B_end_prev == -1) {
//      cov_bases_A += seed_length;
//      cov_bases_B += seed_length;
//    } else {
//      cov_bases_A += std::min((A_end - A_end_prev), seed_length);
//      cov_bases_B += std::min((B_end - B_end_prev), seed_length);
//    }

//    if (std::min((A_end - A_end_prev), seed_length) < 0) {
//      printf ("i = %ld, std::min((A_end - A_end_prev), seed_length) = %ld < 0, A_end = %ld, A_end_prev = %ld, B_end = %ld, B_end_prev = %ld, seed_length = %ld\n", i, std::min((A_end - A_end_prev), seed_length), A_end, A_end_prev, B_end, B_end_prev, seed_length);
//      fflush(stdout);
//    }

    if (cov_bases_A < 0 || cov_bases_B < 0) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("ERROR: Number of covered bases should not be less than 0! i = %ld, std::min((A_end - A_end_prev), seed_length) = %ld < 0, A_end = %ld, A_end_prev = %ld, B_end = %ld, B_end_prev = %ld, seed_length = %ld\n", i, std::min((A_end - A_end_prev), seed_length), A_end, A_end_prev, B_end, B_end_prev, seed_length), "CalcCoveredBases");
    }

    A_end_prev = A_end;
    B_end_prev = B_end;
  }

  *ret_cov_A = cov_bases_A;
  *ret_cov_B = cov_bases_B;

  return 0;
}

OverlapResult Owler::GenerateOverlapResult(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t ref_id, int64_t reference_length, bool ref_reversed, std::string ref_header,
                          int64_t read_id, int64_t read_length, bool read_reversed, std::string read_header) {
  OverlapResult ret;

//  int64_t front_id = lcskpp_indices.front();
//  int64_t back_id = lcskpp_indices.back();
  int64_t front_id = lcskpp_indices.back();
  int64_t back_id = lcskpp_indices.front();

  int64_t A_start = seed_hits[hits_start + back_id].query_pos;
  int64_t A_end = seed_hits[hits_start + front_id].query_pos + 12;
  int64_t B_start = seed_hits[hits_start + back_id].ref_pos;
  int64_t B_end = seed_hits[hits_start + front_id].ref_pos + 12;

  if (ref_reversed) {
    int64_t temp = B_start;
    B_start = reference_length - B_end - 1;
    B_end = reference_length - temp - 1;
  }

  int64_t shared_minmers = lcskpp_indices.size();

//  int64_t covered_bases = shared_minmers * (12 + 1);
  int64_t cov_bases_A = 0, cov_bases_B = 0;

//  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, "Calling function CalcCoveredBases.\n", "GenerateOverlapResult");
  CalcCoveredBases(seed_hits, 13, lcskpp_indices, hits_start, hits_end, &cov_bases_A, &cov_bases_B);

//  float jaccard_score = ((float) covered_bases) / ((float) std::min((A_end - A_start), (B_end - B_start)));
  float jaccard_score = std::max(((float) cov_bases_A) / ((float) (A_end - A_start)), ((float) cov_bases_B) / ((float) (B_end - B_start)));

  ret.read_id = read_id;
  ret.ref_id = ref_id;
  ret.jaccard_score = jaccard_score;
  ret.shared_minmers = shared_minmers;
  ret.read_is_reverse = read_reversed;
  ret.read_start = A_start;
  ret.read_end = A_end;
  ret.read_length = read_length;
  ret.ref_is_reverse = ref_reversed;
  ret.ref_start = B_start;
  ret.ref_end = B_end;
  ret.ref_length = reference_length;

  ret.read_header = read_header;
  ret.ref_header = ref_header;
  ret.cov_bases_read = cov_bases_A;
  ret.cov_bases_ref = cov_bases_B;

  ret.front_id = front_id;
  ret.back_id = back_id;

  return ret;
}

//std::string Owler::OverlapToMHAP(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t ref_id, int64_t reference_length, bool ref_reversed,
//                          int64_t read_id, int64_t read_length) {
//  std::stringstream ret;
//
////  int64_t front_id = lcskpp_indices.front();
////  int64_t back_id = lcskpp_indices.back();
//  int64_t front_id = lcskpp_indices.back();
//  int64_t back_id = lcskpp_indices.front();
//
//  int64_t A_start = seed_hits[hits_start + back_id].query_pos;
//  int64_t A_end = seed_hits[hits_start + front_id].query_pos + 12;
//  int64_t B_start = seed_hits[hits_start + back_id].ref_pos;
//  int64_t B_end = seed_hits[hits_start + front_id].ref_pos + 12;
//
//  if (ref_reversed) {
//    int64_t temp = B_start;
//    B_start = reference_length - B_end - 1;
//    B_end = reference_length - temp - 1;
//  }
//
//  int64_t shared_minmers = lcskpp_indices.size();
//
////  int64_t covered_bases = shared_minmers * (12 + 1);
//  int64_t cov_bases_A = 0, cov_bases_B = 0;
////  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, "Calling function CalcCoveredBases.\n", "OverlapToMHAP");
//  CalcCoveredBases(seed_hits, 13, lcskpp_indices, hits_start, hits_end, &cov_bases_A, &cov_bases_B);
//
////  float jaccard_score = ((float) covered_bases) / ((float) std::min((A_end - A_start), (B_end - B_start)));
//  float jaccard_score = std::max(((float) cov_bases_A) / ((float) (A_end - A_start)), ((float) cov_bases_B) / ((float) (B_end - B_start)));
//
//  ret << read_id << " ";      /// read1_id
//  ret << ref_id << " ";      /// read2_id
//  ret << jaccard_score << " ";      /// Jaccard score
//  ret << shared_minmers << " ";        /// Shared minmers
//  ret << "0 ";        /// A is reverse
//  ret << A_start << " ";
//  ret << A_end << " ";
//  ret << read_length << " ";
//  ret << (ref_reversed ? 1 : 0) << " ";
//  ret << B_start << " ";
//  ret << B_end << " ";
//  ret << reference_length;
//
//  return ret.str();
//}
//
//std::string Owler::OverlapToAFG(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t ref_id, int64_t reference_length, bool ref_reversed,
//                          int64_t read_id, int64_t read_length) {
//  std::stringstream ret;
//
//  int64_t front_id = lcskpp_indices.front();
//  int64_t back_id = lcskpp_indices.back();
//
//  int64_t A_start = seed_hits[hits_start + back_id].query_pos;
//  int64_t A_end = seed_hits[hits_start + front_id].query_pos + 12;
//  int64_t B_start = seed_hits[hits_start + back_id].ref_pos;
//  int64_t B_end = seed_hits[hits_start + front_id].ref_pos + 12;
//
////  ret << read_id << " ";      /// read1_id
////  ret << ref_id << " ";      /// read2_id
////  ret << "0.0 ";      /// Jaccard score
////  ret << "0 ";        /// Shared minmers
////  ret << "0 ";        /// A is reverse
////  ret << A_start << " ";
////  ret << A_end << " ";
////  ret << read_length << " ";
////  ret << (ref_reversed ? 1 : 0) << " ";
////  ret << B_start << " ";
////  ret << B_end << " ";
////  ret << reference_length;
//
//
//
//  std::string adj = (ref_reversed == false) ? OVERLAP_NORMAL : OVERLAP_INNIE;
//  int64_t read1_id = read_id + 1;
//  int64_t read2_id = ref_id + 1;
//  int64_t score = std::min((A_end - A_start), (B_end - B_start));
//
//  //ahg - Ahang. Length of the non-overlapping portion of the first read.
//  //bhg - Bhang. Length of the non-overlapping portion of the second read.
//  /// The position (alignment_start - clip_count_front) would be roughly where alignment of one reads starts on another.
//  /// If this value is > 0, the first read starts within read2, and thus ahang needs to be negative (hence the '-' sign).
//  int64_t ahang = 0; // - (alignment_start - clip_count_front);
//  int64_t bhang = 0; // reference_length - (alignment_end + clip_count_back);
//
//  ret << "{OVL" << "\n";
//  ret << "adj:" << adj << "\n";
//  ret << "rds:" << read1_id << "," << read2_id << "\n";
//  ret << "scr:" << score << "\n";
//  ret << "ahg:" << ahang << "\n";
//  ret << "bhg:" << bhang << "\n";
//  ret << "}";
//
//
//
//  return ret.str();
//}

bool CheckContainment(int64_t A_start, int64_t A_end, int64_t A_length, int64_t B_start, int64_t B_end, int64_t B_length) {
  /// Check if A is contained in B
  if ((B_start - A_start) >= 0 && (B_end + (A_length - A_end)) < B_length)
    return true;

  /// Check if B is contained in A
  if ((A_start - B_start) >= 0 && (A_end + (B_length - B_end)) < A_length)
    return true;

  return false;
}

bool CheckValidOverlap(int64_t A_start, int64_t A_end, int64_t A_length, int64_t B_start, int64_t B_end, int64_t B_length) {
  /// Check if A is contained in B
  if ((B_start - A_start) >= 0 && (B_end + (A_length - A_end)) < B_length) {
    // In this case, the length of the overlap should be compared to the entire length of overhangs (not-overlapping regions of A).
    int64_t non_overlapping_len = A_start + (A_length - A_end + 1);
    int64_t overlapping_len = (A_end - A_start + 1);

    if (non_overlapping_len > overlapping_len/2) {
//      printf ("Failed 1!\n");
//      fflush(stdout);
      return false;
    }
  }
  else if ((A_start - B_start) >= 0 && (A_end + (B_length - B_end)) < A_length) {  /// Check if B is contained in A
    // In this case, the length of the overlap should be compared to the entire length of overhangs (not-overlapping regions of B).
    int64_t non_overlapping_len = B_start + (B_length - B_end + 1);
    int64_t overlapping_len = (B_end - B_start + 1);

    if (non_overlapping_len > overlapping_len/2) {
//      printf ("Failed 2!\n");
//      fflush(stdout);
      return false;
    }
  }

  /// Check if suffix of A is prefix of B.
  if ((B_start - A_start) < 0 && (B_end + (A_length - A_end)) < B_length) {
    int64_t non_overlapping_len = B_start + // This is the part not covered with kmer hits on the prefix of B
                                  (A_length - A_end + 1); // This is the part of A at the end, not covered with kmer hits.
    int64_t overlapping_len = (A_end - A_start + 1);
    if (non_overlapping_len > overlapping_len/2)
//      printf ("Failed 3!\n");
//      fflush(stdout);
      return false;
  }
  else if ((A_start - B_start) < 0 && (A_end + (B_length - B_end)) < A_length) { /// Check if suffix of B is prefix of A.
    int64_t non_overlapping_len = A_start + // This is the part not covered with kmer hits on the prefix of B
                                  (B_length - B_end + 1); // This is the part of A at the end, not covered with kmer hits.
    int64_t overlapping_len = (B_end - B_start + 1);
    if (non_overlapping_len > overlapping_len/2) {
//      printf ("Failed 4!\n");
//      fflush(stdout);
      return false;
    }
  }

  return true;
}

bool CheckAPrefixB(int64_t A_start, int64_t A_end, int64_t A_length, int64_t B_start, int64_t B_end, int64_t B_length) {
  /// Check if A is prefix of B.
//  if ((B_start - A_start) < 0 && (B_end + (A_length - A_end)) < B_length)
  if ((B_start - A_start) < 0 && ((B_length - B_end) > (A_length - A_end)))
    return true;

  return false;
}

bool CheckAInB(int64_t A_start, int64_t A_end, int64_t A_length, int64_t B_start, int64_t B_end, int64_t B_length) {
  /// Check if A is contained in B.
  if (A_start < B_start && (A_length - A_end) < (B_length - B_end))
    return true;

  return false;
}

int Owler::OverlapLength(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t *A_start, int64_t *A_end, int64_t *ret_query_length, int64_t *B_start, int64_t *B_end, int64_t *ret_ref_length) {
  if (lcskpp_indices.size() < 2) {
    *ret_query_length = 0;
    *ret_ref_length = 0;
    return 1;
  }

  int64_t front_id = lcskpp_indices.back();
  int64_t back_id = lcskpp_indices.front();

  *A_start = seed_hits[hits_start + back_id].query_pos;
  *A_end = seed_hits[hits_start + front_id].query_pos + 12;
  *B_start = seed_hits[hits_start + back_id].ref_pos;
  *B_end = seed_hits[hits_start + front_id].ref_pos + 12;

  *ret_query_length = ((*A_end) - (*A_start));
  *ret_ref_length = (*(B_end) - *(B_start));

  return 0;
}

std::string Owler::OverlapToDot(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t ref_id, int64_t reference_length, bool ref_reversed,
                          int64_t read_id, int64_t read_length) {
  std::stringstream ret;

//  int64_t front_id = lcskpp_indices.front();
//  int64_t back_id = lcskpp_indices.back();
  int64_t front_id = lcskpp_indices.back();
  int64_t back_id = lcskpp_indices.front();

  int64_t A_start = seed_hits[hits_start + back_id].query_pos;
  int64_t A_end = seed_hits[hits_start + front_id].query_pos + 12;
  int64_t A_length = read_length;
  int64_t B_start = seed_hits[hits_start + back_id].ref_pos;
  int64_t B_end = seed_hits[hits_start + front_id].ref_pos + 12;
  int64_t B_length = reference_length;

//  ret << read_id << " ";      /// read1_id
//  ret << ref_id << " ";      /// read2_id
//  ret << "0.0 ";      /// Jaccard score
//  ret << "0 ";        /// Shared minmers
//  ret << "0 ";        /// A is reverse
//  ret << A_start << " ";
//  ret << A_end << " ";
//  ret << read_length << " ";
//  ret << (ref_reversed ? 1 : 0) << " ";
//  ret << B_start << " ";
//  ret << B_end << " ";
//  ret << reference_length;

//  ret << read_id << " -- " << ref_id << " [dir=" << "];";

  /// Remove the containment edges from the output, just to make the graph more clear.
  if (CheckValidOverlap(A_start, A_end, A_length, B_start, B_end, B_length)) {
//    if (CheckContainment(A_start, A_end, A_length, B_start, B_end, B_length) == false) {
      if (CheckAPrefixB(A_start, A_end, A_length, B_start, B_end, B_length)) {
  //      ret << "\t" << "Q" << read_id << "_" << A_length << " -> R" << ref_id << "_" << B_length << " [label=\"" << std::min((A_end - A_start), (B_end, B_start)) << "\"]" << ";";
  //      ret << "\t" << "Q" << read_id << "_" << A_length << " -> R" << ref_id << "_" << B_length << ";";
        if (ref_reversed == false)
          ret << "\t" << "Q" << read_id << " -> R" << ref_id << " [style=solid,arrowhead=normal,dir=normal]" << ";";
        else
          ret << "\t" << "Q" << read_id << " -> rev_R" << ref_id << " [style=solid,arrowhead=normal,dir=normal]" << ";";

        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, "\t- Q is the prefix of R\n", "[]");

      }
      else if (CheckAPrefixB(B_start, B_end, B_length, A_start, A_end, A_length)) {
  //      ret << "\t" << "R" << ref_id << "_" << B_length << " -> Q" << read_id << "_" << A_length << " [label=\"" << std::min((A_end - A_start), (B_end, B_start)) << "\"]" << ";";
  //      ret << "\t" << "R" << ref_id << "_" << B_length << " -> Q" << read_id << "_" << A_length << ";";
        if (ref_reversed == false)
          ret << "\t" << "R" << ref_id << " -> Q" << read_id << " [style=solid,arrowhead=normal,dir=normal]" << ";";
        else
          ret << "\t" << "rev_R" << ref_id << " -> Q" << read_id << " [style=solid,arrowhead=normal,dir=normal]" << ";";

        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, "\t- R is the prefix of Q\n", "[]");

      } else {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, "\t- one read is contained in the other\n", "[]");
        if (CheckAInB(A_start, A_end, A_length, B_start, B_end, B_length)) {
          if (ref_reversed == false)
            ret << "\t" << "R" << ref_id << " -> Q" << read_id << " [style=dotted,arrowhead=diamond,dir=normal]" << ";";
          else
            ret << "\t" << "rev_R" << ref_id << " -> Q" << read_id << " [style=dotted,arrowhead=diamond,dir=normal]" << ";";
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, "\t- R is contained in Q\n", "[]");

        } else if (CheckAInB(B_start, B_end, B_length, A_start, A_end, A_length)) {
          if (ref_reversed == false)
            ret << "\t" << "Q" << read_id << " -> R" << ref_id << " [style=dotted,arrowhead=diamond,dir=normal]" << ";";
          else
            ret << "\t" << "Q" << read_id << " -> rev_R" << ref_id << " [style=dotted,arrowhead=diamond,dir=normal]" << ";";
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, "\t- Q is contained in R\n", "[]");

        } else {
          ret << "\t" << "Q" << read_id << " -> R" << ref_id << " [style=dotted,dir=both,arrowhead=diamond]" << ";";
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, "\t- Q and R are completely overlapping\n", "[]");
        }
      }

//    } else {
//  //    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//        LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, "\t- one read is contained in the other\n", "[]");
//  //    }
//    }
  } else {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, "\t- overlap is not valid!\n", "[]");
  }



//  std::string adj = (ref_reversed == false) ? OVERLAP_NORMAL : OVERLAP_INNIE;
//  int64_t read1_id = read_id + 1;
//  int64_t read2_id = ref_id + 1;
//  int64_t score = std::min((A_end - A_start), (B_end - B_start));
//
//  //ahg - Ahang. Length of the non-overlapping portion of the first read.
//  //bhg - Bhang. Length of the non-overlapping portion of the second read.
//  /// The position (alignment_start - clip_count_front) would be roughly where alignment of one reads starts on another.
//  /// If this value is > 0, the first read starts within read2, and thus ahang needs to be negative (hence the '-' sign).
//  int64_t ahang = - (alignment_start - clip_count_front);
//  int64_t bhang = reference_length - (alignment_end + clip_count_back);
//
//  ret << "{OVL" << "\n";
//  ret << "adj:" << adj << "\n";
//  ret << "rds:" << read1_id << "," << read2_id << "\n";
//  ret << "scr:" << score << "\n";
//  ret << "ahg:" << ahang << "\n";
//  ret << "bhg:" << bhang << "\n";
//  ret << "}";



  return ret.str();
}



std::string Owler::OverlapMHAPVerbose(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters, int64_t ref_id, int64_t hits_start, std::vector<int> &lcskpp_indices) {
  std::stringstream ss;

  if (lcskpp_indices.size() < 2)
    return std::string("");

  int64_t front_id = lcskpp_indices.back();
  int64_t back_id = lcskpp_indices.front();

  int64_t A_start = owler_data->seed_hits2[hits_start + back_id].query_pos;
  int64_t A_end = owler_data->seed_hits2[hits_start + front_id].query_pos + 12;
  int64_t B_start = owler_data->seed_hits2[hits_start + back_id].ref_pos;
  int64_t B_end = owler_data->seed_hits2[hits_start + front_id].ref_pos + 12;

  ss << " # ";
  ss << "qry:" << read->get_header() << "\t";
  ss << "ref:" << indexes[0]->get_headers()[ref_id] << "\t";
  ss << "num_anchors:" << lcskpp_indices.size() << "\t";
  ss << "Adiff:" << (A_end - A_start) << "\t";
  ss << "Bdiff:" << (B_end - B_start) << "\n";

  return ss.str();
}

int CalcNonOverlapLength(int64_t A_start, int64_t A_end, int64_t A_len, int64_t B_start, int64_t B_end, int64_t B_len, int64_t *ret_overlap_len, int64_t *ret_nonoverlap_len_start, int64_t *ret_nonoverlap_len_end) {
  int64_t overlap_len = sqrt((float) ((A_end - A_start)*(A_end - A_start) + (B_end - B_start)*(B_end - B_start)) );
  return 0;
}

int Owler::ApplyLCS2(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  std::sort(owler_data->seed_hits2.begin(), owler_data->seed_hits2.end(), seedhits2_refid_less_than_key());
//  int64_t min_num_hits = std::min(100.0f, 0.10f * read->get_sequence_length());
//  int64_t min_num_hits = std::min(33.0f, 0.05f * read->get_sequence_length());
//  int64_t min_num_hits = std::min(50.0f, 0.10f * read->get_sequence_length());
//  int64_t min_num_hits = std::max(50.0f, 0.05f * read->get_sequence_length());
  int64_t min_num_hits = std::min(10.0f, 0.05f * read->get_sequence_length());
//  float max_overhang_percent = 0.10f; // 0.33f;
//  float min_perc_covered_bases = 0.10f;
//  float min_perc_overlap_len = 0.10f;
  float max_overhang_percent = 0.33f;
//  float max_overhang_percent = 0.10f;
  float min_perc_covered_bases = 0.10f;
  float min_perc_overlap_len = 0.0f;
  int64_t seed_length = 13;

  Index *index = indexes[0];

  int64_t read_id = read->get_sequence_id();
  int64_t read_length = read->get_sequence_length();
  int64_t num_references_fwd = index->get_num_sequences_forward();
  std::string read_header = read->get_header();

  int64_t previous_ref_id = 0, current_ref_id = 0, next_ref_id = 0;
  int64_t ref_streak_start = 0;

  int64_t num_output_overlaps = 0;

  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");
  }

  std::vector<OverlapResult> found_overlaps;

  for (int64_t i = 0; i < owler_data->seed_hits2.size(); i++) {
//    if (i > 1000)
//      break;

//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[%ld] owler_data->seed_hits2.size() = %ld\n", i, owler_data->seed_hits2.size()), "[]");

    if (i == 0)
      ref_streak_start = 0;
    current_ref_id = owler_data->seed_hits2[i].ref_id;
    int64_t ref_length = index->get_reference_lengths()[current_ref_id];

    if ((i + 1) == owler_data->seed_hits2.size() || ((i + 1) < owler_data->seed_hits2.size() && owler_data->seed_hits2[i+1].ref_id != current_ref_id)) {
      int64_t read_len = read->get_sequence_length();
      int64_t ref_id = owler_data->seed_hits2[i].ref_id;
      int64_t ref_len = index->get_reference_lengths()[ref_id];
      std::string ref_header = index->get_headers()[ref_id % num_references_fwd];

      /// Check if there are enough hits on a reference to initiate the LCSk process.
      if ((i - ref_streak_start + 1) >= min_num_hits) {
        std::stringstream debug_verbose;

        int64_t ref_streak_end = i;

        int32_t lcskpp_length;
        std::vector<int> raw_lcskpp_indices;
        CalcLCSFromLocalScoresCacheFriendly2_(owler_data, ref_streak_start, ref_streak_end, &lcskpp_length, &raw_lcskpp_indices);

        std::vector<int> lcskpp_indices;
        std::vector<int32_t> cluster_ids;
        int num_svs = FilterAnchorBreakpoints(raw_lcskpp_indices, ref_streak_start, ref_streak_end, seed_length, 0.01f*read->get_sequence_length(), 0.01f, owler_data, indexes, read, parameters, lcskpp_indices, &cluster_ids);
        num_svs = 0;

        int64_t A_start = 0, A_end = 0, query_overlap_length = 0, B_start = 0, B_end = 0, ref_overlap_length = 0;
        OverlapLength(owler_data->seed_hits2, lcskpp_indices, ref_streak_start, ref_streak_end, &A_start, &A_end, &query_overlap_length, &B_start, &B_end, &ref_overlap_length);

        int64_t cov_bases_read = 0, cov_bases_ref = 0;
//        LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, "Calling function CalcCoveredBases.\n", "ApplyLCS2");
        CalcCoveredBases(owler_data->seed_hits2, seed_length, lcskpp_indices, ref_streak_start, ref_streak_end, &cov_bases_read, &cov_bases_ref);

        float size_diff = 1.0f - ((float) std::min(query_overlap_length, ref_overlap_length)) / ((float) std::max(query_overlap_length, ref_overlap_length));
        float perc_covered_bases_read = (query_overlap_length > 0) ? (((float) cov_bases_read) / ((float) query_overlap_length)) : 0.0f;
        float perc_covered_bases_ref = (ref_overlap_length > 0) ? (((float) cov_bases_ref) / ((float) ref_overlap_length)) : 0.0f;


        bool overhang_ok = true;
        /// Test filter - some reads might span over repeat regions with great hits, but only in the middle of those reads. Limit the allowed overhangs compared to the matching overlap length.
//        if (A_start > (query_overlap_length * 0.33f) && (read_length - A_end) > (query_overlap_length * 0.33f))
//          overhang_ok = false;
//        if (B_start > (ref_overlap_length * 0.33f) && (ref_length - B_end) > (ref_overlap_length * 0.33f))
//          overhang_ok = false;

        // Not perfect, because one read can start in the middle of another but still have a large chunk at the beginning uncovered.
//        if (A_start > (query_overlap_length * max_overhang_percent) && (read_length - A_end) > (query_overlap_length * max_overhang_percent) &&
//            B_start > (ref_overlap_length * max_overhang_percent) && (ref_length - B_end) > (ref_overlap_length * max_overhang_percent)) {
//          overhang_ok = false;
//
////          printf ("A_start > (query_overlap_length * max_overhang_percent) = %d, A_start = %ld, (query_overlap_length * max_overhang_percent) = %f\n", (A_start > (query_overlap_length * max_overhang_percent)), A_start, (query_overlap_length * max_overhang_percent));
////          printf ("A_start = %ld\n", A_start);
////          printf ("(query_overlap_length * max_overhang_percent) = %f\n", (query_overlap_length * max_overhang_percent));
////          printf ("(read_length - A_end) = %ld\n", (read_length - A_end));
////          printf ("(query_overlap_length * max_overhang_percent) = %f\n", (query_overlap_length * max_overhang_percent));
////
////          printf ("B_start = %ld\n", B_start);
////          printf ("(ref_overlap_length * max_overhang_percent) = %f\n", (ref_overlap_length * max_overhang_percent));
////          printf ("(ref_length - B_end) = %ld\n", (ref_length - B_end));
////          printf ("(ref_overlap_length * max_overhang_percent) = %f\n", (ref_overlap_length * max_overhang_percent));
////
////          fflush(stdout);
//        }

//        if ((A_start > (query_overlap_length * max_overhang_percent) || (read_length - A_end) > (query_overlap_length * max_overhang_percent)) &&
//            B_start <= (ref_overlap_length * max_overhang_percent) && (ref_length - B_end) <= (ref_overlap_length * max_overhang_percent))
//          overhang_ok = false;
//        if ((B_start > (query_overlap_length * max_overhang_percent) || (read_length - B_end) > (query_overlap_length * max_overhang_percent)) &&
//            A_start <= (ref_overlap_length * max_overhang_percent) && (ref_length - B_end) <= (ref_overlap_length * max_overhang_percent))
//          overhang_ok = false;
//        if ((A_start > (query_overlap_length * max_overhang_percent) &&
//                B_start > (ref_overlap_length * max_overhang_percent) &&
//                (read_length - A_end) <= (query_overlap_length * max_overhang_percent) &&
//                (ref_length - B_end) <= (ref_overlap_length * max_overhang_percent)) ||
//            (A_start <= (query_overlap_length * max_overhang_percent) &&
//                B_start <= (ref_overlap_length * max_overhang_percent) &&
//                (read_length - A_end) > (query_overlap_length * max_overhang_percent) &&
//                (ref_length - B_end) > (ref_overlap_length * max_overhang_percent)))
//          overhang_ok = false;

//        int64_t max_query_overhang = query_overlap_length * max_overhang_percent;
//        int64_t max_ref_overhang = ref_overlap_length * max_overhang_percent;
//        if (((A_start < max_query_overhang || B_start < max_ref_overhang) &&
//              (read_length - A_end) > max_query_overhang && (ref_length - B_end) > max_ref_overhang) ||
//            (((read_length - A_end) <= max_query_overhang || (ref_length - B_end) > max_ref_overhang) &&
//              (A_start >= max_query_overhang && B_start >= max_ref_overhang)))
//          overhang_ok = false;

//        int64_t dist_start = std::min(A_start, B_start);
//        int64_t dist_end = std::min((read_length - A_end), (ref_length - B_end));
//        int64_t max_overhang = std::max((query_overlap_length * max_overhang_percent), ref_overlap_length * max_overhang_percent);
//        if (dist_start > max_overhang || dist_end > max_overhang) {
//          overhang_ok = false;
////          printf ("dist_start = %ld\n", dist_start);
////          printf ("max_overhang = %ld\n", max_overhang);
////          printf ("dist_end = %ld\n", dist_end);
////          printf ("max_overhang = %ld\n", max_overhang);
////          fflush(stdout);
//        }

        int64_t max_overhang_A = query_overlap_length * max_overhang_percent;
        int64_t max_overhang_B = ref_overlap_length * max_overhang_percent;
        if ((A_start > max_overhang_A && B_start > max_overhang_B) ||
            ((read_length - A_end) > max_overhang_A && (ref_length - B_end) > max_overhang_B)) {
          overhang_ok = false;
        }






        /// Testing filter - small overhangs can be a result of indels. Simply checking the overlap length is not enough without alignment, because we do not know how good the alignment is.
        /// For testing purposes - limit the overhang length to 2% of read length. If on both ends the overhang of a read is less than that, it will be called contained.
//        if (A_start < (query_overlap_length * 0.02f) && (read_length - A_end) < (query_overlap_length * 0.02f) && B_start < (ref_overlap_length * 0.02f) && (ref_length - B_end) < (ref_overlap_length * 0.02f))
//          overhang_ok = false;
//        if (abs(A_start - B_start) < std::min(read_length, ref_length) * 0.02f && abs((read_length - A_end) - (ref_length - B_end)) < std::min(read_length, ref_length) * 0.02f)
//          overhang_ok = false;



        if (num_svs == 0 && lcskpp_indices.size() >= 5 && query_overlap_length > min_perc_overlap_len*read_len && ref_overlap_length > min_perc_overlap_len*ref_len && size_diff < parameters->error_rate &&
            (perc_covered_bases_read > min_perc_covered_bases || perc_covered_bases_ref > min_perc_covered_bases) &&
            overhang_ok == true) { // min_num_hits) {


          if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
            int64_t first_lcskpp_id = ref_streak_start + lcskpp_indices.back();
            int64_t last_lcskpp_id = ref_streak_start + lcskpp_indices.front();
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[%ld] /good/ start -> %s, end -> %s, ref_id_fwd = %ld, read_len = %ld, ref_len = %ld\n", (num_output_overlaps + 1), owler_data->seed_hits2[ref_streak_start].VerboseToString().c_str(), owler_data->seed_hits2[i].VerboseToString().c_str(), (owler_data->seed_hits2[i].ref_id % index->get_num_sequences_forward()), read_len, ref_len), "[]");
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\t- start -> %s, (end) %s\n", owler_data->seed_hits2[first_lcskpp_id].VerboseToString().c_str(), owler_data->seed_hits2[last_lcskpp_id].VerboseToString().c_str()), "[]");
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\t- found different reference, [%ld, %ld]\n", ref_streak_start, i), "[]");
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\t- passed filter\n", "[]");
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\t- LCSk performed, lcskpp_indices.size() = %ld\n", lcskpp_indices.size()), "[]");
          }

          if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\t- number of SV detected: %d\n", num_svs), "[]");
          }

          found_overlaps.push_back(GenerateOverlapResult(owler_data->seed_hits2, lcskpp_indices, ref_streak_start, i, (current_ref_id % num_references_fwd), ref_length, (current_ref_id >= num_references_fwd), ref_header,
                                                         read_id, read_length, false, read_header));

//          if (parameters->outfmt == "afg") {
//            if (owler_data->overlap_lines.size() > 0)
//              owler_data->overlap_lines += "\n";
//            owler_data->overlap_lines += OverlapToAFG(owler_data->seed_hits2, lcskpp_indices, ref_streak_start, i, (current_ref_id % num_references_fwd), ref_length, (current_ref_id >= num_references_fwd),
//                                                       read_id + index->get_num_sequences(), read_length);
//
//          } else if (parameters->outfmt == "dot") {
//            if (owler_data->overlap_lines.size() > 0)
//              owler_data->overlap_lines += "\n";
//            owler_data->overlap_lines += OverlapToDot(owler_data->seed_hits2, lcskpp_indices, ref_streak_start, i, (current_ref_id % num_references_fwd), ref_length, (current_ref_id >= num_references_fwd),
//                                                       read_id, read_length);
//
//          } else {
//            if (owler_data->overlap_lines.size() > 0)
//              owler_data->overlap_lines += "\n";
//            owler_data->overlap_lines += OverlapToMHAP(owler_data->seed_hits2, lcskpp_indices, ref_streak_start, i, (current_ref_id % num_references_fwd), ref_length, (current_ref_id >= num_references_fwd),
//                                                       read_id, read_length);
//
////            if (parameters->verbose_sam_output > 0) {
////              owler_data->overlap_lines += OverlapMHAPVerbose(owler_data, indexes, read, parameters, owler_data->seed_hits2[i].ref_id, ref_streak_start, lcskpp_indices);
////            }
//          }

          num_output_overlaps += 1;
        } else {
          if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read,
                                                FormatString("[%ld] /bad/ (num_output_overlaps) current_ref_id=%ld, ref_id_fwd=%ld, overhang_ok = %d, num_svs = %d, lcskpp_indices.size() = %ld, qlen = %ld, rlen = %ld\n\t° query_overlap_length = %ld, min_perc_overlap_len*read_len = %ld\n\t° ref_overlap_length = %ld, min_perc_overlap_len*ref_len = %ld\n\t° size_diff = %f, perc_covered_bases_read = %f, perc_covered_bases_ref = %f, min_perc_covered_bases = %f\n\t° cov_bases_read = %ld, cov_bases_ref = %ld\n",
                                                             (num_output_overlaps + 1), current_ref_id, (current_ref_id % num_references_fwd), overhang_ok, num_svs, lcskpp_indices.size(), read_length, ref_len, query_overlap_length, (int64_t) min_perc_overlap_len*read_len, ref_overlap_length, (int64_t) min_perc_overlap_len*ref_len, size_diff, perc_covered_bases_read, perc_covered_bases_ref, min_perc_covered_bases, cov_bases_read, cov_bases_ref), "[]");
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read,
                                                FormatString("\t° A_start = %ld, A_end = %ld\n", A_start, A_end),"[]");
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read,
                                                FormatString("\t° B_start = %ld, B_end = %ld\n", B_start, B_end),"[]");
//            LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read,
//                                                FormatString("\t° dist_start = %ld, dist_end = %ld, max_overhang = %ld\n", dist_start, dist_end, max_overhang),"[]");

//                                                             , lcskpp_indices.size() >= 5 && query_overlap_length > min_perc_overlap_len*read_len && ref_overlap_length > min_perc_overlap_len*ref_len && size_diff < parameters->error_rate, (perc_covered_bases_read > min_perc_covered_bases || perc_covered_bases_ref > min_perc_covered_bases)), "[]");
          }
        }



      #ifndef RELEASE_VERSION
          if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
            std::string reference_header = index->get_headers()[(current_ref_id % num_references_fwd)];

            WriteHits(FormatString("temp/overlaps/scores-%ld.csv", current_ref_id), owler_data->seed_hits2, ref_streak_start, i, current_ref_id, std::string(read->get_header()), read->get_sequence_length(), reference_header, indexes[0]->get_reference_lengths()[current_ref_id], NULL, NULL);
            WriteHits(FormatString("temp/overlaps/LCS-%ld.csv", current_ref_id), owler_data->seed_hits2, ref_streak_start, i, current_ref_id, std::string(read->get_header()), read->get_sequence_length(), reference_header, indexes[0]->get_reference_lengths()[current_ref_id], &raw_lcskpp_indices, NULL);
            WriteHits(FormatString("temp/overlaps/LCSL1-%ld.csv", current_ref_id), owler_data->seed_hits2, ref_streak_start, i, current_ref_id, std::string(read->get_header()), read->get_sequence_length(), reference_header, indexes[0]->get_reference_lengths()[current_ref_id], &lcskpp_indices, &cluster_ids);
            WriteHits(FormatString("temp/overlaps/double_LCS-%ld.csv", current_ref_id), owler_data->seed_hits2, ref_streak_start, i, current_ref_id, std::string(read->get_header()), read->get_sequence_length(), reference_header, indexes[0]->get_reference_lengths()[current_ref_id], &lcskpp_indices, &cluster_ids);
          }
      #endif

      cluster_ids.clear();
      lcskpp_indices.clear();
      }

      ref_streak_start = (i + 1);
    }
  }

  std::sort(found_overlaps.begin(), found_overlaps.end(), overlapresult_sort_key());

//  std::vector<OverlapResult> final_overlaps
  for (int64_t i=0; i<found_overlaps.size(); i++) {
    if (i > 0 && found_overlaps[i-1].read_id == found_overlaps[i].read_id && found_overlaps[i-1].ref_id == found_overlaps[i].ref_id) {
      continue;
    }
    if (owler_data->overlap_lines.size() > 0)
      owler_data->overlap_lines += "\n";

    /// Generate output lines.
    if (parameters->outfmt == "mhap") {
      owler_data->overlap_lines += found_overlaps[i].GenerateMHAPLine();
    } else if (parameters->outfmt == "paf") {
      owler_data->overlap_lines += found_overlaps[i].GeneratePAFLine();
    } else if (parameters->outfmt == "afg") {
      owler_data->overlap_lines += found_overlaps[i].GenerateAFGLine();
    } else {
      owler_data->overlap_lines += found_overlaps[i].GenerateMHAPLine();
    }
  }

  return 0;
}

// Assumes vertices are sorted.
// If use_l1_filtering is true, then all vertices/anchors that have coordinates further than allowed_dist from the L1 line are filtered out.
// Otherwise, all vertices will be used.
// The L1 line is specified with k = 1 and l parameters (y = k*x + l).
void Owler::CalcLCSFromLocalScoresCacheFriendly2_(OwlerData* owler_data, int64_t ref_hits_start, int64_t ref_hits_end, int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices) {
  uint32_t num_vertices = ref_hits_end - ref_hits_start + 1;

  if (num_vertices <= 0)
    return;

  unsigned __int128 *events = NULL;
  uint64_t *matches_starts = NULL;
  uint64_t *matches_dists_ref = NULL;
  uint64_t *matches_indices = NULL;

  SeedHit2 *seed_hits = &(owler_data->seed_hits2[ref_hits_start]);

  uint32_t n = 0;
  int64_t num_matches = 0;

  int64_t num_events = 0;
  int64_t lcskpp_length = 0;

  uint32_t min_ref = seed_hits[0].ref_pos, min_query = seed_hits[0].query_pos;

  for (uint32_t i=0; i<num_vertices; i++) {
    min_ref = std::min(min_ref, seed_hits[i].ref_pos);
    min_query = std::min(min_query, seed_hits[i].query_pos);
  }

  events = (unsigned __int128 *) malloc(sizeof(unsigned __int128) * num_vertices * 2);
  matches_starts = (uint64_t *) malloc(sizeof(uint64_t) * num_vertices);
  matches_dists_ref = (uint64_t *) malloc(sizeof(uint64_t) * num_vertices);
  matches_indices = (uint64_t *) malloc(sizeof(uint64_t) * num_vertices);

  for (uint32_t i=0; i<num_vertices; i++) {
    uint32_t seed_len = 14; // owler_data->seed_types[owler_data->seed_hits2[i].seed_type].length();
    uint32_t ref_start = seed_hits[i].ref_pos - min_ref;
    uint32_t ref_end = ref_start + seed_len;
    uint32_t query_start = seed_hits[i].query_pos - min_query;
    uint32_t query_end = query_start + seed_len;

    unsigned __int128 event1 = (((unsigned __int128) ref_start) << (8 * 8)) | (((unsigned __int128) query_start) << (4 * 8)) | (((unsigned __int128) (i + num_vertices)));
    events[num_events] = event1;
    num_events += 1;

    unsigned __int128 event2 = (((unsigned __int128) ref_end) << (8 * 8)) | (((unsigned __int128) query_end) << (4 * 8)) | ((((unsigned __int128) i)));
    events[num_events] = event2;
    num_events += 1;

    matches_starts[num_matches] = (uint64_t) (event1 >> (4 * 8));
    matches_dists_ref[num_matches] = seed_len;
    matches_indices[num_matches] = i;

    num_matches += 1;

    n = std::max(n, ref_end);
    n = std::max(n, query_end);
  }

  std::sort(events, (events + num_events));

  // Indexed by column, first:dp value, second:index in matches.
  FenwickMax<std::pair<int, int> > dp_col_max(n);
  std::vector<int> dp(num_matches);
  std::vector<int> recon(num_matches);
  std::vector<int> continues(num_matches, -1);

//    matches bi trebao biti sortiran na pocetku!

  for (int64_t curr = 0; curr < num_matches; curr++) {
    uint64_t G = 0;
    uint64_t G1 = (((matches_starts[curr] >> (4 * 8)) & (0x00000000FFFFFFFF)) - 1);
    uint64_t G2 = (((matches_starts[curr]) & (0x00000000FFFFFFFF)) - 1);
    G = (G1 << (4 * 8)) | G2;

    auto prev = std::lower_bound(matches_starts, (matches_starts + num_matches), G);

    if (prev != (matches_starts + num_matches) && *prev == G) {
      continues[curr] = prev - matches_starts;
    }
  }

  int best_idx = 0;
  lcskpp_length = 0;

  for (int64_t current_event=0; current_event<num_events; current_event++) {
    int64_t raw_idx = (int64_t) ((events[current_event] & (0x00000000FFFFFFFF)));

    int64_t idx = (raw_idx >= ((int64_t) num_matches)) ? (raw_idx - ((int64_t) num_matches)) : (raw_idx);
    bool is_beginning = (raw_idx >= num_matches);
    uint64_t i = (uint64_t) ((events[current_event] >> (8 * 8)) & (0x00000000FFFFFFFF));
    uint64_t j = (uint64_t) ((events[current_event] >> (4 * 8)) & (0x00000000FFFFFFFF));
    int primary_diagonal = n - 1 + i - j;

    if (is_beginning) { // begin
      std::pair<int, int> prev_dp = dp_col_max.get(j);
      uint64_t k_length = matches_dists_ref[idx];
      dp[idx] = k_length;      // k
      recon[idx] = -1;

      if (prev_dp.first > 0) {
        dp[idx] = prev_dp.first + k_length;
        recon[idx] = prev_dp.second;
      }
    } else {
      if (continues[idx] != -1) {
        if (dp[continues[idx]] + 1 > dp[idx]) {
          dp[idx] = dp[continues[idx]] + 1;
          recon[idx] = continues[idx];
        }
      }

      dp_col_max.update(j, std::make_pair(dp[idx], idx));

      if (dp[idx] > lcskpp_length) {
        lcskpp_length = dp[idx];
        best_idx = idx;
      }
    }
  }

  ret_lcskpp_indices->clear();
  ret_lcskpp_indices->reserve(num_matches);

  if (best_idx != -1 && recon.size() > 0) {
    ret_lcskpp_indices->push_back(matches_indices[best_idx]);

    for (int i1 = best_idx; i1 != -1; i1 = recon[i1]) {
      if (recon[i1] != -1) {
        ret_lcskpp_indices->push_back(matches_indices[recon[i1]]);
      }
    }
  }

  *ret_lcskpp_length = lcskpp_length;

  if (events)
    free(events);
  events = NULL;
  if (matches_starts)
    free(matches_starts);
  matches_starts = NULL;
  if (matches_dists_ref)
    free(matches_dists_ref);
  matches_dists_ref = NULL;
  if (matches_indices)
    free(matches_indices);
  matches_indices = NULL;
}





struct ClusterAndIndices {
  Range query;
  Range ref;
  int32_t num_anchors = 0;
  int32_t coverage = 0;
  std::vector<int> lcskpp_indices;
};

bool Owler::CheckDistanceTooBig(OwlerData* owler_data, int64_t index_last, int64_t index_current, float error_rate) { // , const ProgramParameters* parameters) {
  int64_t seed_length = 12;
  int64_t distance_query = (owler_data->seed_hits2[index_current].query_pos + seed_length) - owler_data->seed_hits2[index_last].query_pos;
  int64_t distance_ref = (owler_data->seed_hits2[index_current].ref_pos + seed_length) - owler_data->seed_hits2[index_last].ref_pos;
  float max_length = ((float) std::max(distance_query, distance_ref));
  float min_length = ((float) std::min(distance_query, distance_ref));
  if ((min_length == 0 && max_length != 0) || (min_length > 0 && (max_length / min_length - 1.0f) > error_rate)) {
    return true;
  }

  return false;
}

bool Owler::CheckDistanceStep(OwlerData* owler_data, int64_t index_first, int64_t index_last, int64_t index_current, float max_diff) { // , const ProgramParameters* parameters) {
  int64_t seed_length = 12;
  int64_t distance_query = (owler_data->seed_hits2[index_current].query_pos + seed_length) - owler_data->seed_hits2[index_first].query_pos;
  int64_t distance_query_before = (owler_data->seed_hits2[index_last].query_pos + seed_length) - owler_data->seed_hits2[index_first].query_pos;

  int64_t distance_ref = (owler_data->seed_hits2[index_current].ref_pos + seed_length) - owler_data->seed_hits2[index_first].ref_pos;
  int64_t distance_ref_before = (owler_data->seed_hits2[index_last].ref_pos + seed_length) - owler_data->seed_hits2[index_first].ref_pos;

  float diff_query = ((float) distance_query) / ((float) distance_query_before);
  float diff_ref = ((float) distance_ref) / ((float) distance_ref_before);

  if (diff_query > max_diff || diff_ref > max_diff)
    return true;

  return false;
}

int Owler::FilterAnchorBreakpoints(const std::vector<int> &lcskpp_indices, int64_t ref_hits_start, int64_t ref_hits_end, int64_t seed_length, int64_t min_cluster_length, float min_cluster_coverage, OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters, std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {
//  int64_t min_cluster_length = 0;
//  int64_t min_covered_bases = std::max(30.0f, read->get_sequence_length() * 0.02f);
//  printf ("lcskpp_indices.size() = %ld\n", lcskpp_indices.size());
//  fflush(stdout);

  std::vector<ClusterAndIndices *> clusters;
  ClusterAndIndices *new_cluster = NULL;
  int64_t last_nonskipped_i = lcskpp_indices.size() + 1;
  for (int64_t i=(lcskpp_indices.size() - 1); i >= 0; i--) {
    /// Skip anchors which might be too erroneous.
    int64_t current_lcskp_index = lcskpp_indices.at(i) + ref_hits_start;
    if (CheckDistanceTooBig(owler_data, current_lcskp_index, current_lcskp_index, parameters->error_rate / 2.0f) == true) {
      continue;
    }

    if (last_nonskipped_i > lcskpp_indices.size()) {

    } else {
      /// This is going to work, because last_nonskipped_i will be set the second iteration of the loop. The value of i starts counting from int64_t i=(lcskpp_indices.size() - 1).
      int64_t previous_lcskp_index = lcskpp_indices.at(last_nonskipped_i) + ref_hits_start;

      bool wrong_to_previous1 = CheckDistanceTooBig(owler_data, previous_lcskp_index, current_lcskp_index, parameters->error_rate / 2.0f);
      bool wrong_to_previous2 = (new_cluster->lcskpp_indices.size() < 2) ? false :
                                (CheckDistanceTooBig(owler_data, new_cluster->lcskpp_indices[new_cluster->lcskpp_indices.size()-2] + ref_hits_start, current_lcskp_index, parameters->error_rate / 2.0f));
//      bool wrong_by_distance = CheckDistanceStep(owler_data, new_cluster->lcskpp_indices.front() + ref_hits_start, previous_lcskp_index, current_lcskp_index, 1.5f);
      bool wrong_by_distance = false;

      if ((wrong_to_previous1 == true && wrong_to_previous2 == true) || (new_cluster->lcskpp_indices.size() > 1 && wrong_by_distance == true)) {
//      if ((wrong_to_previous1 == true && wrong_to_previous2 == true) || (wrong_by_distance == true)) {
        /// In this case, the new point is a general outlier to the previous LCSk, because it doesn't fit nesither to the previous point, nor to the point before that.
        if (new_cluster != NULL) {
          int64_t cov_bases_read = 0, cov_bases_ref = 0;
//          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, "Calling function CalcCoveredBases (1).\n", "FilterAnchorBreakpoints");
          CalcCoveredBases(owler_data->seed_hits2, seed_length, new_cluster->lcskpp_indices, ref_hits_start, ref_hits_end, &cov_bases_read, &cov_bases_ref);
          new_cluster->coverage = std::max(cov_bases_read, cov_bases_ref);

          int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * min_cluster_coverage;

          if ((new_cluster->query.end - new_cluster->query.start + 1) >= min_cluster_length && new_cluster->coverage >= min_covered_bases) {
            clusters.push_back(new_cluster);
          } else {
//            printf ("Tu sam 1! (new_cluster->query.end - new_cluster->query.start + 1) = %ld, min_cluster_length = %ld, new_cluster->coverage = %ld, min_covered_bases = %ld, (new_cluster->query.end - new_cluster->query.start + 1) = %ld\n", (new_cluster->query.end - new_cluster->query.start + 1), min_cluster_length, new_cluster->coverage, min_covered_bases, (new_cluster->query.end - new_cluster->query.start + 1));
//            printf ("wrong_to_previous1 = %d, wrong_to_previous2 = %d, wrong_by_distance = %d\n", wrong_to_previous1, wrong_to_previous2, wrong_by_distance);
//            fflush(stdout);
            delete new_cluster;
          }
          new_cluster = NULL;
        }
      } else if (wrong_to_previous1 == true && wrong_to_previous2 == false) {
        /// In this case, the previous point was an outlier, because the new point fits better to the one before the previous one. Overwrite the previous entry in new_cluster.

        new_cluster->query.end = owler_data->seed_hits2[current_lcskp_index].query_pos + 12 - 1;
        new_cluster->ref.end = owler_data->seed_hits2[current_lcskp_index].ref_pos + 12 - 1;

        /// This should not change, as we remove 12 bases and add 12 bases.
//          new_cluster->coverage -= local_score->get_registry_entries().covered_bases_queries[previous_lcskp_index];
//          new_cluster->coverage += local_score->get_registry_entries().covered_bases_queries[current_lcskp_index];
        new_cluster->lcskpp_indices[new_cluster->lcskpp_indices.size()-1] = current_lcskp_index - ref_hits_start;

        if (new_cluster->lcskpp_indices.size() == 1) {
          new_cluster->query.start = owler_data->seed_hits2[current_lcskp_index].query_pos;
          new_cluster->ref.start = owler_data->seed_hits2[current_lcskp_index].ref_pos;
        }
        last_nonskipped_i = i;

        continue;
      }
    }

    if (new_cluster == NULL) {
      new_cluster = new ClusterAndIndices;
      new_cluster->query.start = owler_data->seed_hits2[current_lcskp_index].query_pos;
      new_cluster->ref.start = owler_data->seed_hits2[current_lcskp_index].ref_pos;
    }

    new_cluster->query.end = owler_data->seed_hits2[current_lcskp_index].query_pos + 12 - 1;
    new_cluster->ref.end = owler_data->seed_hits2[current_lcskp_index].ref_pos + 12 - 1;
    new_cluster->num_anchors += 1;
//    new_cluster->coverage += 12;
    new_cluster->lcskpp_indices.push_back(current_lcskp_index - ref_hits_start);

    last_nonskipped_i = i;
  }
  if (new_cluster != NULL) {
    int64_t cov_bases_read = 0, cov_bases_ref = 0;
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, "Calling function CalcCoveredBases (2).\n", "FilterAnchorBreakpoints");
    CalcCoveredBases(owler_data->seed_hits2, seed_length, new_cluster->lcskpp_indices, ref_hits_start, ref_hits_end, &cov_bases_read, &cov_bases_ref);
    new_cluster->coverage = std::max(cov_bases_read, cov_bases_ref);

    int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * min_cluster_coverage;

    if ((new_cluster->query.end - new_cluster->query.start + 1) >= min_cluster_length && new_cluster->coverage >= min_covered_bases) {
      clusters.push_back(new_cluster);
    } else {
      delete new_cluster;
    }
    new_cluster = NULL;
  }

  if (ret_cluster_ids) {
    ret_cluster_ids->clear();
  }

//  printf ("clusters.size() = %ld\n", clusters.size());
//  fflush(stdout);

  int ret_val = 0;
  /// Check if the leftover clusters are linear and that only outlier anchors are filtered. This is important for overlapping.
  for (int64_t i=1; i<clusters.size(); i++) {
    int64_t current_lcskp_index = clusters[i]->lcskpp_indices.front() + ref_hits_start;
    int64_t previous_lcskp_index = clusters[i-1]->lcskpp_indices.back() + ref_hits_start;

    bool wrong_to_previous1 = CheckDistanceTooBig(owler_data, previous_lcskp_index, current_lcskp_index, parameters->error_rate);
    if (wrong_to_previous1 == true) {
      ret_val += 1;
    }
  }

  ret_filtered_lcskpp_indices.clear();
//  std::vector<int> cluster_indices;

  for (int64_t i=0; i<clusters.size(); i++) {
//    int64_t cluster_length = clusters[i]->query.end - clusters[i]->query.start + 1;
//    if (cluster_length >= min_cluster_length && clusters[i]->coverage >= min_covered_bases) {
    ret_filtered_lcskpp_indices.insert(ret_filtered_lcskpp_indices.end(), clusters[i]->lcskpp_indices.begin(), clusters[i]->lcskpp_indices.end());

    /// Create indices for debugging purposes (so we can differentiate clusters).
    if (ret_cluster_ids) {
      std::vector<int32_t> cluster_indices(clusters[i]->lcskpp_indices.size(), i);
      ret_cluster_ids->insert(ret_cluster_ids->end(), cluster_indices.begin(), cluster_indices.end());
    }
//    }

    if (clusters[i])
      delete clusters[i];
  }

//  int num_clusters = clusters.size();

  clusters.clear();

//  return num_clusters;
  return ret_val;
}

int Owler::FilterAnchorBreakpointsExperimental(const std::vector<int> &lcskpp_indices, int64_t ref_hits_start, int64_t ref_hits_end, int64_t seed_length, int64_t min_cluster_length, float min_cluster_coverage, OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters, std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {
  std::vector<ClusterAndIndices *> clusters;
  ClusterAndIndices *new_cluster = NULL;
  int64_t last_nonskipped_i = lcskpp_indices.size() + 1;
  for (int64_t i=(lcskpp_indices.size() - 1); i >= 0; i--) {
    /// Skip anchors which might be too erroneous.
    int64_t current_lcskp_index = lcskpp_indices.at(i) + ref_hits_start;
    if (CheckDistanceTooBig(owler_data, current_lcskp_index, current_lcskp_index, parameters->error_rate / 2.0f) == true)
      continue;

    if (last_nonskipped_i > lcskpp_indices.size()) {

    } else {
      /// This is going to work, because last_nonskipped_i will be set the second iteration of the loop. The value of i starts counting from int64_t i=(lcskpp_indices.size() - 1).
      int64_t previous_lcskp_index = lcskpp_indices.at(last_nonskipped_i) + ref_hits_start;

      bool wrong_to_previous1 = CheckDistanceTooBig(owler_data, previous_lcskp_index, current_lcskp_index, parameters->error_rate / 2.0f);
      bool wrong_to_previous2 = (new_cluster->lcskpp_indices.size() < 2) ? false :
                                (CheckDistanceTooBig(owler_data, new_cluster->lcskpp_indices[new_cluster->lcskpp_indices.size()-2] + ref_hits_start, current_lcskp_index, parameters->error_rate / 2.0f));
//      bool wrong_by_distance = CheckDistanceStep(owler_data, new_cluster->lcskpp_indices.front(), previous_lcskp_index, current_lcskp_index, 1.5f);
      bool wrong_by_distance = false;

      if ((wrong_to_previous1 == true && wrong_to_previous2 == true) || wrong_by_distance == true) {
        /// In this case, the new point is a general outlier to the previous LCSk, because it doesn't fit nesither to the previous point, nor to the point before that.
        if (new_cluster != NULL) {
          int64_t cov_bases_read = 0, cov_bases_ref = 0;
//          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, "Calling function CalcCoveredBases (1).\n", "FilterAnchorBreakpointsExperimental");
          CalcCoveredBases(owler_data->seed_hits2, seed_length, new_cluster->lcskpp_indices, ref_hits_start, ref_hits_end, &cov_bases_read, &cov_bases_ref);
          new_cluster->coverage = std::max(cov_bases_read, cov_bases_ref);

          int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * min_cluster_coverage;

          if ((new_cluster->query.end - new_cluster->query.start + 1) >= min_cluster_length && new_cluster->coverage >= min_covered_bases) {
            clusters.push_back(new_cluster);
          } else {
            delete new_cluster;
          }
          new_cluster = NULL;
        }
      } else if (wrong_to_previous1 == true && wrong_to_previous2 == false) {
        /// In this case, the previous point was an outlier, because the new point fits better to the one before the previous one. Overwrite the previous entry in new_cluster.
        new_cluster->query.end = owler_data->seed_hits2[current_lcskp_index].query_pos + 12 - 1;
        new_cluster->ref.end = owler_data->seed_hits2[current_lcskp_index].ref_pos + 12 - 1;

        /// This should not change, as we remove 12 bases and add 12 bases.
        new_cluster->lcskpp_indices[new_cluster->lcskpp_indices.size()-1] = current_lcskp_index - ref_hits_start;

        if (new_cluster->lcskpp_indices.size() == 1) {
          new_cluster->query.start = owler_data->seed_hits2[current_lcskp_index].query_pos;
          new_cluster->ref.start = owler_data->seed_hits2[current_lcskp_index].ref_pos;
        }
        last_nonskipped_i = i;

        continue;
      }
    }

    if (new_cluster == NULL) {
      new_cluster = new ClusterAndIndices;
      new_cluster->query.start = owler_data->seed_hits2[current_lcskp_index].query_pos;
      new_cluster->ref.start = owler_data->seed_hits2[current_lcskp_index].ref_pos;
    }

    new_cluster->query.end = owler_data->seed_hits2[current_lcskp_index].query_pos + 12 - 1;
    new_cluster->ref.end = owler_data->seed_hits2[current_lcskp_index].ref_pos + 12 - 1;
    new_cluster->num_anchors += 1;
    new_cluster->lcskpp_indices.push_back(current_lcskp_index - ref_hits_start);

    last_nonskipped_i = i;
  }
  if (new_cluster != NULL) {
    int64_t cov_bases_read = 0, cov_bases_ref = 0;
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, "Calling function CalcCoveredBases (2).\n", "FilterAnchorBreakpointsExperimental");
    CalcCoveredBases(owler_data->seed_hits2, seed_length, new_cluster->lcskpp_indices, ref_hits_start, ref_hits_end, &cov_bases_read, &cov_bases_ref);
    new_cluster->coverage = std::max(cov_bases_read, cov_bases_ref);

    int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * min_cluster_coverage;

    if ((new_cluster->query.end - new_cluster->query.start + 1) >= min_cluster_length && new_cluster->coverage >= min_covered_bases) {
      clusters.push_back(new_cluster);
    } else {
      delete new_cluster;
    }
    new_cluster = NULL;
  }

  if (ret_cluster_ids) {
    ret_cluster_ids->clear();
  }






  //////////////////////
  /// Filter shady anchors by thresholding the coverage by some statistics on the clusters.
  /// Concretely, for each cluster the percentage of covered bases is calculated. Using the median and the standard deviation
  /// we set the threshold for the minimum number of covered bases in a cluster.
  /// If below this threshold, the cluster is dismissed.
  std::vector<float> cluster_coverages;
  for (int64_t i=0; i<clusters.size(); i++) {
    int64_t cluster_length = clusters[i]->query.end - clusters[i]->query.start + 1;
    int64_t covered_bases = clusters[i]->coverage;
    cluster_coverages.push_back(((float) covered_bases) / ((float) cluster_length));
  }
  std::sort(cluster_coverages.begin(), cluster_coverages.end());
  float median = 0.0f;
  int64_t num_clusters = cluster_coverages.size();
  if (num_clusters > 0) {
    median = ((num_clusters % 2) == 0) ? ((cluster_coverages[num_clusters/2] + cluster_coverages[num_clusters/2+1]) / 2.0f) :
                                         cluster_coverages[(int64_t) floor(((float) num_clusters) / 2.0f)];
  }
  float mean = 0.0f;
  for (int64_t i=0; i<cluster_coverages.size(); i++) {
    mean += cluster_coverages[i];
  }
  if (cluster_coverages.size() > 0)
    mean /= cluster_coverages.size();
  float std = 0.0f, std_med = 0.0f;
  for (int64_t i=0; i<cluster_coverages.size(); i++) {
    std += (cluster_coverages[i] - mean) * (cluster_coverages[i] - mean);
    std_med += (cluster_coverages[i] - median) * (cluster_coverages[i] - median);
  }
  if (cluster_coverages.size() > 1) {
    std /= (cluster_coverages.size() - 1);
    std_med /= (cluster_coverages.size() - 1);
  }
  std = sqrt(std);
  std_med = sqrt(std_med);
  float min_relative_cluster_coverage = median - std_med - 0.001;  /// The 0.001 is subtracted to avoid the numerical error.

  std::vector<ClusterAndIndices *> filtered_clusters;
  for (int64_t i=0; i<clusters.size(); i++) {
      int64_t cluster_length = clusters[i]->query.end - clusters[i]->query.start + 1;
      int64_t covered_bases = clusters[i]->coverage;
      float cluster_coverage = ((float) covered_bases) / ((float) cluster_length);

      if (cluster_coverage >= min_relative_cluster_coverage) {
        filtered_clusters.push_back(clusters[i]);
      } else {
        delete clusters[i];
      }
  }



  int ret_val = 0;
  /// Check if the leftover clusters are linear and that only outlier anchors are filtered. This is important for overlapping.
  for (int64_t i=1; i<filtered_clusters.size(); i++) {
    int64_t current_lcskp_index = filtered_clusters[i]->lcskpp_indices.front() + ref_hits_start;
    int64_t previous_lcskp_index = filtered_clusters[i-1]->lcskpp_indices.back() + ref_hits_start;
    bool wrong_to_previous1 = CheckDistanceTooBig(owler_data, previous_lcskp_index, current_lcskp_index, parameters->error_rate);
    if (wrong_to_previous1 == true) {
      ret_val += 1;
    }
  }

  ret_filtered_lcskpp_indices.clear();
  for (int64_t i=0; i<filtered_clusters.size(); i++) {
    ret_filtered_lcskpp_indices.insert(ret_filtered_lcskpp_indices.end(), filtered_clusters[i]->lcskpp_indices.begin(), filtered_clusters[i]->lcskpp_indices.end());
    /// Create indices for debugging purposes (so we can differentiate clusters).
    if (ret_cluster_ids) {
      std::vector<int32_t> cluster_indices(filtered_clusters[i]->lcskpp_indices.size(), i);
      ret_cluster_ids->insert(ret_cluster_ids->end(), cluster_indices.begin(), cluster_indices.end());
    }
    if (filtered_clusters[i])
      delete filtered_clusters[i];
  }


  clusters.clear();
  filtered_clusters.clear();

  return ret_val;
}
