/*
 * core_graphmap.cc
 *
 *  Created on: Mar 19, 2015
 *      Author: isovic
 */

#include "graphmap/graphmap.h"



int GraphMap::GraphMap_(ScoreRegistry* local_score, Index *index_read, MappingData* mapping_data, const Index* index, const Index* index_secondary, const SingleSequence* read, const ProgramParameters* parameters) {
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("Entered function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024)), "GraphMap_");

  uint64_t readlength = read->get_sequence_length();

  #ifndef RELEASE_VERSION
    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n"), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Region:"), "Inspecting region.is_split == false");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("ref_id = %ld, region_id = %ld, region_votes = %ld\n", local_score->get_region().reference_id, local_score->get_region().region_index, local_score->get_region().region_votes), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("rname: %s\n", index->get_headers()[local_score->get_region().reference_id % index->get_num_sequences_forward()].c_str()), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("read_length = %ld\n", read->get_sequence_length()), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("start = %ld, end = %ld\n", local_score->get_region().start, local_score->get_region().end), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("is_split = %s\n", ((local_score->get_region().is_split == true) ? ("true") : ("false"))), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("split_start = %ld, split_end = %ld\n", local_score->get_region().split_start, local_score->get_region().split_end), "[]");
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n\n"), "[]");
    }
  #endif

  const int8_t *data_ptr = NULL;
  int8_t *data_copy = NULL;
  int64_t data_start = 0;
  int64_t data_end = 0;

  // Set the data pointers and start/end locations of the regions.
  // This part takes care of the split regions.
  if (local_score->get_region().is_split == false) {
    local_score->Reserve((local_score->get_region().end - local_score->get_region().start) * 2);
    data_ptr = index_->get_data();
    data_start = local_score->get_region().start;
    data_end = local_score->get_region().end - parameters->k_graph + 1;

  } else {
    int64_t region_length_joined = 0, start_offset = 0, position_of_ref_end = 0;
    ConcatenateSplitRegion(index_, (Region &) local_score->get_region(), &data_copy, &region_length_joined, &start_offset, &position_of_ref_end);
    local_score->Reserve((region_length_joined * 2));
    data_ptr = data_copy;
    data_start = 0;
    data_end = region_length_joined - parameters->k_graph + 1;
  }

  for (uint64_t i = data_start; i < data_end; i++) {  // i+=parameters->kmer_step) {
    ProcessKmerCacheFriendly_((int8_t *) &(data_ptr[i]), i, local_score, mapping_data, index_read, index, index_secondary, read, parameters);
    mapping_data->iteration += 1;
  }

  // Ensures that the previous information stored in the vertices will not influence the current bin.
  mapping_data->iteration += parameters->num_links * 2;
  // This is used to avoid the overflow of the iteration counter. We use half the uint64_t range, just to make sure.
  if (mapping_data->iteration > ITERATION_RESET_LIMIT) {
//    for (int64_t current_vertex = 0; current_vertex < vertices_.num_vertices;
//        current_vertex++)
    mapping_data->vertices.EraseValues();
    mapping_data->iteration = 0;
  }

  if (data_copy)
    delete[] data_copy;
  data_copy = NULL;
  data_ptr = NULL;

  return 0;
}

int GraphMap::ProcessKmerCacheFriendly_(int8_t *kmer, int64_t kmer_start_position, ScoreRegistry *local_score, MappingData* mapping_data, Index *index_read, const Index* index, const Index* index_secondary, const SingleSequence* read, const ProgramParameters* parameters) {

  int64_t k = parameters->k_graph;
  int64_t num_links = parameters->num_links;

  uint64_t hits_start = -1;
  uint64_t num_hits = 0;
  int64_t *hits = NULL;
  int ret_search = ((IndexHash *) index_read)->FindAllRawPositionsOfIncrementalSeed(kmer, (uint64_t) k, (uint64_t) parameters->max_num_hits, &hits, &hits_start, &num_hits);

  if (ret_search == 1) {      // There are no hits for the current kmer.
    return 1;
  } else if (ret_search == 2) {  // There are too many seeds for the current kmer.
    // We ignore this case, and we'll use all the hits.
  } else if (ret_search > 2) {        // Something other went wrong.
    return 2;
  }

  // Sorting ensures the correct order of processing vertices.
  // Hits are sorted in descending order, because for each vertex, l previous vertices
  // are checked. Sorting solves the problem of repeats that occur very near each other
  // (for short kmers very often, i.e. any homopolymer run longer than k), and remove the
  // need for a backbufffer, thus reduced memory consumption and faster execution.
//  std::vector<int64_t> sorted_index_vector;
//  sorted_index_vector.resize(num_hits);
//  for (int64_t i=0; i<num_hits; i++) {
//    sorted_index_vector[i] = hits[i + hits_start];
//  }
//  std::sort(sorted_index_vector.begin(), sorted_index_vector.end(), std::greater<int64_t>());

  int64_t num_vertices = mapping_data->vertices.num_vertices;

  int64_t *hits_start_ptr = &hits[hits_start];

  for (int64_t i = 0; i < (num_hits); i++) {
    // Each hit position is a location on the read. Reference position is passed through function parameter kmer_start.
    int64_t hit = hits_start_ptr[i];
    int64_t position = num_vertices - hit - 1;
    int64_t best_vertex_idx = -1;

    // We have reversed the order of coordinates, so 'hit' variable tells us how far we are to the end of the vertex array.
    int64_t max_j = std::min(num_links, hit) + position;
    int64_t steps_away_query = 0;
    int64_t best_vertex_length = -1;
    for (int64_t j = (position + 1); j <= max_j; j++) {
      int64_t timestamp_diff = mapping_data->iteration - mapping_data->vertices.timestamps[j];
      // The timestamp_diff < mapping_data->iteration prevents negative mapping_data->vertices.num_kmers[j] values.
      if (timestamp_diff <= mapping_data->iteration && timestamp_diff <= num_links && timestamp_diff > 0) {
        int64_t vertex_length = mapping_data->vertices.num_kmers[j];
        if (vertex_length > 0 && (vertex_length > best_vertex_length)) {
          best_vertex_idx = j;
          best_vertex_length = vertex_length;
          steps_away_query = j - position;
        }
      }
    }

    if (best_vertex_idx == -1) {
      mapping_data->vertices.timestamps[position] = mapping_data->iteration;
      mapping_data->vertices.reference_starts[position] = kmer_start_position;
      mapping_data->vertices.reference_ends[position] = kmer_start_position + k;
      mapping_data->vertices.query_starts[position] = hit;
      mapping_data->vertices.query_ends[position] = hit + k;
      mapping_data->vertices.num_kmers[position] = 1;
      mapping_data->vertices.covered_bases_queries[position] = k;
      mapping_data->vertices.covered_bases_references[position] = k;
      mapping_data->vertices.registry_numbers[position] = -1;

    } else {
      mapping_data->vertices.reference_starts[position] = mapping_data->vertices.reference_starts[best_vertex_idx];
      mapping_data->vertices.query_starts[position] = mapping_data->vertices.query_starts[best_vertex_idx];
      mapping_data->vertices.registry_numbers[position] = mapping_data->vertices.registry_numbers[best_vertex_idx];
      mapping_data->vertices.covered_bases_queries[position] = mapping_data->vertices.covered_bases_queries[best_vertex_idx];
      mapping_data->vertices.covered_bases_references[position] = mapping_data->vertices.covered_bases_references[best_vertex_idx];
      mapping_data->vertices.num_kmers[position] = mapping_data->vertices.num_kmers[best_vertex_idx];

      mapping_data->vertices.timestamps[position] = mapping_data->iteration;
      mapping_data->vertices.reference_ends[position] = kmer_start_position + k;
      mapping_data->vertices.query_ends[position] = hit + k;
      mapping_data->vertices.num_kmers[position] += 1;
      int64_t steps_away_reference = kmer_start_position - (mapping_data->vertices.reference_ends[best_vertex_idx] - k);
      mapping_data->vertices.covered_bases_queries[position] += ((steps_away_query < k) ? steps_away_query : k);  // Check if hit overlaps the existing path, or is a little offset.
      mapping_data->vertices.covered_bases_references[position] += ((steps_away_reference < k) ? steps_away_reference : k);  // Check if hit overlaps the existing path, or is a little offset.

      // Put some constraints on the length of an anchor.
//      if (mapping_data->vertices.covered_bases_queries[position] > (parameters->k_graph + (parameters->num_links / parameters->k_graph)) &&
//          mapping_data->vertices.covered_bases_references[position] > (parameters->k_graph + (parameters->num_links / parameters->k_graph))) {
//      if (mapping_data->vertices.covered_bases_queries[position] >= (parameters->k_graph*3) &&
//          mapping_data->vertices.covered_bases_references[position] >= (parameters->k_graph*3)) {
      if (mapping_data->vertices.covered_bases_queries[position] >= parameters->min_num_anchor_bases &&
          mapping_data->vertices.covered_bases_references[position] >= parameters->min_num_anchor_bases) {
        local_score->Register(mapping_data->vertices, position);
      }
    }
  }

  return 0;
}
