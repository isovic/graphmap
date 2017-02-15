/*
 * owler_experimental.cc
 *
 *  Created on: Jul 10, 2015
 *      Author: isovic
 */

#include "owler/owler.h"

int Owler::CollectSeedHitsExperimental(OwlerData* owler_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  bool test_verbose = false;

  if (owler_data == NULL || read == NULL || parameters == NULL)
    return 1;
  if (indexes.size() == 0 || (indexes.size() > 0 && indexes[0] == NULL))
    return 2;

  int64_t readlength = read->get_sequence_length();

  /// Initialize the data structures to hold the results.
  owler_data->Init((SingleSequence*) read, indexes);

  SingleSequence *read_2bitpacked = read->ConvertDataFormatAndReturn(kDataFormat2BitPacked2);

  uint64_t MASK_A = 0xFFF;
  uint64_t MASK_B[] = {MASK_A<<6*2, MASK_A<<((7*2)), MASK_A<<((8*2))};

  if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    printf ("MASK_A = %X\n", MASK_A);
    fflush(stdout);
  }

  uint64_t seed_full = 0x0;

  if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    printf ("sizeof(seed_full) = %d\n", sizeof(seed_full));
    fflush(stdout);
  }

  for (int64_t i = 0; i < readlength && readlength >= sizeof(seed_full); i += parameters->kmer_step) {
    auto index = indexes[0];

    /// Initialize the full seed with 8 bytes of the 2bit packed sequence (32 bases).
    if (i == 0) {
      for (int64_t j = 0; j < 8; j++) {
        seed_full |= ((uint64_t) ((uint8_t) read_2bitpacked->get_data()[j])) << (j * 8);
      }
    }

    uint64_t seed_left_part = (seed_full & MASK_A);
    uint64_t keys[] = {seed_left_part, seed_left_part, seed_left_part};
    keys[0] = seed_left_part | ((uint64_t) (seed_full & MASK_B[0]));
    keys[1] = seed_left_part | (((uint64_t) seed_full & MASK_B[1]) >> (1 * 2));
    keys[2] = seed_left_part | (((uint64_t) seed_full & MASK_B[2]) >> (2 * 2));

    uint64_t hits_start = 0, num_hits = 0;
    int64_t *hits = NULL;
    for (int64_t key_id = 0; key_id < 3; key_id++) {
      int ret_search = index->FindAllRawPositionsOfSeedKey(keys[key_id], 12, parameters->max_num_hits, &hits, &hits_start, &num_hits);

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
        int64_t reference_index = position & 0x00000000FFFFFFFF;
        int64_t position_local = position >> 32;

        int64_t reference_length = index->get_reference_lengths()[reference_index];
        int64_t reference_start = index->get_reference_starting_pos()[reference_index];
        int64_t reference_end = reference_start + reference_length;

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
        owler_data->overlaps[reference_index].seed_hits.push_back(SeedHit((uint32_t) i, (uint32_t) position_local, 0));
      }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)
    }

    /// Shift the seed by one base, to prepare it for the next round.
    seed_full = seed_full >> 2;
    /// Check if a full byte is already removed, and reload.
    int64_t byte_index = (i + 1)/4 + 7;
    if ((i + 1) % 4 == 0 && byte_index < read_2bitpacked->get_data_length()) {

      seed_full |= (((uint64_t) read_2bitpacked->get_data()[byte_index]) << (7 * 8));
    }
  }

  if (read_2bitpacked)
    delete read_2bitpacked;

  return 0;
}

int Owler::CollectSeedHitsExperimentalCalcSubseedsFast(OwlerData* owler_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  bool test_verbose = false;

  if (owler_data == NULL || read == NULL || parameters == NULL)
    return 1;
  if (indexes.size() == 0 || (indexes.size() > 0 && indexes[0] == NULL))
    return 2;

  int64_t readlength = read->get_sequence_length();

  /// Initialize the data structures to hold the results.
  owler_data->Init((SingleSequence*) read, indexes);

  SingleSequence *read_2bitpacked = read->ConvertDataFormatAndReturn(kDataFormat2BitPacked2);

  int64_t SHAPE_LENGTH = 12;
  uint64_t MASK_A = 0xFFF;
  uint64_t MASK_B[] = {MASK_A<<6*2, MASK_A<<((7*2)), MASK_A<<((8*2))};

  if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    printf ("MASK_A = %X\n", MASK_A);
    fflush(stdout);
  }

  uint64_t seed_full = 0x0;

  if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    printf ("sizeof(seed_full) = %d\n", sizeof(seed_full));
    fflush(stdout);
  }

//  printf ("%s\n", read->get_data());
//  fflush(stdout);

  for (int64_t i = 0; i < (readlength - SHAPE_LENGTH) && readlength >= sizeof(seed_full); i += parameters->kmer_step) {
    auto index = indexes[0];

    /// Initialize the full seed with 8 bytes of the 2bit packed sequence (32 bases).
    if (i == 0) {
      for (int64_t j = 0; j < 8; j++) {
        seed_full |= ((uint64_t) ((uint8_t) read_2bitpacked->get_data()[j])) << (j * 8);
      }
    }

    uint64_t seed_left_part = (seed_full & MASK_A);
    uint64_t keys[] = {seed_left_part, seed_left_part, seed_left_part};
    keys[0] = seed_left_part | ((uint64_t) (seed_full & MASK_B[0]));
    keys[1] = seed_left_part | (((uint64_t) seed_full & MASK_B[1]) >> (1 * 2));
    keys[2] = seed_left_part | (((uint64_t) seed_full & MASK_B[2]) >> (2 * 2));

    uint64_t hits_start = 0, num_hits = 0;
    int64_t *hits = NULL;
    for (int64_t key_id = 0; key_id < 3; key_id++) {
      int ret_search = index->FindAllRawPositionsOfSeedKey(keys[key_id], 12, parameters->max_num_hits, &hits, &hits_start, &num_hits);

      // Check if there is too many hits (or too few).
      if (ret_search == 1) {
        owler_data->num_seeds_with_no_hits += 1;
      } else if (ret_search == 2) {
        owler_data->num_seeds_over_limit += 1;
//        printf ("Tu sam 1!\n");
//        printf ("parameters->max_num_hits = %ld\n", parameters->max_num_hits);
//        fflush(stdout);
        continue;
      } else if (ret_search > 2) {
        owler_data->num_seeds_errors += 1;
      }

      /// Counting kmers in regions of bin_size on the genome
      for (int64_t j1 = hits_start; j1 < (hits_start + num_hits); j1++) {
        int64_t position = hits[j1];

        /// Find the index of the reference that was hit. This also includes the reverse sequences.
        /// Reverse sequences are considered the same as any other reference sequence.
        int64_t reference_index = ((uint64_t) position) & 0x00000000FFFFFFFF;
        int64_t position_local = (int64_t) (((uint64_t) position) >> 32);

        /////////////////////
        ///// This would handle self-overlapping, and only compares the read to uper-half of the matrix.
        ///// Not used in this function because for self-overlapping a faster function can be used (check 'CollectSeedHitsExperimentalSubseededIndex').
        /////////////////////
//        if ((reference_index % index->get_num_sequences_forward()) <= read_id)
//          continue;

        int64_t reference_length = index->get_reference_lengths()[reference_index];
        int64_t reference_start = index->get_reference_starting_pos()[reference_index];
        int64_t reference_end = reference_start + reference_length;

        if (reference_index < 0) {
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index.\n"), "SelectRegionsWithHoughAndCircular");
//          continue;
        }
        /// Don't count self hits
        if (index->get_headers()[reference_index % index->get_num_sequences_forward()] == ((std::string) read->get_header())) {
//          continue;
        }
        /// Count unique hits for a pair of reads.
        if (owler_data->last_update[reference_index] < (i + 1)) {
          owler_data->num_unique_hits[reference_index] += 1;
        }

        owler_data->last_update[reference_index] = (i + 1);
        SeedHit seed_hit;
//        owler_data->overlaps[reference_index].seed_hits.push_back(SeedHit((uint32_t) i, (uint32_t) position_local, 0));
        owler_data->seed_hits2.push_back(SeedHit2((uint32_t) i, (uint32_t) position_local, reference_index));
//        owler_data->seed_hits2.push_back(SeedHit2((uint32_t) query_pos, (uint32_t) position_local, reference_index));

      }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)
    }

    /// Shift the seed by one base, to prepare it for the next round.
    seed_full = seed_full >> 2;
    /// Check if a full byte is already removed, and reload.
    int64_t byte_index = (i + 1)/4 + 7;
    if ((i + 1) % 4 == 0 && byte_index < read_2bitpacked->get_data_length()) {

      seed_full |= (((uint64_t) read_2bitpacked->get_data()[byte_index]) << (7 * 8));
    }
  }

  if (read_2bitpacked)
    delete read_2bitpacked;

  return 0;
}

int Owler::CollectSeedHitsExperimentalSubseededIndex(OwlerData* owler_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  bool test_verbose = false;

  if (owler_data == NULL || read == NULL || parameters == NULL)
    return 1;
  if (indexes.size() == 0 || (indexes.size() > 0 && indexes[0] == NULL))
    return 2;

  auto index = indexes[0];
  int64_t read_id = read->get_sequence_absolute_id();
  int64_t readlength = read->get_sequence_length();
  /// Initialize the data structures to hold the results.
  owler_data->Init((SingleSequence*) read, indexes);
  SubIndex *read_subindex = (SubIndex *) index->get_read_subindex()[read_id];

//  for (int64_t i = 0; i < index->get_subindex_counts()[read_id]; i++) {
//    int64_t key = read_subindex[i].key;
//    uint32_t query_pos = read_subindex[i].position;
//    printf ("key = %ld,\tquery_pos = %d\n", key, query_pos);
//  }
//  fflush(stdout);

  for (int64_t i = 0; i < index->get_subindex_counts()[read_id]; i++) {
    int64_t key = read_subindex[i].key;
    uint32_t query_pos = read_subindex[i].position;
    uint64_t hits_start = 0, num_hits = 0;
    int64_t *hits = NULL;

    int ret_search = index->FindAllRawPositionsOfSeedKey(key, 12, parameters->max_num_hits, &hits, &hits_start, &num_hits);

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
      int64_t reference_index = position & 0x00000000FFFFFFFF;
      int64_t position_local = (((uint64_t) position) >> 32);
      int64_t reference_index_fwd = reference_index % index->get_num_sequences_forward();

      /////////////////////
      ///// This handles self-overlapping, and only compares the read to uper-half of the matrix.
      /////////////////////
      if ((reference_index_fwd) <= read_id)
        continue;

      int64_t reference_length = index->get_reference_lengths()[reference_index];
      int64_t reference_start = index->get_reference_starting_pos()[reference_index];
      int64_t reference_end = reference_start + reference_length;

      if (reference_index < 0) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index.\n"), "SelectRegionsWithHoughAndCircular");
        continue;
      }
      /// Don't count self hits
      if (index->get_headers()[reference_index_fwd] == ((std::string) read->get_header())) {
        continue;
      }
      /// Count unique hits for a pair of reads.
      if (owler_data->last_update[reference_index] < (i + 1)) {
        owler_data->num_unique_hits[reference_index] += 1;
      }

      owler_data->last_update[reference_index] = (i + 1);
      SeedHit seed_hit;
      owler_data->seed_hits2.push_back(SeedHit2((uint32_t) query_pos, (uint32_t) position_local, reference_index));
    }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)
  }

  return 0;
}
