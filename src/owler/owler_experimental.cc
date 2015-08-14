/*
 * owler_experimental.cc
 *
 *  Created on: Jul 10, 2015
 *      Author: isovic
 */

#include "owler/owler.h"



std::string num_to_binary(uint64_t num) {
  std::stringstream ret;
//  std::string ret;

//  while (num) {
  for (int64_t i=0; i<sizeof(num)*8; i++) {
    if (num & 1)
      ret << "1";
    else
      ret << "0";
    num >>= 1;

    if ((i + 1) % 8 == 0 && ((i + 1) < (sizeof(num) * 8)))
      ret << "_";
  }

  std::string retstr = ret.str();
  std::string rev_ret(retstr.rbegin(), retstr.rend());
  return rev_ret;
}

//int Owler::CollectSeedHitsExperimental(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
//  bool test_verbose = false;
////  bool test_verbose = true;
//
//  if (owler_data == NULL || read == NULL || parameters == NULL)
//    return 1;
//  if (indexes.size() == 0 || (indexes.size() > 0 && indexes[0] == NULL))
//    return 2;
//
//  int64_t readlength = read->get_sequence_length();
//
//  /// Initialize the data structures to hold the results.
//  owler_data->Init((SingleSequence*) read, indexes);
////  SingleSequence read_2bit;
////  read_2bit.CopyFrom(*read);
////  read_2bit.ConvertDataFormat(kDataFormat2BitPacked);
//  SingleSequence *read_2bitpacked = read->ConvertDataFormatAndReturn(kDataFormat2BitPacked2);
//
////  printf ("Original:\n");
////  printf ("%s\n", read->GetSubstring(0, readlength).c_str());
////  printf ("\n");
////  printf ("Converted:\n");
////  printf ("%s\n", read_2bitpacked->GetSubstring(0, readlength).c_str());
////  printf ("\n");
////  printf ("In bits:\n");
////  for (int i=0; i<read_2bitpacked->get_data_length(); i++) {
////    printf ("%X ", (uint8_t) read_2bitpacked->get_data()[i]);
////  }
////  printf ("\n\n");
////  printf ("In numbers:\n");
////  for (int i=0; i<read_2bitpacked->get_data_length(); i++) {
////    printf ("%d ", (uint8_t) read_2bitpacked->get_data()[i]);
////  }
////  printf ("\n\n");
////  fflush(stdout);
////
//  uint64_t MASK_A = 0xFFF;
////  uint64_t MASK_B0 = 0x3F << 6;
////  uint64_t MASK_B1 = 0x3F;
//  uint64_t MASK_B[] = {MASK_A<<6*2, MASK_A<<((7*2)), MASK_A<<((8*2))};
//
//  if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//    printf ("MASK_A = %X\n", MASK_A);
//    fflush(stdout);
//  }
//
//  uint64_t seed_full = 0x0;
//
//  if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//    printf ("sizeof(seed_full) = %d\n", sizeof(seed_full));
//    fflush(stdout);
//  }
//
//  for (int64_t i = 0; i < readlength && readlength >= sizeof(seed_full); i += parameters->kmer_step) {
//    IndexSpacedHashFast *index = (IndexSpacedHashFast *) indexes[0];
//
//    /// Initialize the full seed with 8 bytes of the 2bit packed sequence (32 bases).
//    if (i == 0) {
//      for (int64_t j = 0; j < 8; j++) {
//        seed_full |= ((uint64_t) ((uint8_t) read_2bitpacked->get_data()[j])) << (j * 8);
//      }
//
//      if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//        printf ("Initialized:\n");
//        printf ("[%ld] seed_full = %lX\n", i, seed_full);
//        printf ("      seed_full = %s\n", num_to_binary(seed_full).c_str());
//        printf ("\n\n");
//        fflush(stdout);
//      }
//    }
//
////Cini mi se da je problem u tome sto kljuc negdje pogresno generiram, pa mi je dohvat krivi.
//
//    uint64_t seed_left_part = (seed_full & MASK_A);
//    uint64_t keys[] = {seed_left_part, seed_left_part, seed_left_part};
//    keys[0] = seed_left_part | ((uint64_t) (seed_full & MASK_B[0]));
//    keys[1] = seed_left_part | (((uint64_t) seed_full & MASK_B[1]) >> (1 * 2));
//    keys[2] = seed_left_part | (((uint64_t) seed_full & MASK_B[2]) >> (2 * 2));
//
//    if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//      uint64_t hash_key[] = {0, 0, 0, 0, 0};
//      hash_key[0] = index->GenerateHashKeyFromShape((int8_t *) &(read->get_data()[i]), "111111111111", 12);
//      hash_key[1] = index->GenerateHashKeyFromShape((int8_t *) &(read->get_data()[i]), "1111110111111", 13);
//      hash_key[2] = index->GenerateHashKeyFromShape((int8_t *) &(read->get_data()[i]), "11111100111111", 14);
//      hash_key[3] = index->GenerateHashKeyFromShape((int8_t *) &(read->get_data()[i]), "11111111111111", 14);
//      hash_key[4] = index->GenerateHashKeyFromShape((int8_t *) &(read->get_data()[i]), "11111111111111111111111111111111", 32);
//
//      printf ("[%ld] seed_full = %lX\t\t%lX\n", i, seed_full, hash_key[3]);
//      printf ("      seed_full = %s\n", num_to_binary(seed_full).c_str());
//      printf ("                  %s\n", num_to_binary(hash_key[3]).c_str());
//      printf ("                  %s\n", num_to_binary(hash_key[4]).c_str());
//      printf ("      keys[0] =   %X\t\t%X\n", keys[0], hash_key[0]);
//      printf ("      keys[1] =   %X\t\t%X\n", keys[1], hash_key[1]);
//      printf ("      keys[2] =   %X\t\t%X\n", keys[2], hash_key[2]);
//      printf ("\n");
//      printf ("Qry: %s\n", GetSubstring((char *) &(read->get_data()[i]), 16).c_str());
//      fflush(stdout);
//      if (i > 25)
//        break;
//    }
//
//    uint64_t hits_start = 0, num_hits = 0;
//    int64_t *hits = NULL;
//    for (int64_t key_id = 0; key_id < 3; key_id++) {
//      int ret_search = index->FindAllRawPositionsOfSeedKey(keys[key_id], 12, parameters->max_num_hits, &hits, &hits_start, &num_hits);
//      if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//  //      if (num_hits > 0) {
//          printf (", num_hits[%ld] = %ld (%X, %d) (seed_key = %ld)", key_id, num_hits, keys[key_id], ret_search, key_id);
//          fflush(stdout);
//  //        exit(1);
//  //      }
//      }
//
//      // Check if there is too many hits (or too few).
//      if (ret_search == 1) {
//        owler_data->num_seeds_with_no_hits += 1;
//      } else if (ret_search == 2) {
//        owler_data->num_seeds_over_limit += 1;
//      } else if (ret_search > 2) {
//        owler_data->num_seeds_errors += 1;
//      }
//
//      /// Counting kmers in regions of bin_size on the genome
//      for (int64_t j1 = hits_start; j1 < (hits_start + num_hits); j1++) {
//        int64_t position = hits[j1];
//
//        /// Find the index of the reference that was hit. This also includes the reverse sequences.
//        /// Reverse sequences are considered the same as any other reference sequence.
//        int64_t reference_index = position & 0x00000000FFFFFFFF;
//        int64_t position_local = position >> 32;
//
////        int64_t reference_index = index->RawPositionToReferenceIndexWithReverse(position);
//        int64_t reference_length = index->get_reference_lengths()[reference_index];
//        int64_t reference_start = index->get_reference_starting_pos()[reference_index];
//        int64_t reference_end = reference_start + reference_length;
////        int64_t position_local = position - reference_start;
//
////        if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
////          printf ("position = %ld\n", position);
////          printf ("position_local = %ld\n", position_local);
////          printf ("reference_index = %ld\n", reference_index);
////          printf ("is_reverse = %d\n", (int) (reference_index >= index->get_num_sequences_forward()));
////
////          printf ("\n");
////          printf ("Ref: %s\n", GetSubstring((char *) &(index->get_data()[position]), 16).c_str());
////          printf ("Qry: %s\n", GetSubstring((char *) &(read->get_data()[i]), 16).c_str());
////          printf ("\n");
////          printf ("Ref: %s %s  %s\n", GetSubstring((char *) (&index->get_data()[position]), 7-1).c_str(),
////                                     GetSubstring((char *) &(index->get_data()[position+7-1]), 1).c_str(),
////                                     GetSubstring((char *) &(index->get_data()[position+7]), 7-1).c_str());
////          printf ("Qry: %s    %s\n", GetSubstring((char *) (&(read->get_data()[i])), 7-1).c_str(),
////                                     GetSubstring((char *) &((read->get_data()[i+7-1])), 7-1).c_str());
////          printf ("Qry: %s %s  %s\n", GetSubstring((char *) (&(read->get_data()[i])), 7-1).c_str(),
////                                     GetSubstring((char *) &((read->get_data()[i+7-1])), 1).c_str(),
////                                     GetSubstring((char *) &((read->get_data()[i+7])), 7-1).c_str());
////          printf ("Qry: %s %s %s\n", GetSubstring((char *) (&(read->get_data()[i])), 7-1).c_str(),
////                                     GetSubstring((char *) &((read->get_data()[i+7-1])), 2).c_str(),
////                                     GetSubstring((char *) &((read->get_data()[i+8])), 7-1).c_str());
////          printf ("\n\n");
////
////  //        printf ("position = %ld\n", position);
////  //        fflush(stdout);
////        }
//
//        if (reference_index < 0) {
//          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index.\n"), "SelectRegionsWithHoughAndCircular");
//          continue;
//        }
//        /// Don't count self hits
//        if (index->get_headers()[reference_index % index->get_num_sequences_forward()] == ((std::string) read->get_header())) {
//          continue;
//        }
//        /// Count unique hits for a pair of reads.
//        if (owler_data->overlaps[reference_index].last_update < (i + 1)) {
//          owler_data->overlaps[reference_index].num_unique_hits += 1;
//        }
//
//        owler_data->overlaps[reference_index].last_update = (i + 1);
//        SeedHit seed_hit;
//        owler_data->overlaps[reference_index].seed_hits.push_back(SeedHit((uint32_t) i, (uint32_t) position_local, 0));
//      }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)
//    }
//
//    /// Shift the seed by one base, to prepare it for the next round.
//    seed_full = seed_full >> 2;
//    /// Check if a full byte is already removed, and reload.
////    if ((i + 7) < readlength && (i + 7) % 8 == 7) {
//    int64_t byte_index = (i + 1)/4 + 7;
//    if ((i + 1) % 4 == 0 && byte_index < read_2bitpacked->get_data_length()) {
//
//      if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//        printf ("\nAdding new byte to the end of seed_full.\n");
//        printf ("Qry: %s\n", GetSubstring((char *) &(read->get_data()[i]), 40).c_str());
//        printf ("Qry: %s\n", GetSubstring((char *) &(read->get_data()[i]), 32).c_str());
//        printf ("Qry: %s\n", GetSubstring((char *) &(read->get_data()[i + 4*7]), 4).c_str());
//        printf ("%X\n", read_2bitpacked->get_data()[byte_index]);
//        fflush(stdout);
//      }
//
//      seed_full |= (((uint64_t) read_2bitpacked->get_data()[byte_index]) << (7 * 8));
//    }
//
////    if (num_hits > 0) {
//    if (test_verbose == true && parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//      printf ("\n");
//      fflush(stdout);
//    }
////    }
//  }
//
//  if (read_2bitpacked)
//    delete read_2bitpacked;
//
//  return 0;
//
////      int64_t start_shift = i
//
//
//
////    int8_t *seed = (int8_t *) &(read->get_data()[i]);
////
////    int32_t seed_shape_id = 0;
////
////    /// Loop through all indexes to collect all different hits.
////    for (int64_t j = 0; j < indexes.size(); j++) {
////      seed_shape_id = j;
////
////      Index *index = indexes[j];
////      if (index == NULL)
////        continue;
////
////      int64_t k = (int64_t) ((IndexSpacedHash *) index)->get_shape_index_length();
////      if ((i + k) >= readlength)
////        continue;
////
////      uint64_t hits_start = 0, num_hits = 0;
////      int64_t *hits = NULL;
////
////      int ret_search = index->FindAllRawPositionsOfSeed(seed, k, parameters->max_num_hits, &hits, &hits_start, &num_hits);
////
////      // Check if there is too many hits (or too few).
////      if (ret_search == 1) {
////        owler_data->num_seeds_with_no_hits += 1;
////      } else if (ret_search == 2) {
////        owler_data->num_seeds_over_limit += 1;
////      } else if (ret_search > 2) {
////        owler_data->num_seeds_errors += 1;
////      }
////
////      /// Counting kmers in regions of bin_size on the genome
////      for (int64_t j1 = hits_start; j1 < (hits_start + num_hits); j1++) {
////        int64_t position = hits[j1];
////
////        /// Find the index of the reference that was hit. This also includes the reverse sequences.
////        /// Reverse sequences are considered the same as any other reference sequence.
////        int64_t reference_index = index->RawPositionToReferenceIndexWithReverse(position);
////        int64_t reference_length = index->get_reference_lengths()[reference_index];
////        int64_t reference_start = index->get_reference_starting_pos()[reference_index];
////        int64_t reference_end = reference_start + reference_length;
////        int64_t position_local = position - reference_start;
////
////        if (reference_index < 0) {
////          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index.\n"), "SelectRegionsWithHoughAndCircular");
////          continue;
////        }
////        /// Don't count self hits
////        if (index->get_headers()[reference_index % index->get_num_sequences_forward()] == ((std::string) read->get_header())) {
////          continue;
////        }
////        /// Count unique hits for a pair of reads.
////        if (owler_data->overlaps[reference_index].last_update < (i + 1)) {
////          owler_data->overlaps[reference_index].num_unique_hits += 1;
////        }
////
////        owler_data->overlaps[reference_index].last_update = (i + 1);
////        SeedHit seed_hit;
////        owler_data->overlaps[reference_index].seed_hits.push_back(SeedHit((uint32_t) i, (uint32_t) position_local, seed_shape_id));
////      }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)
////
////      if (index->IsManualCleanupRequired("FindAllRawPositionsOfSeed") == 0) {
////        if (hits)
////          free(hits);
////        hits = NULL;
////      }
////      hits = NULL;
////    }
////  }
////
////  return 0;
//}


int Owler::CollectSeedHitsExperimental(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
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
//    IndexSpacedHashFast *index = (IndexSpacedHashFast *) indexes[0];
    IndexOwler *index = (IndexOwler *) indexes[0];

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
          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index.\n"), "SelectRegionsWithHoughAndCircular");
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

int Owler::CollectSeedHitsExperimentalCalcSubseedsFast(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
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
//    IndexSpacedHashFast *index = (IndexSpacedHashFast *) indexes[0];
    IndexOwler *index = (IndexOwler *) indexes[0];

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
          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index.\n"), "SelectRegionsWithHoughAndCircular");
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

int Owler::CollectSeedHitsExperimentalSubseededIndex(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  bool test_verbose = false;

  if (owler_data == NULL || read == NULL || parameters == NULL)
    return 1;
  if (indexes.size() == 0 || (indexes.size() > 0 && indexes[0] == NULL))
    return 2;

  IndexOwler *index = (IndexOwler *) indexes[0];
  int64_t read_id = read->get_sequence_absolute_id();
  int64_t readlength = read->get_sequence_length();
  /// Initialize the data structures to hold the results.
  owler_data->Init((SingleSequence*) read, indexes);
  SubIndex *read_subindex = (SubIndex *) index->get_read_subindex()[read_id];

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
      int64_t position_local = position >> 32;

      /////////////////////
      ///// This handles self-overlapping, and only compares the read to uper-half of the matrix.
      /////////////////////
      if ((reference_index % index->get_num_sequences_forward()) <= read_id)
        continue;

      int64_t reference_length = index->get_reference_lengths()[reference_index];
      int64_t reference_start = index->get_reference_starting_pos()[reference_index];
      int64_t reference_end = reference_start + reference_length;

      if (reference_index < 0) {
        LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index.\n"), "SelectRegionsWithHoughAndCircular");
        continue;
      }
      /// Don't count self hits
      if (index->get_headers()[reference_index % index->get_num_sequences_forward()] == ((std::string) read->get_header())) {
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
