/*
 * process_read.cc
 *
 *  Created on: Jun 6, 2016
 *      Author: isovic
 */

#include "owler2/owler2.h"

uint128_t OWLER2_MASK_LOWER_64 = ((uint128_t) 0x0000000000000000FFFFFFFFFFFFFFFF);
uint128_t OWLER2_MASK_UPPER_64 = ((uint128_t) 0x0000000000000000FFFFFFFFFFFFFFFF) << 64;

int Owler2::ProcessRead(const ProgramParameters& parameters, const IndexGappedMinimizer* index, const SingleSequence* read, const std::vector<CompiledShape> &lookup_shapes, OwlerResult *owler_result) {
//  const auto &ls = index->get_lookup_shapes();

  std::vector<std::string> shvec = {"1111110111111"};
  std::vector<CompiledShape> ls = CompileShapes(shvec);

  for (int32_t i=0; i<ls.size(); i++) {
    printf ("Lookup shape: %s\n", ls[i].shape.c_str());
  }

  std::vector<uint128_t> minimizers;
  minimizers.resize(ls.size() * read->get_sequence_length() * 2);     // Allocate maximum space for the seeds. Factor 2 is for the reverse complement.
  int64_t num_minimizers = 0;
  IndexGappedMinimizer::CollectMinimizers(read->get_data(), read->get_quality(), read->get_sequence_length(), -1, 0, ls, 5, minimizers, &num_minimizers);

  printf ("num_minimizers = %ld\n", num_minimizers);

  int64_t max_count = ceil(index->get_count_cutoff());

  printf ("max_count = %ld\n", max_count);

//  SeedHashType hash_;
//  hash_.set_empty_key(empty_hash_key);
//  SeedHashValue a;
//  a.start = 123; a.num = 456;
//  hash_[1] = a;
//  a.start = 789; a.num = 012;
//  hash_[2] = a;
//  printf ("hash[1] = (%ld, %ld)\n", hash_[1].start, hash_[1].num);
//  printf ("hash[2] = (%ld, %ld)\n", hash_[2].start, hash_[2].num);
//  printf ("hash[3] = (%ld, %ld)\n", hash_[3].start, hash_[3].num);
//  exit(1);



  int64_t num_valid_hits = 0;
  std::vector<uint128_t> hits;
  hits.resize(read->get_sequence_length(), 0);
  int64_t last_hit_id = 0;

  for (int64_t i=0; i<num_minimizers; i++) {
    // Reverse complements are already indexed, so skip them here.
//    if (IS_HIT_REVERSE(minimizers[i])) { continue; }
//    printf ("\n");
//    printf ("%s\t%ld\t%ld\t%ld\n", SeedToString(((uint64_t) ((minimizers[i] >> 64) & 0x0FFFFFFFF)), 12).c_str(), ((uint64_t) ((minimizers[i] >> 32) & 0x0FFFFFFFF)), ((uint64_t) ((minimizers[i]) & 0x0FFFFFFFF)), GET_REAL_POS_FROM_HIT(((uint64_t) ((minimizers[i]) & 0x0FFFFFFFF))));
    if ((minimizers[i] & (((uint128_t) 1) << 31)) != 0) { continue; }

    uint64_t key = GET_KEY_FROM_HIT(minimizers[i]);
    uint64_t pos = GET_REAL_POS_FROM_HIT(minimizers[i]);
    uint128_t coded_pos = minimizers[i] & (0x0FFFFFFFFFFFFFFFF); // << 64;

//    printf ("[j = %d/%d] key = %s\tpos = %ld\n", (i+1), num_minimizers, SeedToString(key, 12).c_str(), pos);
//    fflush(stdout);

    uint128_t *seeds = NULL;
    int64_t num_seeds = 0;
    if (index->GetHits(key, &seeds, &num_seeds)) {
//      printf ("%s skip\n", SeedToString(key, 12).c_str());
      continue;
    }

//    printf ("%X\t%ld\n", key, num_seeds);
    if (max_count > 0 && num_seeds > max_count) {
//      printf ("%s\t%ld\tskip\n", SeedToString(key, 12).c_str(), num_seeds);
      continue;
    }

//    printf ("[%ld] %s\t%ld\t%ld\t%s\n", num_valid_hits, SeedToString(key, 12).c_str(), pos, num_seeds, SeedToString(key, 12).c_str());

    if ((last_hit_id + num_seeds) >= hits.size()) {
      hits.resize(hits.size() + read->get_sequence_length());
//      printf ("resize\n");
    }

//    memcpy(&hits[last_hit_id], seeds, num_seeds*sizeof(uint128_t));
    for (int32_t j=0; j<num_seeds; j++) {
//      hits[last_hit_id + j] = (hits[last_hit_id + j] & OWLER2_MASK_LOWER_64) | coded_pos;
      hits[last_hit_id + j] = (seeds[j] << 64) | coded_pos;
//      uint64_t t = hits[last_hit_id + j] >> 64;
//      printf ("[j = %d/%d] low 64 bits: pos = %ld, hits[%ld] = %ld, high 64 bits: hits[%ld] = %ld, seed pos: %ld\n", (j+1), num_seeds, pos, last_hit_id + j, GET_REAL_POS_FROM_HIT(hits[last_hit_id + j]), last_hit_id + j, GET_REAL_POS_FROM_HIT(hits[last_hit_id + j]>>64), GET_REAL_POS_FROM_HIT(seeds[j]));
//      printf ("%ld\n", GET_REAL_POS_FROM_HIT(t));
    }
    last_hit_id += num_seeds;
//    printf ("last_hit_id = %ld\n", last_hit_id);
//    printf ("\n");

    num_valid_hits += 1;
  }

//  for (int64_t i=0; i<last_hit_id; i++) {
//    uint64_t qid = GET_SEQ_ID_FROM_HIT(hits[i] >> 64);
//    uint64_t qstart = GET_REAL_POS_FROM_HIT((hits[i] >> 64));
//    uint64_t rid = GET_SEQ_ID_FROM_HIT(hits[i]);
//    uint64_t rstart = GET_REAL_POS_FROM_HIT(hits[i]);
//    printf ("qid = %ld, qstart = %ld, rid = %ld, rstart = %ld\n", qid, qstart, rid, rstart);
////    printf ("%ld\n", hits[i]);
//  }

//  printf ("last_hit_id = %ld\n", last_hit_id);
//  printf ("read_length = %ld\n", read->get_sequence_length());

  std::sort(hits.begin(), hits.begin() + last_hit_id);

  for (int64_t i=0; i<last_hit_id; i++) {
    uint64_t rid = GET_SEQ_ID_FROM_HIT(hits[i] >> 64);
//    uint64_t rstart = GET_REAL_POS_FROM_HIT((hits[i] >> 64));
    uint64_t rstart = GET_POS_FROM_HIT_WITH_REV((hits[i] >> 64));
    uint64_t qid = GET_SEQ_ID_FROM_HIT(hits[i]);
//    uint64_t qstart = GET_REAL_POS_FROM_HIT(hits[i]);
    uint64_t qstart = GET_POS_FROM_HIT_WITH_REV()(hits[i]);
    printf ("qid = %ld, qfwd = %d, qstart = %ld, qlen = %ld, rid = %ld, rfwd = %d, rstart = %ld, rlen = %ld\n", qid, IS_HIT_REVERSE(hits[i]), qstart, index->get_reference_lengths()[qid], rid, IS_HIT_REVERSE(hits[i] >> 64), rstart, index->get_reference_lengths()[rid]);
//    printf ("%ld\n", hits[i]);
  }

  exit(1);

  return 0;
}
