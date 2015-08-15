/*
 * index_spaced_hash.h
 *
 *  Created on: July 11, 2015
 *      Author: ivan
 */

#ifndef INDEX_SPACED_HASH_FAST_H_
#define INDEX_SPACED_HASH_FAST_H_

#include <vector>
#include "divsufsort64.h"
#include "index/index.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"

#define MASK_REF_ID       ((uint64_t) 0x00000000FFFFFFFF)
#define MASK_SEED_POS     ((uint64_t) 0xFFFFFFFF00000000)



class IndexSpacedHashFast : public Index {
 public:
  IndexSpacedHashFast();
  ~IndexSpacedHashFast();

  IndexSpacedHashFast(uint32_t shape_type);

  void Clear();
  int FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **entire_sa, uint64_t *start_hit, uint64_t *num_hits) const;

  int64_t GenerateHashKeyFromShape(int8_t *seed, const char *shape, int64_t shape_length) const;
  int64_t CalcNumHashKeysFromShape(const char *shape, int64_t shape_length) const;
  void CountKmersFromShape(int8_t *sequence_data, int64_t sequence_length, const char *shape, int64_t shape_length, int64_t **ret_kmer_counts, int64_t *ret_num_kmers) const;

  void Verbose(FILE *fp) const;
  std::string VerboseToString() const;

  int IsManualCleanupRequired(std::string function_name) const;

  int InitShapes(std::string shape_for_indexing, std::vector<std::string> &shapes_for_search);

  int FindAllRawPositionsOfSeedKey(int64_t hash_key, int64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *ret_start_hit, uint64_t *ret_num_hits) const;

  char* get_shape_index() const;
  void set_shape_index(char* shapeIndex);
  int64_t get_shape_index_length() const;
  void set_shape_index_length(int64_t shapeIndexLength);

  /// This class overrides the RawPositionToReferenceIndexWithReverse with a much faster implementation. In the IndexSpacedHashFast, the reference id is already stored in the seed hit position.
//  int64_t RawPositionToReferenceIndexWithReverse(int64_t raw_position) const;
//  int64_t RawPositionConverter(int64_t raw_position, int64_t query_length, int64_t *ret_absolute_position=NULL, int64_t *ret_relative_position=NULL, SeqOrientation *ret_orientation=NULL, int64_t *ret_reference_index_with_reverse=NULL) const;
  int64_t RawPositionConverter2(int64_t raw_position, int64_t query_length, int64_t *ret_absolute_position=NULL, int64_t *ret_relative_position=NULL, SeqOrientation *ret_orientation=NULL, int64_t *ret_reference_index_with_reverse=NULL) const;

//  int get_k() const;
//  void set_k(int k);
//  const std::vector<std::vector<int64_t> >& get_kmer_hash() const;
//  void set_kmer_hash(const std::vector<std::vector<int64_t> >& kmerHash);

 private:
//  std::vector<std::vector<int64_t> > kmer_hash_;
  int64_t **kmer_hash_array_;
  int64_t *kmer_counts_;
  int64_t num_kmers_;
//  int64_t k_;
  char *shape_index_;
  int64_t shape_index_length_;
  int64_t *all_kmers_;
  int64_t all_kmers_size_;
  std::vector<std::string> shapes_lookup_;

  int CreateIndex_(int8_t *data, uint64_t data_length);
  int SerializeIndex_(FILE *fp_out);
  int DeserializeIndex_(FILE *fp_in);

  int InitShapesPredefined(uint32_t shape_type);
};

#endif /* INDEX_SPACED_HASH_H_ */
