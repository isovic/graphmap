/*
 * index_spaced_hash.h
 *
 *  Created on: Feb 5, 2015
 *      Author: ivan
 */

#ifndef INDEX_SPACED_HASH_H_
#define INDEX_SPACED_HASH_H_

#include <vector>
#include "index/index.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"



class IndexSpacedHash : public Index {
 public:
  IndexSpacedHash();
  ~IndexSpacedHash();

  IndexSpacedHash(uint32_t shape_type);

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
