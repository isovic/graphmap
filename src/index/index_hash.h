/*
 * index_hash.h
 *
 *  Created on: Oct 21, 2014
 *      Author: ivan
 */

#ifndef INDEX_HASH_H_
#define INDEX_HASH_H_

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <vector>
#include "index/index.h"
#include "sequences/single_sequence.h"

typedef std::vector<int64_t> HashIndexVector;
typedef std::vector<int64_t> SAIndexVector;

class IndexHash: public Index {
 public:
  IndexHash();
  ~IndexHash();

  void Init();
  void Clear();

  int GenerateFromSingleSequenceOnlyForward(const SingleSequence &sequence);
  int FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **hits, uint64_t *start_hit, uint64_t *num_hits) const;
  int FindAllRawPositionsOfSeedHash(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, HashIndexVector **hits, uint64_t *start_hit, uint64_t *num_hits, int64_t *kmer_hash_key);
  inline int64_t GenerateHashKey(const int8_t *seed, uint64_t seed_length) const;
  void CountKmers(const SingleSequence &sequence, int k, int64_t **ret_kmer_counts, int64_t *ret_num_kmers) const; // , int k, std::vector<int64_t> &ret_kmer_counts);
  inline int64_t UpdateHashKey(const int8_t *seed, uint64_t seed_length, int64_t hash_key) const;

  int LoadFromFile(std::string index_path);
  int GenerateFromFile(std::string sequence_file_path);
  int GenerateFromSequenceFile(SequenceFile &sequence_file);
  int GenerateFromSingleSequence(SingleSequence &sequence);     // TODO: Implement this method in other derived types (other than IndexSA).
  int LoadOrGenerate(std::string reference_path, std::string out_index_path, bool verbose=false);
  int StoreToFile(std::string output_index_path);

  int64_t RawPositionConverter(int64_t raw_position, int64_t query_length, int64_t *ret_absolute_position=NULL, int64_t *ret_relative_position=NULL, SeqOrientation *ret_orientation=NULL, int64_t *ret_reference_index_with_reverse=NULL) const;
  int64_t RawPositionToReferenceIndexWithReverse(int64_t raw_position) const;

  void Verbose(FILE *fp) const;
  std::string VerboseToString() const;

  int IsManualCleanupRequired(std::string function_name) const;

  int FindAllRawPositionsOfIncrementalSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **hits, uint64_t *start_hit, uint64_t *num_hits);



  const int8_t* get_data() const;
  uint64_t get_data_length() const;
  uint64_t get_data_length_forward() const;
  const std::vector<std::string>& get_headers() const;
  int get_k() const;
  void set_k(int k);

  uint64_t get_num_sequences() const;
  uint64_t get_num_sequences_forward() const;
  const std::vector<uint64_t>& get_reference_lengths() const;
  const std::vector<uint64_t>& get_reference_starting_pos() const;

 private:
  uint64_t num_sequences_;
  uint64_t data_length_;
  uint64_t num_sequences_forward_;
  uint64_t data_length_forward_;
  std::vector<uint64_t> reference_starting_pos_;
  std::vector<uint64_t> reference_lengths_;
  int8_t *data_;
  int64_t *suffix_array_;
  std::vector<std::string> headers_;

  uint64_t num_threads_;
  uint64_t data_ptr_;

//  std::vector<HashIndexVector> kmer_hash_;
  int64_t **kmer_hash_;
  int64_t kmer_hash_size_;
  int64_t *kmer_hash_counts_;
  int64_t kmer_hash_counts_size_;
  int64_t kmer_hash_last_key_;
  bool kmer_hash_last_key_initialized_;

  int64_t *all_kmers_;
  int64_t all_kmers_size_;

//  int64_t **kmer_hash1_;
//  int64_t kmer_hash_size_;
  int k_;

  int SerializeIndex_(FILE *fp_out);
  int DeserializeIndex_(FILE *fp_in);
  int CreateIndex_(int8_t *data, uint64_t data_length);
};

#endif /* INDEX_HASH_H_ */
