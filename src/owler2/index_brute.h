/*
 * index_brute.h
 *
 *  Created on: Mar 3, 2016
 *      Author: isovic
 */

#ifndef SRC_INDEX_BRUTE_H_
#define SRC_INDEX_BRUTE_H_

#include <vector>
#include <map>
#include <stdint.h>
// #include "index.h"
#include <stdint.h>
#include "sequences/sequence_file.h"
#include "utility/utility_general.h"
#include "compiled_shape.h"

int compare_int128(const void *a, const void *b);
#define MASK_LOWER_64 ((((unsigned __int128) 1) << 64) - 1)

#define islt(a,b) (a->key < b->key)



struct HashLocation {
  int32_t seq_id = 0;
  int32_t position = 0;
  int32_t key = 0;
  SeqOrientation orientation = kForward;
  char unused[3];
//  int32_t unused = 0; /// 4 bytes to fill the chunk to 128bits (power of two).

  bool operator<(const HashLocation& x) const {
    return this->key < x.key || (this->key == x.key && this->seq_id > x.seq_id);
//    return ((this->key == x.key) ? (this->seq_id > x.seq_id) : (this->key < x.key));
  }
  bool operator>(HashLocation& x) {
    return this->key > x.key;
  }
};

struct HashLocationBrute {
  int32_t seq_id = 0;
  float weight = 0.0f;
};

int compare_hash_pos(const void *a, const void *b);
int compare_int128(const void *a, const void *b);

//typedef std::vector<uint64_t> SeedListType;
#define SeedListType std::vector<uint64_t>

class IndexBrute {
 public:
  IndexBrute();
  ~IndexBrute();
  IndexBrute(std::string seq_path, std::vector<std::string> &shapes, float min_avg_seed_qv, bool index_reverse_strand, std::string index_desig, int64_t batch_size, int32_t num_threads);

  void Clear();

  int Create(std::string seq_path, std::vector<std::string> &shapes, float min_avg_seed_qv, bool index_reverse_strand, std::string index_desig, int64_t batch_size, int32_t num_threads);
  int CreateFromSequenceFile(const SequenceFile &seqs, float min_avg_seed_qv, const std::vector<CompiledShape> &compiled_shapes, bool index_reverse_strand, std::string index_desig, int32_t num_threads);
  int CreateFromSequenceFileSeparateReads(const SequenceFile &seqs, float min_avg_seed_qv, const std::vector<CompiledShape> &compiled_shapes, bool index_reverse_strand, std::string index_desig, int64_t batch_size, int32_t num_threads);
  int CreateFromSequenceFile(const SequenceFile &seqs, float min_avg_seed_qv, const std::vector<std::string> &shapes, bool index_reverse_strand, std::string index_desig, int32_t num_threads);
  int CreateFromSequence(SingleSequence &seq, std::vector<CompiledShape> &compiled_shapes, bool index_reverse_strand, std::string index_desig, int64_t batch_size, int32_t num_threads);

//  static std::vector<CompiledShape> CompileShapes(const std::vector<std::string> &shapes);

  /// Takes a buffer of bases (max 32 bases in 64-bits), 2bit packed, and extracts those inclusive ones (defined by a shape).
  /// Parameters:
  ///  @buffer an integer containing 2-bit packed values of the sequence, max. 32 bases from the starting position.
  ///  @shape a string specifying the shape to be extracted from the buffer. Specified with '1' as inclusive bases and '0' as don't care bases.
  static uint64_t CreateSeedFromShape(const CompiledShape &compiled_shape, uint64_t bases2bit);

  const std::vector<uint64_t *>& get_hash() const;
  void set_hash(const std::vector<uint64_t *>& hash);
  const std::vector<int64_t>& get_hash_counts() const;
  void set_hash_counts(const std::vector<int64_t>& hashCounts);
  const std::vector<uint64_t>& get_seed_list() const;
  void set_seed_list(const std::vector<uint64_t>& seedList);
  const std::vector<std::string>& get_headers() const;
  void set_headers(const std::vector<std::string>& headers);
  int64_t get_num_sequences() const;
  void set_num_sequences(int64_t numSequences);
  int64_t get_num_sequences_forward() const;
  void set_num_sequences_forward(int64_t numSequencesForward);
  const std::vector<int64_t>& get_reference_lengths() const;
  void set_reference_lengths(const std::vector<int64_t>& referenceLengths);
  const std::vector<float>& get_key_occ() const;
  void set_key_occ(const std::vector<float>& keyOcc);
  const std::map<int32_t, std::vector<int32_t> >& get_seeds_to_reads() const;
  void set_seeds_to_reads(const std::map<int32_t, std::vector<int32_t> >& seedsToReads);
  const std::vector<int64_t>& get_seq_seed_starts() const;
  void set_seq_seed_starts(const std::vector<int64_t>& seqSeedStarts);
  const std::vector<int64_t>& get_seq_seed_counts() const;
  void set_seq_seed_counts(const std::vector<int64_t>& seqSeedCounts);
  int64_t get_count_cutoff() const;
  void set_count_cutoff(int64_t countCutoff);
  const std::vector<std::vector<int8_t> >& get_data() const;
  void set_data(const std::vector<std::vector<int8_t> >& data);

 private:
//  uint64_t **hash_table_;
//  std::vector<std::vector<uint64_t> > hash_;
  std::vector<uint64_t> seed_list_;
  std::vector<uint64_t *> hash_;
  std::vector<int64_t> hash_counts_;
  std::vector<float> key_occ_;
  std::vector<std::vector<int8_t> > data_;

  int64_t num_sequences_;
  int64_t num_sequences_forward_;
  std::vector<int64_t> reference_lengths_;
  std::vector<std::string> headers_;
  // Seeds are extracted for each sequence separately, but are stored in a giant array. Each sequence 'i' is designated to belong to a part of that array, starting with seq_seed_starts_[i] position in seed_list_.
  std::vector<int64_t> seq_seed_starts_;
  std::vector<int64_t> seq_seed_counts_;

  int64_t count_cutoff_;

  int ProcessKmerSpectrum_(int8_t *seqdata, int8_t *seqqual, int64_t seqlen, float min_avg_seed_qv, uint64_t seq_id, const std::vector<CompiledShape> &compiled_shapes, uint64_t *seed_list);
  int ProcessKmerSpectrumSeparateReads_(int8_t *seqdata, int8_t *seqqual, int64_t seqlen, float min_avg_seed_qv, uint64_t seq_id, const std::vector<CompiledShape> &compiled_shapes, uint64_t *seed_list);

};

#endif /* SRC_INDEX_BRUTE_H_ */
