/*
 * index_spaced_hash_fast.h
 *
 *  Created on: July 11, 2015
 *      Author: ivan
 */

#ifndef INDEX_SPACED_HASH_FAST_H_
#define INDEX_SPACED_HASH_FAST_H_

#include <vector>
#include <algorithm>
#include <map>
#include "index/index.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"

#define MASK_REF_ID       ((uint64_t) 0x00000000FFFFFFFF)
#define MASK_SEED_POS     ((uint64_t) 0xFFFFFFFF00000000)
#define MASK_32_BIT       ((uint64_t) 0x00000000FFFFFFFF)



struct SeedHit3 {
  int32_t ref_id;
  int32_t y;
  int32_t x;
};
struct seed_hit3_compare
{
    inline bool operator() (const SeedHit3& op1, const SeedHit3& op2) {
      if (op1.ref_id != op2.ref_id)
        return (op1.ref_id < op2.ref_id);
      else if (op1.y != op2.y)
        return (op1.y < op2.y);
      else
        return (op1.x < op2.x);
    }
};

class MaskOperation {
 public:
  uint64_t mask;
  int32_t shift;

  MaskOperation(uint64_t m_mask, int32_t m_shift) {
    mask = m_mask;
    shift = m_shift;
  }
};

class CompiledSeed {
 public:
  std::vector<MaskOperation> mask_ops;
  std::string shape;

  CompiledSeed(std::string m_shape);
  void Generate(std::string m_shape);
  uint64_t Apply(uint64_t full_seed);
};

template <typename T>
void OrderedSortArray(const T* values, int64_t values_size, std::vector<size_t> &indices) {
    indices.resize(values_size);
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
}

class IndexSpacedHashFast : public Index {
 public:
  IndexSpacedHashFast();
  ~IndexSpacedHashFast();

  IndexSpacedHashFast(uint32_t shape_type);

  void Clear();
  int LoadOrGenerateTranscriptome(std::string reference_path, std::string gtf_path, std::string out_index_path, bool verbose);
  int GenerateTranscriptomeFromFile(std::string sequence_file_path, std::string gtf_path);

  int FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **entire_sa, uint64_t *start_hit, uint64_t *num_hits) const;
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

  int CalcAllKeysFromSequence(const SingleSequence *read, int64_t kmer_step, std::vector<int64_t> &ret_hash_keys, std::vector<int64_t> &ret_key_counts);
  int LookUpHashKeys(int64_t bin_size, const SingleSequence *read, const std::vector<int64_t> &hash_keys, const std::vector<int64_t> &key_counts, std::vector<SeedHit3> &ret_hits);
  void CalcPercentileHits(double percentile, int64_t *ret_count, int64_t *ret_max_seed_count=NULL);

  // Experimental function, does not copy the hits but only returns the pointers to the buckets.
  int FindAllRawPositionsOfSeedNoCopy(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, std::vector<int64_t *> &ret_hits, std::vector<uint64_t> &ret_num_hits) const;

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

  std::vector<CompiledSeed> compiled_seeds_;

  int CreateIndex_(int8_t *data, uint64_t data_length);
  int SerializeIndex_(FILE *fp_out);
  int DeserializeIndex_(FILE *fp_in);

  int InitShapesPredefined(uint32_t shape_type);
  int64_t GenerateHashKeyFromShape(int8_t *seed, const char *shape, int64_t shape_length) const;
  int64_t CalcNumHashKeysFromShape(const char *shape, int64_t shape_length) const;
  void CountKmersFromShape(int8_t *sequence_data, int64_t sequence_length, const char *shape, int64_t shape_length, int64_t **ret_kmer_counts, int64_t *ret_num_kmers) const;
  void CalcPercentileHits_(int64_t *seed_counts, int64_t num_seeds, double percentile, int64_t *ret_count, int64_t *ret_max_seed_count=NULL);

  // Creates a transcriptome from a given reference sequence and a path to a file with gene annotations.
  // Parameters:
  // @annotations_path Path to a GFF file (or another supported format) which contains the annotations of exonic regions.
  // @references A SequenceFile object which contains reference sequences already loaded from disk.
  // @transcripts A SequenceFile which will contain the generated transcriptomes.
  // @return 0 if everything went fine (C-style).
  int MakeTranscript_(std::string annotations_path, const SequenceFile &references, SequenceFile &transcripts);
  int GenerateExons(std::string annotations_path, std::map<std::string, std::vector<std::pair<std::string, char>>> &seqToTrans, std::map<std::string, std::vector<std::pair<int64_t, int64_t>>> &transToExons);
  std::string trim(std::string s);
  std::vector<std::string> split(std::string s, char c);
  std::string getSequenceName(SingleSequence &seq);
  std::string getTID(std::string s);
  void outputSeq(char *header, size_t headerLen, const int8_t *seq, size_t seqLen);
};

#endif /* INDEX_SPACED_HASH_H_ */
