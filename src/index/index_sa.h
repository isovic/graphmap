/*
 * index_sa.h
 *
 *  Created on: Jun 4, 2014
 *      Author: ivan
 */

#ifndef INDEX_SA_H_
#define INDEX_SA_H_

#include <vector>
#include "libs/libdivsufsort-2.0.1-64bit/divsufsort64.h"
#include "index/index.h"
#include "log_system/log_system.h"

class IndexSA : public Index {
 public:
  IndexSA();
  ~IndexSA();

  void Clear();
  int FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **entire_sa, uint64_t *start_hit, uint64_t *num_hits) const; //, std::vector<int64_t> &return_positions);

  // Converts the raw position of a query to the real position on the original sequence. This is required in cases when the index has been
  // constructed from both the forward and the reverse complement sequences. Since the sequences are truncated into a single data array,
  // this function can be used to identify to which sequence the raw position belongs to, and whether it is the forward of the reverse strand.
  // @param raw_position the raw position returned by the suffix array search function.
  // @param query_length the length of the original query that was found on the raw position (e.g. seed length). In case the hit was on the reverse strand, the query_length is reduced from the relative hit position to point to the correct beginning of the hit.
  // @param ret_relative_position the position of the query within the bounds of the corresponding reference sequence (in range of [0, sequence_length> ).
  // @param ret_absolute_position the absolute position of the query on the suffix array data (forward and reverse sequences included).
  // @param ret_orientation kForward if the hit landed on the forward strand, kReverse otherwise.
  // @param ret_reference_index_with_reverse this is only for testing purposes. The meaning of this value is the following. The index is created from the original sequences (forward), followed with their reverse. This parameter does not return the same value necessarily as the ret_reference_index parameter, but the absolute index value. For example, if there is 5 reference sequences, then there are 10 sequences in the index (forward + reverse). Parameter ret_reference_index returns a value in range [0, 5> every time, while this parameter returns a value between [0, 10>.
  // @return returns a value lesser than 0 if something went wrong, otherwise it returns the index of the sequence to which the hit belongs to.
//  int64_t RawPositionConverter(saidx64_t raw_position, int64_t query_length, int64_t *ret_absolute_position=NULL, int64_t *ret_relative_position=NULL, SeqOrientation *ret_orientation=NULL, int64_t *ret_reference_index_with_reverse=NULL);
//  int64_t RawPositionToReferenceIndexWithReverse(saidx64_t raw_position);

  void Verbose(FILE *fp) const;
  std::string VerboseToString() const;

  int IsManualCleanupRequired(std::string function_name) const;

 private:
  int64_t *suffix_array_;

  int CreateIndex_(int8_t *data, uint64_t data_length);
  int SerializeIndex_(FILE *fp_out);
  int DeserializeIndex_(FILE *fp_in);
};

#endif /* INDEX_SA_H_ */
