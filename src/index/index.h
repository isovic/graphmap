/*
 * index.h
 *
 *  Created on: 23 May, 2014
 *      Author: Ivan Sovic
 */

#ifndef INDEX_H_
#define INDEX_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include "sequences/sequence_file.h"
#include "utility/utility_general.h"


#define INDEX_VERSION     ((int64_t) -3)

//class IndexChunk {
// public:
//  std::string designator;
//  int64_t dlen;
//  int8_t *chunk_data;
//
//  IndexChunk() {
//    designator = std::string("");
//    dlen = 0;
//    chunk_data = NULL;
//  }
//
//  ~IndexChunk() {
//    if (chunk_data)
//      free(chunk_data);
//    chunk_data = NULL;
//  }
//
//  void Clear() {
//    designator = std::string("");
//    dlen = 0;
//    if (chunk_data)
//      free(chunk_data);
//    chunk_data = NULL;
//  }
//
//  static int GetChunk(FILE *fp_in, IndexChunk &ret_chunk) {
//    ret_chunk.Clear();
//
//    char chunk_designator[9];
//    if (fread(chunk_designator, sizeof(char), 8, fp_in) != 8) {
//      return 1;
//    }
//    chunk_designator[8] = '\0';
//    ret_chunk.designator = std::string(chunk_designator);
//
//    if (fread(&ret_chunk.dlen, sizeof(ret_chunk.dlen), 1, fp_in) != 1) {
//      return 2;
//    }
//
//    ret_chunk.chunk_data = (int8_t *) calloc(ret_chunk.dlen, sizeof(int8_t));
//    if (fread(&ret_chunk.chunk_data, sizeof(int8_t), ret_chunk.dlen, fp_in) != 1) {
//      return 3;
//    }
//
//    return 0;
//  }
//
//  static int WriteChunk(FILE *fp_out, std::string out_designator, int64_t out_dlen, int8_t *out_data) {
//    if (out_designator.size() > 8) {
//      out_designator = out_designator.substr(0, 8);
//    } else if (out_designator.size() < 8) {
//      int num_empty_chars = 8 - out_designator.size();
//      out_designator += std::string(num_empty_chars, ' ');
//    }
//
//    fwrite(((int8_t *) out_designator.c_str()), sizeof(int8_t), 8, fp_out);
//    fwrite(&out_dlen, sizeof(out_dlen), 1, fp_out);
//    fwrite(out_data, sizeof(out_data), out_dlen, fp_out);
//
//    return 0;
//  }
//};

class Index {
 public:
  Index();
  virtual ~Index();

  virtual void Clear() = 0;

  virtual int LoadFromFile(std::string index_path);
  virtual int GenerateFromFile(std::string sequence_file_path);
  virtual int GenerateFromSequenceFile(const SequenceFile &sequence_file);
  virtual int GenerateFromSingleSequence(const SingleSequence &sequence);
  virtual int GenerateFromSingleSequenceOnlyForward(const SingleSequence &sequence);
  virtual int LoadOrGenerate(std::string reference_path, std::string out_index_path, bool verbose=false);
  virtual int StoreToFile(std::string output_index_path);

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
  virtual int64_t RawPositionConverter(int64_t raw_position, int64_t query_length, int64_t *ret_absolute_position=NULL, int64_t *ret_relative_position=NULL, SeqOrientation *ret_orientation=NULL, int64_t *ret_reference_index_with_reverse=NULL) const;
  virtual int64_t RawPositionToReferenceIndexWithReverse(int64_t raw_position) const;

  virtual int FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *start_hit, uint64_t *num_hits) const = 0;

  // Some implementations of the class might require manuall deallocation of memory (e.g. ret_hits from the FindAllRawPositionsOfSeed).
  // Returns 1 if no allocation is needed (C-style).
  virtual int IsManualCleanupRequired(std::string function_name) const = 0;



  virtual const int8_t* get_data() const;
  virtual uint64_t get_data_length() const;
  virtual uint64_t get_data_length_forward() const;
  virtual const std::vector<std::string>& get_headers() const;
  virtual uint64_t get_num_sequences_forward() const;
  virtual const std::vector<uint64_t>& get_reference_lengths() const;
  virtual const std::vector<uint64_t>& get_reference_starting_pos() const;

  virtual void Verbose(FILE *fp) const = 0;
  virtual std::string VerboseToString() const = 0;

 protected:
  uint64_t num_sequences_;
  uint64_t data_length_;
  uint64_t num_sequences_forward_;
  uint64_t data_length_forward_;
  std::vector<uint64_t> reference_starting_pos_;
  std::vector<uint64_t> reference_lengths_;
  int8_t *data_;
  std::vector<std::string> headers_;
  uint64_t data_ptr_;

 private:
  virtual int Serialize_(FILE *fp_out);
  virtual int Deserialize_(FILE *fp_in);
  virtual int DeprecatedDeserialize_(FILE *fp_in);
  virtual int SerializeIndex_(FILE *fp_out) = 0;
  virtual int DeserializeIndex_(FILE *fp_in) = 0;
  virtual int CreateIndex_(int8_t *data, uint64_t data_length) = 0;

  virtual int InsertHeaders_(const SequenceFile& sequence_file);
  virtual int InsertSequencesIntoData_(const SequenceFile& sequence_file);
  virtual int InsertReverseSequencesIntoData_(const SequenceFile& sequence_file);
  virtual int InsertSingleHeader_(const SingleSequence *sequence);
  virtual int InsertSingleSequenceIntoData_(const SingleSequence *sequence);
  virtual int InsertReverseSingleSequenceIntoData_(const SingleSequence *sequence);
  virtual int InsertHeader_(const char *header, uint64_t header_length);
  virtual int InsertSequence_(const int8_t *sequence_data, uint64_t sequence_length);
};

#endif /* INDEX_H_ */
