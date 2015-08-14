/*
 * Copyright 2014, Ivan Sovic.
 * All rights reserved.
 *
 * SequenceFile.h
 *
 *  Created on: 15 May, 2014
 *      Author: Ivan Sovic
 */

#ifndef SEQUENCEFILE_H_
#define SEQUENCEFILE_H_

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <zlib.h>
//#include "bwa/bwa.h"
//#include "bwa/bwamem.h"
//#include "bwa/kvec.h"
//#include "sequences/utils.h"
#include "sequences/single_sequence.h"
#include "utility/utility_conversion-inl.h"

#include "sequences/kseq.h"
//KSEQ_DECLARE(gzFile)
KSEQ_INIT(gzFile, gzread)

// Type defined to hold sequences in the SequenceFile container.
typedef std::vector<SingleSequence *> SequenceVector;

// This is used for parsing and storing sequences. Supports
// parsing from FASTA and FASTQ files. Files can be parsed
// entirely, or in batches, where the size of a batch can
// be limited either by the number of sequences to load, or
// by the size in MegaBytes to occupy in memory.
// Sample usage 1:
//    SequenceFile sequences1("test1.fasta");
//    sequences1.Verbose(stdout);
//
// Sample usage 2:
//    SequenceFile sequences2;
//    sequences2.OpenFileForBatchLoading("test2.fasta");
//    sequences2.LoadNextBatchNSequences(15);
//    sequences2.LoadNextBatchInMegabytes(2048);
//    sequences2.CloseFileAfterBatchLoading();
//
class SequenceFile {
 public:
  SequenceFile();

  // Initializes the object, and calls LoadAllFromFastaOrFastq.
  SequenceFile(std::string file_path);

  // Initializes the object, and calls OpenFileForBatchLoading and
  // LoadNextBatchNSequences respectively.
  SequenceFile(std::string file_path, uint64_t num_seqs_to_load);

  ~SequenceFile();

  // Sets all member variables to initial values, and frees any memory
  // if it was allocated, followed by setting all pointers to NULL value.
  void Clear();

  // Clears only data related to sequences (concretely, the sequences_
  // vector and the current_data_size_ value).
  void ClearOnlyData();

  // Adds a new sequence to the container.
  // Inputs:
  //    sequence  - Pointer to a SingleSequence object holding a new sequence.
  void AddSequence(SingleSequence *sequence);

  // Given the path to a FASTA or a FASTQ file, this function opens the file
  // handlers to enable loading of data from the file.
  // Input:
  //    file_path - Path to the FASTA/FASTQ file.
  // Return:
  //    Returns 0 if successful.
  int OpenFileForBatchLoading(std::string file_path);

  // This function closes the file handles opened with OpenFileForBatchLoading
  // function. Must always be used in pair with OpenFileForBatchLoading.
  // Return:
  //    Returns 0 if successful.
  int CloseFileAfterBatchLoading();

  // Given the path to a FASTA or a FASTQ file, this function loads all
  // sequences present in the file into this object.
  // Input:
  //    file_path - Path to the FASTA/FASTQ file.
  // Return:
  //    Returns 0 if successful.
  int LoadAllFromFastaOrFastq(std::string file_path, bool randomize_non_acgt_bases=false);
  int LoadAllFromFastaOrFastqAsBatch(bool randomize_non_acgt_bases=false);

  // This function loads only a part of the sequences present in the
  // file into this object. Before using LoadNextBatchNSequences a file
  // needs to explicitly be opened with the OpenFileForBatchLoading function.
  // This function loads N sequences at once. When called the next time,
  // it will load another batch of N sequences, until EOF is reached.
  // After all batches have been loaded from file, users must call
  // the CloseFileAfterBatchLoading function themselves.
  // Input:
  //    num_seqs_to_load  - Number of sequences to load in a single batch.
  //                        If bigger than total number of sequences in the
  //                        input file, the rest of the file until EOF
  //                        will be loaded.
  // Return:
  //    Returns 0 if successful, -1 if no more sequences can be loaded
  //    (i.e. EOF), and otherwise if unsuccessful.
  int LoadNextBatchNSequences(uint64_t num_seqs_to_load, bool randomize_non_acgt_bases=false);

  // This function loads only a part of the sequences present in the
  // file into this object. Before using LoadNextBatchInMegabytes a file
  // needs to explicitly be opened with the OpenFileForBatchLoading function.
  // This function loads only a fixed amount of MB of sequences at once.
  // When called the next time, it will load another batch fixed size,
  // until EOF is reached.
  // After all batches have been loaded from file, users must call
  // the CloseFileAfterBatchLoading function themselves.
  // Input:
  //    megabytes_to_load - Maximum amount of memory (in MB) to load in a single
  //                        batch. The loaded batch will exceed the value of
  //                        this parameter by the size of the last loaded
  //                        sequence. If the parameter's value is bigger
  //                        than the file size, the entire file will be loaded.
  // Return:
  //    Returns 0 if successful, -1 if no more sequences can be loaded
  //    (i.e. EOF), and otherwise if unsuccessful.
  int LoadNextBatchInMegabytes(uint64_t megabytes_to_load, bool randomize_non_acgt_bases=false);

  // Calculates the size of the sequences that it currently occupies in memory.
  // Calculation of size includes the sum of lengths of the headers, the data
  // and the quality scores of all sequences. It does not include other member
  // variables of this class.
  // Inputs:
  //    memory_unit - specifies the memory unit in which the return value
  //                  should be given (byte, kB, MB or GB). These are
  //                  predefined as constants: MEMORY_UNIT_BYTE,
  //                  MEMORY_UNIT_KILOBYTE, MEMORY_UNIT_MEGABYTE and
  //                  MEMORY_UNIT_GIGABYTE.
  // Return:
  //    Returns the size of the sequences occupying this container.
  uint64_t CalculateTotalSize(int32_t memory_unit=MEMORY_UNIT_BYTE);

  // Returns the sum of number of bases of all the sequences in the file.
  uint64_t GetNumberOfBases();

  // Outputs the contents of this object to the stream given by file pointer.
  // Inputs:
  //    fp  - file pointer to an open file. Can also be stdout and stderr.  void Verbose(FILE *fp);
  void Verbose(FILE *fp);

  const SequenceVector& get_sequences() const;
  void set_sequences(const SequenceVector& sequences);

  uint64_t get_current_batch_starting_sequence_id() const {
    return current_batch_starting_sequence_id_;
  }

  void set_current_batch_starting_sequence_id(uint64_t currentBatchStartingSequenceId) {
    current_batch_starting_sequence_id_ = currentBatchStartingSequenceId;
  }

 private:
  SequenceVector sequences_;  // Vector holding all the sequences in the file (or in a batch).
  std::string open_file_path_;  // Path to the sequences file that is currently opened (i.e. during batch loading).
  kseq_t *bwa_seq_;  // Variable used by BWA's functions for file parsing.
  gzFile bwa_fp_;  // File pointer for an opened FASTA/FASTQ file.
  uint64_t current_batch_id_;  // ID of the current batch (ordinal number).
  uint64_t current_batch_starting_sequence_id_;  // Absolute ID of the first sequence in this object. If all sequences loaded at once is equal to 0, otherwise to the absolute ID of the starting sequence in the batch.
  uint64_t current_data_size_;  // When new sequences are added to this object
                                // using the AddSequence, their size (in bytes)
                                // is automatically calculated. Used for batch
                                // loading of fixed size of sequences.

};

#endif /* SEQUENCEFILE_H_ */
