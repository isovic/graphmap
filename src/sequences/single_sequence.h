/*
 * Copyright 2014, Ivan Sovic.
 * All rights reserved.
 *
 * SingleSequence.h
 *
 *  Created on: 15 May, 2014
 *      Author: Ivan Sovic
 */

#ifndef SINGLESEQUENCE_H_
#define SINGLESEQUENCE_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include "utility/utility_conversion-inl.h"

enum DataFormat {
  kDataFormatAscii = 0,
  kDataFormat2BitSparse,
  kDataFormat2BitPacked,  /// Bases are stored in the 2-bit format (A = 00, C = 01, G = 10, T = 11). Byte is packed from MSB for the first base to LSB for the last base. I.e. the first base occupies the 2 MSB bits of a byte, and the last base two LSB bits of a byte.
  kDataFormat2BitPacked2  /// The same as kDataFormat2BitPacked, but bases in a byte are stored in a normal ascending order. I.e. the first base occupies 2 LSB bits of a byte, and the last base 2 MSB bits of a byte.
};

// A container class that holds a single genomic sequence, for example a
// single entry in a FASTA/FASTQ file. A SingleSequence is defined with
// three components, the same as in a FASTQ files: header, sequence and
// quality. This class provides options for conversion to/from three
// different formats for sequence data representation (only the sequence,
// not the header and quality values):
//    1. plain ASCII format (bases are chars 'A', 'C', 'T', 'G', 'N', etc.)
//    2. 2 bit format, but a sparse one, meaning that one base still takes
//       one byte, but instead of using ASCII values for coding bases,
//       bases are now encoded as A = 0, C = 1, G = 2 and T = 3. Inverse
//       of a base can be obtained by expression (3 - base). BWA uses this
//       representation for index searching.
//    3. 2 bit packed (dense) format, where 4 bases are packed in one byte.
//       Encoding of bases is the same as in the sparse 2 bit format:
//       A = 0, C = 1, G = 2 and T = 3. This is useful to reduce the size
//       of a sequence in memory (or on disk) when it is not continuously
//       needed (i.e. after indexing), but only occasionally.
//
// Sample usage:
//    int64_t sequence_id = 0;
//    std::string header = "This is a FASTA header!";
//    std::string data = "ACCTCTGCAAAAAC";
//    SingleSequence sequence;
//    sequence.InitHeaderAndDataFromAscii(header.c_str(), header.size(), (int8_t *) data.c_str(), data.size(), sequence_id);
//
class SingleSequence {
 public:
  SingleSequence();
  ~SingleSequence();

  // Sets all member variables to initial values, and frees any memory
  // if it was allocated, followed by setting all pointers to NULL value.
  void Clear();

  // Sets only the header of a SingleSequence object.
  // Inputs:
  //    header  - C-style formatted string.
  //    header_length - length of the header string (in case it is
  //                    not \0 terminated.
  //  Return:
  //    Returns 0 if successful.
  int InitHeader(char *header, uint32_t header_length);

  // Sets only the sequence data of a SingleSequence object. If this function
  // is used for data initialization, the SingleSequence object will be marked
  // that it is in the ASCII format.
  // Inputs:
  //    data  - array of chars (casted to (int8_t*)) with bases
  //            encoded in ASCII format.
  //    sequence_length - length of the sequence (number of bases).
  //  Return:
  //    Returns 0 if successful.
  int InitDatafromAscii(int8_t *data, uint64_t sequence_length);

  // Sets only the sequence data of a SingleSequence object. If this function
  // is used for data initialization, the SingleSequence object will be marked
  // that it is in the 2 bit sparse format.
  // Inputs:
  //    data  - array with bases encoded in 2 bit sparse
  //            format (1 base per byte).
  //    sequence_length - length of the sequence (number of bases).
  //  Return:
  //    Returns 0 if successful.
  int InitDataFrom2BitSparse(int8_t *data, uint64_t sequence_length);

  // Sets only the sequence data of a SingleSequence object. If this function
  // is used for data initialization, the SingleSequence object will be marked
  // that it is in the 2 bit packed format (4 bases per byte).
  // Inputs:
  //    data  - array with bases encoded in 2 bit packed (dense)
  //            format (4 bases per byte).
  //    data_length - length of the packed sequence (different than the
  //                  sequence length). Concretelly,
  //                  data_length = ceil(sequence_length / 4)
  //    sequence_length - length of the sequence (number of bases).
  //  Return:
  //    Returns 0 if successful.
  int InitDataFrom2BitPacked(int8_t *data, uint64_t data_length,
                             uint64_t sequence_length);

  // Sets only the quality values of a SingleSequence object.
  // Inputs:
  //    quality  - Array of quality values, for example Phred scores.
  //    quality_length - length of the quality array.
  //  Return:
  //    Returns 0 if successful.
  int InitQuality(int8_t *quality, uint64_t quality_length);

  // Initializes both header and data of a SingleSequence object, setting
  // the object to ASCII format. The function is implemented by calling the
  // InitHeader and the InitDatafromAscii functions.
  // Inputs:
  //    header  - C-style formatted string.
  //    header_length - length of the header string (in case it is
  //                    not \0 terminated.
  //    data  - array of chars (casted to (int8_t*)) with bases
  //            encoded in ASCII format.
  //    sequence_length - length of the sequence (number of bases).
  //    sequence_id - if this sequence object is part of an array of
  //                  sequences (i.e. a sequence file), sequence_id can
  //                  be used to mark its ordinal number in the array.
  //                  This parameter is not necessary, and has a default
  //                  value of -2.
  //  Return:
  //    Returns 0 if successful.
  int InitHeaderAndDataFromAscii(char *header, uint32_t header_length,
                                 int8_t *data, uint64_t data_length,
                                 int64_t sequence_id = -2, int64_t sequence_absolute_id = -2);

  // Initializes both header and data of a SingleSequence object, setting
  // the object to 2 bit sparse format. The function is implemented by calling the
  // InitHeader and the InitDataFrom2BitSparse functions.
  // Inputs:
  //    header  - C-style formatted string.
  //    header_length - length of the header string (in case it is
  //                    not \0 terminated.
  //    data  - array with bases encoded in 2 bit sparse
  //            format (1 base per byte).
  //    sequence_length - length of the sequence (number of bases).
  //    sequence_id - if this sequence object is part of an array of
  //                  sequences (i.e. a sequence file), sequence_id can
  //                  be used to mark its ordinal number in the array.
  //                  This parameter is not necessary, and has a default
  //                  value of -2.
  //  Return:
  //    Returns 0 if successful.
  int InitHeaderAndDataFrom2BitSparse(char *header, uint32_t header_length,
                                      int8_t *data, uint64_t sequence_length,
                                      int64_t sequence_id = -2, int64_t sequence_absolute_id = -2);

  // Initializes both header and data of a SingleSequence object, setting
  // the object to 2 bit packed format. The function is implemented by calling the
  // InitHeader and the InitDataFrom2BitPacked functions.
  // Inputs:
  //    header  - C-style formatted string.
  //    header_length - length of the header string (in case it is
  //                    not \0 terminated.
  //    data  - array with bases encoded in 2 bit packed (dense)
  //            format (4 bases per byte).
  //    data_length - length of the packed sequence (different than the
  //                  sequence length). Concretelly,
  //                  data_length = ceil(sequence_length / 4)
  //    sequence_length - length of the sequence (number of bases).
  //    sequence_id - if this sequence object is part of an array of
  //                  sequences (i.e. a sequence file), sequence_id can
  //                  be used to mark its ordinal number in the array.
  //                  This parameter is not necessary, and has a default
  //                  value of -2.
  //  Return:
  //    Returns 0 if successful.
  int InitHeaderAndDataFrom2BitPacked(char *header, uint32_t header_length,
                                      int8_t *data, uint64_t data_length,
                                      uint64_t sequence_length,
                                      int64_t sequence_id = -2, int64_t sequence_absolute_id = -2);

  // Initializes both header and data of a SingleSequence object, setting
  // the object to ASCII format. The function is implemented by calling the
  // InitHeader, InitDatafromAscii and InitQuality functions. Please refer to
  // descriptions of these functions for more details. Additionally,
  // sequence_id can be used as the ordinal number of this sequence, e.g. if
  // the sequence object is part of an array of sequences
  // (i.e. a sequence file), sequence_id can be used to mark its ordinal number
  // in the array. This parameter is not necessary, and has a default value of -2.
  int InitAllFromAscii(char *header, uint32_t header_length, int8_t *data,
                       int8_t *quality, uint64_t sequence_length,
                       int64_t sequence_id = -2, int64_t sequence_absolute_id = -2);

  // Initializes both header and data of a SingleSequence object, setting
  // the object to 2 bit sparse format. The function is implemented by calling
  // InitHeader, InitDatafrom2BitSparse and InitQuality functions. Please refer
  // to descriptions of these functions for more details. Additionally,
  // sequence_id can be used as the ordinal number of this sequence, e.g. if
  // the sequence object is part of an array of sequences
  // (i.e. a sequence file), sequence_id can be used to mark its ordinal number
  // in the array. This parameter is not necessary, and has a default value of -2.
  int InitAllFrom2BitSparse(char *header, uint32_t header_length, int8_t *data,
                            int8_t *quality, uint64_t sequence_length,
                            int64_t sequence_id = -2, int64_t sequence_absolute_id = -2);

  // Initializes both header and data of a SingleSequence object, setting
  // the object to 2 bit packed format. The function is implemented by calling
  // InitHeader, InitDatafrom2BitPacked and InitQuality functions. Please refer
  // to descriptions of these functions for more details. Additionally,
  // sequence_id can be used as the ordinal number of this sequence, e.g. if
  // the sequence object is part of an array of sequences
  // (i.e. a sequence file), sequence_id can be used to mark its ordinal number
  // in the array. This parameter is not necessary, and has a default value -2.
  int InitAllFrom2BitPacked(char *header, uint32_t header_length, int8_t *data,
                            int8_t *quality, uint64_t data_length,
                            uint64_t sequence_length, int64_t sequence_id = -2, int64_t sequence_absolute_id = -2);

  // Calculates the size of the sequence that it currently occupies in memory.
  // Calculation of size includes the length of the header, the data and the
  // quality scores. It does not include other member variables of this class.
  // Inputs:
  //    memory_unit - specifies the memory unit in which the return value
  //                  should be given (byte, kB, MB or GB). These are
  //                  predefined as constants: MEMORY_UNIT_BYTE,
  //                  MEMORY_UNIT_KILOBYTE, MEMORY_UNIT_MEGABYTE and
  //                  MEMORY_UNIT_GIGABYTE.
  // Return:
  //    Returns the size of the sequence occupying this container.
  uint64_t CalculateTotalSize(int32_t memory_unit = MEMORY_UNIT_BYTE);

  // Converts the sequence data between ASCII, 2 bit sparse and 2 bit dense
  // formats. Sequence data of this object is replaced with the converted
  // value.
  // Inputs:
  //    new_data_format - specifies the data format to which data should
  //                      be converted to. Can be equal to kDataFormatAscii,
  //                      kDataFormat2BitSparse or kDataFormat2BitPacked. If
  //                      conversion to the same format is attempted, warning
  //                      will be issued and logged.
  int ConvertDataFormat(DataFormat new_data_format);

  // Converts the sequence data between ASCII, 2 bit sparse and 2 bit dense
  // formats. The converted value is returned as a new object. User should
  // free the memory of the returned object themself to avoid memory leaks.
  // (use delete return_value).
  // Inputs:
  //    new_data_format - specifies the data format to which data should
  //                      be converted to. Can be equal to kDataFormatAscii,
  //                      kDataFormat2BitSparse or kDataFormat2BitPacked. If
  //                      conversion to the same format is attempted, warning
  //                      will be issued and logged.
  SingleSequence* ConvertDataFormatAndReturn(DataFormat new_data_format) const;

  // Returns the number of bases in the stored sequence.
  uint64_t GetNumberOfBases() const;

  // Returns the number of non ACTG (or actg) bases in the sequence.
  uint64_t CalcNumberNBases() const;

  // Checks if the initialized data contains any non-[ACTG] bases and changes
  // them with randomly (uniformly) chosen bases.
  int RandomizeNonACGTBases();

  // Reverses and complements the sequence data, and returns the new the new
  // sequence.
  int8_t* GetReverseComplement() const;

  // Reverses and complements the sequence data, and returns the new the new
  // sequence as a std::string.
  std::string GetReverseComplementAsString() const;

  // Reverses the quality scores, if they are loaded.
  int8_t* GetReverseQuality() const;

  // Reverses the quality scores, if they are loaded.
  std::string GetReverseQualityAsString() const;

  // Reverses and complements the sequence data, and replaces the data already
  // stored in this object.
  int ReverseComplement();

  // Copy function instead of operator= or a copy constructor. Done so
  // in order to stand by Google CppLint convention.
  void CopyFrom(const SingleSequence &op1);

  // Extracts a substring from the sequence, and converts it to ASCII
  // format if needed, and returns the substring as std::string.
  std::string GetSubstring(uint64_t start=0, uint64_t end=0) const;

  // Outputs the contents of this object to the stream given by file pointer.
  // Inputs:
  //    fp  - file pointer to an open file. Can also be stdout and stderr.
  void Verbose(FILE *fp) const;

  void set_header(char *header);
  char *get_header() const;
  void set_data(int8_t *data);
  const int8_t* get_data() const;
  void set_quality(int8_t *quality);
  const int8_t* get_quality() const;
  DataFormat get_data_format() const;
  uint64_t get_data_length() const;
  void set_data_length(uint64_t dataLength);
  uint32_t get_header_length() const;
  void set_header_length(uint32_t headerLength);
  uint64_t get_quality_length() const;
  void set_quality_length(uint64_t qualityLength);
  int64_t get_sequence_id() const;
  void set_sequence_id(int64_t sequenceId);
  uint64_t get_sequence_length() const;
  void set_sequence_length(uint64_t sequenceLength);

  int64_t get_sequence_absolute_id() const {
    return sequence_absolute_id_;
  }

  void set_sequence_absolute_id(int64_t sequenceAbsoluteId) {
    sequence_absolute_id_ = sequenceAbsoluteId;
  }

 private:
  // Allocates memory needed for storing the data, and copies values given by
  // parameters. If the data has been previously allocated, the memory is
  // freed and replaced with newly allocated data.
  // Inputs:
  //    data  - array holding the data, in any of the supported data formats.
  //    data_length - length of the data array. If the data is in 2bit packed
  //                  format, then data_length = ceil(sequence_length / 4),
  //                  otherwise data_length = sequence_length;
  //    sequence_length - length of the sequence (number of bases).
  int InitData_(int8_t *data, uint64_t data_length, uint64_t sequence_length);

  // Converts the sequence data from ASCII to any of supported data formats.
  // Sequence data of this object is replaced with the converted
  // value.
  // Inputs:
  //    new_data_format - specifies the data format to which data should
  //                      be converted to. Can be equal to kDataFormatAscii,
  //                      kDataFormat2BitSparse or kDataFormat2BitPacked. If
  //                      conversion to the same format is attempted, warning
  //                      will be issued and logged.
  int ConvertDataFormatFromAscii_(DataFormat new_data_format);

  // Converts the sequence data from all supported data formats to ASCII.
  // ASCII format is used as the reference format for conversion between other
  // formats. Sequence data of this object is replaced with the converted
  // value.
  int ConvertDataFormatToAscii_();

  char *header_;  // C-style string that holds the header of a sequence (i.e. FASTA or FASTQ headers).
  int8_t *data_;  // Sequence data which can be stored in any of the supported data formats.
  int8_t *quality_;  // Quality scores, e.g. if FASTQ files were loaded. Otherwise, equal to NULL.

  uint32_t header_length_;  // Length of the header, should be equal to strlen(header_).
  uint64_t data_length_;  // Length of the data container. If data is in ASCII format, then data_length_ equals sequence_length_, otherwise the number of bytes the data is stored in. I.e. in 2bit format, it would be ceil(sequence_length_ / 4).
  uint64_t sequence_length_;  // Length of the genomic sequence, given in the number of bases.
  uint64_t quality_length_;  // Length of the quality score, should be equal to sequence_length.

  int64_t sequence_id_;  // ID of this sequence object. E.g the ordinal number of the sequence in a batch.
  int64_t sequence_absolute_id_;  // ID of the sequence in the sequence file.

  DataFormat data_format_;  // Specifies whether data is in ASCII, 2bit sparse or 2bit packed format.
};

#endif /* SINGLESEQUENCE_H_ */
