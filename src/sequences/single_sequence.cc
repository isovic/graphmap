/*
 * Copyright 2014, Ivan Sovic.
 * All rights reserved.
 *
 * SingleSequence.cc
 *
 *  Created on: 15 May, 2014
 *      Author: Ivan Sovic
 */

#include <random>
//#include "bwa/kseq.h"
//#include "bwa/bwa.h"
//#include "bwa/bwamem.h"
//#include "bwa/kvec.h"
//#include "bwa/utils.h"
//#include "bwa/utils.h"
#include "sequences/single_sequence.h"
#include "log_system/log_system.h"
#include <string.h>

SingleSequence::SingleSequence() {
  header_ = NULL;
  data_ = NULL;
  quality_ = NULL;

  header_length_ = 0;
  data_length_ = 0;
  sequence_length_ = 0;
  quality_length_ = 0;

  sequence_id_ = -1;
  sequence_absolute_id_ = -1;

  data_format_ = kDataFormatAscii;
}

SingleSequence::~SingleSequence() {
  Clear();
}

void SingleSequence::Clear() {
  if (header_)
    delete[] header_;
  header_ = NULL;

  if (data_)
    delete[] data_;
  data_ = NULL;

  if (quality_)
    delete[] quality_;
  quality_ = NULL;

  header_length_ = 0;
  data_length_ = 0;
  sequence_length_ = 0;
  quality_length_ = 0;

  data_format_ = kDataFormatAscii;
}

void SingleSequence::set_header(char* header) {
  header_ = header;
}

char* SingleSequence::get_header() const {
  return header_;
}

void SingleSequence::set_data(int8_t* data) {
  data_ = data;
}

const int8_t* SingleSequence::get_data() const {
  return data_;
}

void SingleSequence::set_quality(int8_t* quality) {
  quality_ = quality;
}

const int8_t* SingleSequence::get_quality() const {
  return quality_;
}

int SingleSequence::InitHeader(char *header, uint32_t header_length) {
  if (header_)
    delete[] header_;
  header_ = NULL;


//  uint32_t trimmed_header_length = 0;
//  for (trimmed_header_length=0; trimmed_header_length<header_length; trimmed_header_length++)
//    if (header[trimmed_header_length] == ' ')
//      break;
//  header_length = trimmed_header_length;

  header_ = (char *) new char[header_length + 1];

  if (header_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: header_."));
    return 1;
  }

  // memcpy is not used because it was updated recently and requires new GLIBC symbols,
  // thus binaries won't run on older versions of systems.
  memmove(header_, header, header_length);

//  for (uint32_t i=0; i<header_length; i++)
//    if (header_[i] == ' ')
//      header_[i] = '_';

  header_length_ = header_length;
  header_[header_length] = '\0';

  return 0;
}

int SingleSequence::InitQuality(int8_t* quality, uint64_t quality_length) {
  if (quality_)
    delete[] quality_;
  quality_ = NULL;

  quality_ = (int8_t *) new int8_t[quality_length + 1];

  if (quality_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: quality_."));
    return 1;
  }

  // memcpy is not used because it was updated recently and requires new GLIBC symbols,
  // thus binaries won't run on older versions of systems.
  memmove(quality_, quality, quality_length);
  quality_[quality_length] = '\0';
  quality_length_ = quality_length;

  return 0;
}

int SingleSequence::InitData_(int8_t* data, uint64_t data_length,
                              uint64_t sequence_length) {
  if (data_)
    delete[] data_;
  data_ = NULL;

  data_ = (int8_t *) new int8_t[data_length + 1];

  if (data_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: data_."));
    return 1;
  }

  // memcpy is not used because it was updated recently and requires new GLIBC symbols,
  // thus binaries won't run on older versions of systems.
  memmove(data_, data, data_length);
  data_[data_length] = (int8_t) '\0';
  data_length_ = data_length;
  sequence_length_ = sequence_length;

  return 0;
}

int SingleSequence::InitDatafromAscii(int8_t* data, uint64_t sequence_length) {
  data_format_ = kDataFormatAscii;
  return InitData_(data, sequence_length, sequence_length);
}

int SingleSequence::InitDataFrom2BitPacked(int8_t* data, uint64_t data_length,
                                           uint64_t sequence_length) {
  data_format_ = kDataFormat2BitPacked;
  return InitData_(data, data_length, sequence_length);
}

int SingleSequence::InitDataFrom2BitSparse(int8_t* data,
                                           uint64_t sequence_length) {
  data_format_ = kDataFormat2BitSparse;
  return InitData_(data, sequence_length, sequence_length);
}

int SingleSequence::InitHeaderAndDataFromAscii(char *header,
                                               uint32_t header_length,
                                               int8_t* data,
                                               uint64_t data_length,
                                               int64_t sequence_id, int64_t sequence_absolute_id) {
  if (sequence_id != -2)
    sequence_id_ = sequence_id;
  if (sequence_absolute_id != -2)
    sequence_absolute_id_ = sequence_absolute_id;

  return (InitHeader(header, header_length)
      | InitDatafromAscii(data, data_length));
}

int SingleSequence::InitHeaderAndDataFrom2BitSparse(char *header,
                                                    uint32_t header_length,
                                                    int8_t* data,
                                                    uint64_t sequence_length,
                                                    int64_t sequence_id, int64_t sequence_absolute_id) {
  if (sequence_id != -2)
    sequence_id_ = sequence_id;
  if (sequence_absolute_id != -2)
    sequence_absolute_id_ = sequence_absolute_id;

  return (InitHeader(header, header_length)
      | InitDataFrom2BitSparse(data, sequence_length));
}

int SingleSequence::InitHeaderAndDataFrom2BitPacked(char *header,
                                                    uint32_t header_length,
                                                    int8_t* data,
                                                    uint64_t data_length,
                                                    uint64_t sequence_length,
                                                    int64_t sequence_id, int64_t sequence_absolute_id) {
  if (sequence_id != -2)
    sequence_id_ = sequence_id;
  if (sequence_absolute_id != -2)
    sequence_absolute_id_ = sequence_absolute_id;

  return (InitHeader(header, header_length)
      | InitDataFrom2BitPacked(data, data_length, sequence_length));
}

int SingleSequence::InitAllFromAscii(char *header, uint32_t header_length,
                                     int8_t *data, int8_t *quality,
                                     uint64_t sequence_length,
                                     int64_t sequence_id, int64_t sequence_absolute_id) {
  if (sequence_id != -2)
    sequence_id_ = sequence_id;
  if (sequence_absolute_id != -2)
    sequence_absolute_id_ = sequence_absolute_id;

  return (InitHeader(header, header_length)
      | InitDatafromAscii(data, sequence_length)
      | InitQuality(quality, sequence_length));
}

int SingleSequence::InitAllFrom2BitPacked(char *header, uint32_t header_length,
                                          int8_t *data, int8_t *quality,
                                          uint64_t data_length,
                                          uint64_t sequence_length,
                                          int64_t sequence_id, int64_t sequence_absolute_id) {
  if (sequence_id != -2)
    sequence_id_ = sequence_id;
  if (sequence_absolute_id != -2)
    sequence_absolute_id_ = sequence_absolute_id;

  return (InitHeader(header, header_length)
      | InitDataFrom2BitPacked(data, data_length, sequence_length)
      | InitQuality(quality, data_length));
}

int SingleSequence::InitAllFrom2BitSparse(char *header, uint32_t header_length,
                                          int8_t *data, int8_t *quality,
                                          uint64_t sequence_length,
                                          int64_t sequence_id, int64_t sequence_absolute_id) {
  if (sequence_id != -2)
    sequence_id_ = sequence_id;
  if (sequence_absolute_id != -2)
    sequence_absolute_id_ = sequence_absolute_id;

  return (InitHeader(header, header_length)
      | InitDataFrom2BitSparse(data, sequence_length)
      | InitQuality(quality, sequence_length));
}

uint64_t SingleSequence::CalculateTotalSize(int32_t memory_unit) {
  uint64_t total_sequence_size = 0;

  total_sequence_size = header_length_ * sizeof(*header_)
      + data_length_ * sizeof(*data_) + quality_length_ * sizeof(*quality_);

  if (memory_unit == MEMORY_UNIT_BYTE)
    total_sequence_size = total_sequence_size / ((uint64_t) 1);
  else if (memory_unit == MEMORY_UNIT_KILOBYTE)
    total_sequence_size = total_sequence_size / ((uint64_t) 1024);
  else if (memory_unit == MEMORY_UNIT_MEGABYTE)
    total_sequence_size = total_sequence_size
        / (((uint64_t) 1024) * ((uint64_t) 1024));
  else if (memory_unit == MEMORY_UNIT_GIGABYTE)
    total_sequence_size = total_sequence_size
        / (((uint64_t) 1024) * ((uint64_t) 1024) * ((uint64_t) 1024));
  else
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Memory unit not recognized! Returning value in bytes."));

  return total_sequence_size;
}

void SingleSequence::Verbose(FILE *fp) const {

  fprintf(fp, "Sequence ID: %ld\n", sequence_id_);

  if (data_format_ == kDataFormatAscii)
    fprintf(fp, "Data format: ASCII\n");
  else if (data_format_ == kDataFormat2BitSparse)
    fprintf(fp, "Data format: 2Bit Sparse\n");
  else if (data_format_ == kDataFormat2BitPacked)
    fprintf(fp, "Data format: 2Bit Packed\n");

  fprintf(fp, "Header: ");
  for (uint64_t current_char = 0; current_char < header_length_; current_char++)
    fprintf(fp, "%c", header_[current_char]);
  fprintf(fp, "\n");

  fprintf(fp, "Sequence: ");
  for (uint64_t current_base = 0; current_base < data_length_; current_base++) {
    if (data_format_ == kDataFormatAscii) {
      fprintf(fp, "%c", ((char) data_[current_base]));
    } else if (data_format_ == kDataFormat2BitSparse) {
      fprintf(fp, "%d", data_[current_base]);
    } else if (data_format_ == kDataFormat2BitPacked) {
      fprintf(fp, "%X_", (uint8_t) data_[current_base]);
    }

    if ((current_base & ((uint64_t) 0x03)) == 0x03)
      fprintf(fp, " ");
  }
  fprintf(fp, "\n");

  fprintf(fp, "Quality: ");
  for (uint64_t current_quality = 0; current_quality < quality_length_;
      current_quality++)
    fprintf(fp, "%c", (char) quality_[current_quality]);
  fprintf(fp, "\n");

  fflush(fp);
}

int SingleSequence::ConvertDataFormat(DataFormat new_data_format) {
  if (data_format_ == new_data_format) {
//    LogSystem::GetInstance().Log(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Targeted data format is the same as current. No changes will be made."));
    return 0;
  }

  ConvertDataFormatToAscii_();
  ConvertDataFormatFromAscii_(new_data_format);

  return 0;
}

// TODO: Implement ConvertAndReturnDataFormat function!!
SingleSequence* SingleSequence::ConvertDataFormatAndReturn(DataFormat new_data_format) const {
  if (data_format_ == new_data_format) {
    LogSystem::GetInstance().Log(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Targeted data format is the same as current. No changes will be made."));
    return NULL;
  }

  SingleSequence *converted_sequence = new SingleSequence();
  converted_sequence->CopyFrom((*this));

  converted_sequence->ConvertDataFormatToAscii_();
  converted_sequence->ConvertDataFormatFromAscii_(new_data_format);

  return converted_sequence;
}

int SingleSequence::ConvertDataFormatFromAscii_(DataFormat new_data_format) {

  if (new_data_format == kDataFormatAscii) {
    data_format_ = new_data_format;

  } else if (new_data_format == kDataFormat2BitSparse) {

    int8_t *new_data = new int8_t[sequence_length_ + 1];
    if (new_data == NULL) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: new_data."));
      return 1;
    }

    // Convert ASCII format to 2bit.
    for (uint64_t current_base = 0; current_base < data_length_;
        current_base++) {
      new_data[current_base] =
          (int8_t) kBaseToBwa[(int32_t) data_[current_base]];  // Convert from ASCII to 2bit encoding using a look-up table.
    }

    new_data[sequence_length_] = '\0';

    delete[] data_;
    data_ = new_data;
    data_length_ = sequence_length_;
    data_format_ = new_data_format;

  } else if (new_data_format == kDataFormat2BitPacked) {

    uint64_t new_data_length = (sequence_length_ >> 2)
        + (((sequence_length_ & 3) == 0) ? 0 : 1);
    int8_t *new_data = new int8_t[new_data_length + 1];
    if (new_data == NULL) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: new_data."));
      return 1;
    }

    // Convert ASCII format to 2bit packed.
    // Explanation:
    // (current_base & 3) Takes only the last two bits of a base's index.
    //                    This is done because we are packing 4 bases to a
    //                    byte, so this is effectively modulo by 4.
    // ((3 - (current_base & 3))  2bit packed bases are stored in such a way
    //                            that the first base of a sequence is stored
    //                            as the Most Significant Bit (MSB) in the byte.
    //                            By subtracting the base in-byte index from 3
    //                            we are basically reversing the order of bases
    //                            in the byte. This value will be used for
    //                            for shifting once it is multiplied by 2.
    // ((3 - (current_base & 3)) << 1)  Multiplication by two, because we need
    //                                  the number of bits by which to shift.
    // The final expression just converts the ASCII base to 2bit, shifts it to
    // the left for calculated number of bits, and performs OR on the current
    // byte of the new_data to change the appropriate bits.
    // new_data[current_base >> 2]  Means that we are dividing the current base
    //                              by 4, thus packing base data in 4 times
    //                              less space.
    memset(new_data, 0, new_data_length);
    for (uint64_t current_base = 0; current_base < sequence_length_;
        current_base++) {
      new_data[current_base >> 2] |=
          kBaseToBwa[(int32_t) data_[current_base]]
              << ((3 - (current_base & 3)) << 1);
    }

    new_data[new_data_length] = '\0';

    delete[] data_;
    data_ = new_data;
    data_length_ = new_data_length;
    data_format_ = new_data_format;
  } else if (new_data_format == kDataFormat2BitPacked2) {

    uint64_t new_data_length = (sequence_length_ >> 2)
        + (((sequence_length_ & 3) == 0) ? 0 : 1);
    int8_t *new_data = new int8_t[new_data_length + 1];
    if (new_data == NULL) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: new_data."));
      return 1;
    }

    memset(new_data, 0, new_data_length);
    for (uint64_t current_base = 0; current_base < sequence_length_;
        current_base++) {
      new_data[current_base >> 2] |=
          kBaseToBwa[(int32_t) data_[current_base]]
              << (((current_base & 3)) << 1);
    }

    new_data[new_data_length] = '\0';

    delete[] data_;
    data_ = new_data;
    data_length_ = new_data_length;
    data_format_ = new_data_format;
  }
  return 0;
}

int SingleSequence::ConvertDataFormatToAscii_() {
  if (data_format_ == kDataFormatAscii) {
    data_format_ = kDataFormatAscii;
    return 0;  // Do nothing, the data is already in the right format.

  } else if (data_format_ == kDataFormat2BitPacked) {

    int8_t *new_data = new int8_t[sequence_length_ + 1];
    if (new_data == NULL) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: new_data."));
      return 1;
    }

    // Extract all bases from the 2bit dense (packed) sequence format, and
    // convert bases from 2bit encoding to ASCII.
    // Explanation:
    // data_[current_base >> 2] because there are 4 bases per byte
    //                          (data_ is of int8_t type), so divide index by 4.
    // (3 - (current_base & 3)) this is used for shifting the bases to the
    //                          beginning of the byte. (current_base & 3) is
    //                          the ordinal number of a base within the byte,
    //                          and (3 - (current_base & 3)) is the number
    //                          of 2bits that need to be shifted to the right.
    //                          This is done because the first base of a
    //                          sequence is stored as the Most Significant Bit
    //                          (MSB) in the 2bit packed representation.
    // ((3 - (current_base & 3)) << 1)  multiply by two, that is, convert
    //                                  the 2bits to number of bits by
    //                                  by which the packed byte needs
    //                                  to be shifted.
    // And finally, the entire expression has the & 3 at the end to
    // mask the remaining bits and remove junk.
    for (uint64_t current_base = 0; current_base < sequence_length_;
        current_base++) {
      new_data[current_base] = (int8_t) kBwaToBase[(data_[current_base >> 2]
          >> ((3 - (current_base & 3)) << 1) & 3)];
    }

    new_data[sequence_length_] = '\0';

    delete[] data_;
    data_ = new_data;
    data_length_ = sequence_length_;
    data_format_ = kDataFormatAscii;

    return 0;

  } else if (data_format_ == kDataFormat2BitSparse) {

    int8_t *new_data = new int8_t[sequence_length_ + 1];
    if (new_data == NULL) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: new_data."));
      return 1;
    }

    // Convert 2bit format to ASCII.
    for (uint64_t current_base = 0; current_base < data_length_;
        current_base++) {
      new_data[current_base] =
          (int8_t) kBwaToBase[(int32_t) data_[current_base]];  // Use the look-up table to convert from 2bit encoding to ASCII base names.
    }
    new_data[sequence_length_] = '\0';

    delete[] data_;
    data_ = new_data;
    data_length_ = sequence_length_;
    data_format_ = kDataFormatAscii;

    return 0;
  }

  LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Unknown targeted data format."));

  return 1;
}

void SingleSequence::CopyFrom(const SingleSequence &op1) {
  Clear();

  this->InitHeader(op1.header_, op1.header_length_);
  this->InitData_(op1.data_, op1.data_length_, op1.sequence_length_);
  this->InitQuality(op1.quality_, op1.quality_length_);

  this->sequence_id_ = op1.sequence_id_;
  this->data_format_ = op1.data_format_;
}

int8_t * SingleSequence::GetReverseComplement() const {
  int8_t *new_data = new int8_t[data_length_ + 1];
  if (new_data == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: new_data."));
    return NULL;
  }

  if (data_format_ == kDataFormat2BitPacked) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_NOT_IMPLEMENTED, ""));
    return NULL;
  }
  else {
    for (uint64_t current_base = 0; current_base < data_length_; current_base++) {
      if (data_format_ == kDataFormatAscii) {
        new_data[current_base] = kBaseComplement[data_[data_length_ - current_base - 1]];
      }
      else if (data_format_ == kDataFormat2BitSparse) {
        new_data[current_base] = 3 - data_[data_length_ - current_base - 1];
      }
    }
  }

  new_data[data_length_] = '\0';

  return new_data;
}

std::string SingleSequence::GetReverseComplementAsString() const {
  std::string ret;
  int8_t *new_data = new int8_t[data_length_ + 1];
  if (new_data == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: new_data."));
    return NULL;
  }

  if (data_format_ == kDataFormat2BitPacked) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_NOT_IMPLEMENTED, ""));
    return NULL;
  }
  else {
    for (uint64_t current_base = 0; current_base < data_length_; current_base++) {
      if (data_format_ == kDataFormatAscii) {
        new_data[current_base] = kBaseComplement[data_[data_length_ - current_base - 1]];
      }
      else if (data_format_ == kDataFormat2BitSparse) {
        new_data[current_base] = 3 - data_[data_length_ - current_base - 1];
      }
    }
  }

  new_data[data_length_] = '\0';

  ret = std::string((char *) new_data);

  if (new_data)
    delete[] new_data;
  new_data = NULL;

  return ret;
}

int8_t * SingleSequence::GetReverseQuality() const {
  if (quality_ == NULL || quality_length_ == 0) {
    return NULL;
  }

  int8_t *new_quality = new int8_t[quality_length_ + 1];
  if (new_quality == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: new_quality."));
    return NULL;
  }

  for (uint64_t current_base = 0; current_base < quality_length_; current_base++) {
    new_quality[current_base] = quality_[quality_length_ - current_base - 1];
  }

  new_quality[quality_length_] = '\0';

  return new_quality;
}

std::string SingleSequence::GetReverseQualityAsString() const {
  if (quality_ == NULL || quality_length_ == 0) {
    return std::string("");
  }

  std::string ret;

  int8_t *new_quality = new int8_t[quality_length_ + 1];
  if (new_quality == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: new_quality."));
    return std::string("");
  }

  for (uint64_t current_base = 0; current_base < quality_length_; current_base++) {
    new_quality[current_base] = quality_[quality_length_ - current_base - 1];
  }

  new_quality[quality_length_] = '\0';

  ret = std::string((char *) new_quality);

  if (new_quality)
    delete[] new_quality;
  new_quality = NULL;

  return ret;
}

int SingleSequence::ReverseComplement() {
  int8_t *new_data = GetReverseComplement();
  int8_t *new_quality = GetReverseQuality();

  if (new_data == NULL)
    return 1;

  if (data_)
    delete[] data_;
  data_ = new_data;

  if (quality_)
    delete[] quality_;
  quality_ = new_quality;

  return 0;
}

DataFormat SingleSequence::get_data_format() const {
  return data_format_;
}

uint64_t SingleSequence::get_data_length() const {
  return data_length_;
}

void SingleSequence::set_data_length(uint64_t dataLength) {
  data_length_ = dataLength;
}

uint32_t SingleSequence::get_header_length() const {
  return header_length_;
}

void SingleSequence::set_header_length(uint32_t headerLength) {
  header_length_ = headerLength;
}

uint64_t SingleSequence::get_quality_length() const {
  return quality_length_;
}

void SingleSequence::set_quality_length(uint64_t qualityLength) {
  quality_length_ = qualityLength;
}

int64_t SingleSequence::get_sequence_id() const {
  return sequence_id_;
}

void SingleSequence::set_sequence_id(int64_t sequenceId) {
  sequence_id_ = sequenceId;
}

uint64_t SingleSequence::get_sequence_length() const {
  return sequence_length_;
}

void SingleSequence::set_sequence_length(uint64_t sequenceLength) {
  sequence_length_ = sequenceLength;
}

int SingleSequence::RandomizeNonACGTBases() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, 3);

  printf ("Randomizing N bases!\n");
  fflush(stdout);

	if (data_format_ == kDataFormat2BitPacked)
		return 1;

	for (uint64_t current_base = 0; current_base < data_length_; current_base++) {
		if (data_format_ == kDataFormatAscii) {
			if (kBaseToBwa[((int8_t) data_[current_base])] > 3) {
				data_[current_base] = kBwaToBase[dis(gen)];
			}
		}
		else if (data_format_ == kDataFormat2BitSparse) {
			if (((int8_t) data_[current_base]) > 3)
				data_[current_base] = dis(gen);
		}
	}

	return 0;
}

std::string SingleSequence::GetSubstring(uint64_t start, uint64_t end) const {
  std::string ret = "";

  if (start > sequence_length_)
    start = 0;

  if (end == 0 || end > sequence_length_)
    end = sequence_length_;

  if (data_format_ == kDataFormat2BitPacked) {
    for (uint64_t current_base = start; current_base < sequence_length_ && current_base < end; current_base++) {
      ret += ((char) kBwaToBase[(data_[current_base >> 2]
          >> ((3 - (current_base & 3)) << 1) & 3)]);
    }
  }
  else {
    for (uint64_t current_base = start; current_base < sequence_length_ && current_base < end; current_base++) {
      if (data_format_ == kDataFormatAscii) {
        ret += ((char) data_[current_base]);
      }
      else if (data_format_ == kDataFormat2BitSparse) {
        ret += ((char) kBwaToBase[data_[current_base]]);
      }
    }
  }

  return ret;
}

uint64_t SingleSequence::GetNumberOfBases() const {
  return sequence_length_;
}

uint64_t SingleSequence::CalcNumberNBases() const {
  uint64_t ret = 0;

  if (data_format_ == kDataFormat2BitPacked) {
    for (uint64_t current_base = 0; current_base < sequence_length_; current_base++) {
      ret += (1 - kIsBase[kBwaToBase[(data_[current_base >> 2]
          >> ((3 - (current_base & 3)) << 1) & 3)]]);
    }
  }
  else {
    for (uint64_t current_base = 0; current_base < sequence_length_; current_base++) {
      if (data_format_ == kDataFormatAscii) {
        ret += (1 - kIsBase[data_[current_base]]);
      }
      else if (data_format_ == kDataFormat2BitSparse) {
        ret += (1 - kIsBase[kBwaToBase[data_[current_base]]]);
      }
    }
  }

  return ret;
}
