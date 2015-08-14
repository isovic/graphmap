/*
 * index.cpp
 *
 *  Created on: 23 May, 2014
 *      Author: Ivan Sovic
 */

#include <index/index.h>
#include "log_system/log_system.h"

Index::Index() {
  reference_starting_pos_.clear();
  reference_lengths_.clear();
  headers_.clear();
  data_length_ = 0;
  data_length_forward_ = 0;
  num_sequences_forward_ = 0;
  data_ptr_ = 0;
  num_sequences_ = 0;
  data_ = NULL;
}

Index::~Index() {
  if (data_)
    delete[] data_;
  data_ = NULL;
}

int Index::LoadFromFile(std::string index_path) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("Loading index from file.\n"), "LoadFromFile");
    Clear();
    FILE *fp = fopen(index_path.c_str(), "r");
    if (fp == NULL) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "Path: '%s'", index_path.c_str()));
      return 1;
    }
    int ret_deserialize = Deserialize_(fp);
    fclose(fp);
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("Index loaded.\n"), "LoadFromFile");
    return ret_deserialize;
}

int Index::GenerateFromFile(std::string sequence_file_path) {
  Clear();
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Loading reference from file to generate index.\n"), "GenerateFromFile");
  SequenceFile sequences(sequence_file_path);
  return GenerateFromSequenceFile(sequences);
}

int Index::GenerateFromSequenceFile(const SequenceFile& sequence_file) {
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Generating index from SequenceFile.\n"), "GenerateFromSequenceFile");

  Clear();

  clock_t time_start = clock();

  uint64_t total_data_length = 0;
  for (SequenceVector::const_iterator sequence_iterator = sequence_file.get_sequences().begin();
      sequence_iterator != sequence_file.get_sequences().end(); sequence_iterator++) {
    total_data_length += (*sequence_iterator)->get_sequence_length();
  }

  uint64_t mem_to_alloc = (total_data_length + sequence_file.get_sequences().size())*2;  // Special sign '!' will be added after every base, and there will be twice as many sequences because of reverse complements.

  data_ = new int8_t[mem_to_alloc];

  if (data_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: data_."));
    return 1;
  }

  data_length_ = mem_to_alloc;
  data_ptr_ = 0;

  InsertHeaders_(sequence_file);
  InsertSequencesIntoData_(sequence_file);
  data_length_forward_ = data_ptr_;
  num_sequences_forward_ = sequence_file.get_sequences().size();
  InsertReverseSequencesIntoData_(sequence_file);

  CreateIndex_(data_, data_length_);

  return 0;
}

int Index::GenerateFromSingleSequence(const SingleSequence &sequence) {
  Clear();

  clock_t time_start = clock();

  uint64_t total_data_length = 0;
  total_data_length += sequence.get_sequence_length();

  uint64_t mem_to_alloc = (total_data_length + 1)*2;  // Special sign '!' will be added after every base, and there will be twice as many sequences because of reverse complements.

  data_ = new int8_t[mem_to_alloc];

  if (data_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: data_."));
    return 1;
  }

  data_length_ = mem_to_alloc;
  data_ptr_ = 0;

  InsertSingleSequenceIntoData_(&sequence);
  data_length_forward_ = data_ptr_;
  num_sequences_forward_ = 1;
  InsertReverseSingleSequenceIntoData_(&sequence);

  CreateIndex_(data_, data_length_);

  headers_.clear();
  InsertHeader_(sequence.get_header(), sequence.get_header_length());

  return 0;
}

int Index::GenerateFromSingleSequenceOnlyForward(const SingleSequence &sequence) {
  Clear();

  uint64_t total_data_length = 0;
  total_data_length += sequence.get_sequence_length();

  uint64_t mem_to_alloc = (total_data_length + 1);  // Special sign '!' will be added after every base, and there will be twice as many sequences because of reverse complements.

  data_ = new int8_t[mem_to_alloc];

  if (data_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: data_."));
    return 1;
  }

  data_length_ = mem_to_alloc;
  data_ptr_ = 0;

  InsertSingleSequenceIntoData_(&sequence);
  data_length_forward_ = data_ptr_;
  num_sequences_forward_ = 1;

  CreateIndex_(data_, data_length_);

  headers_.clear();
  InsertHeader_(sequence.get_header(), sequence.get_header_length());

  return 0;
}

int Index::LoadOrGenerate(std::string reference_path, std::string out_index_path, bool verbose) {
  FILE *fp=NULL;
  fp = fopen(out_index_path.c_str(), "r");
  if (fp != NULL) {
    fclose(fp);
    int ret_load_from_file = LoadFromFile(out_index_path);

    if (ret_load_from_file == 0) {
//      Index needs to be rebuilt. It was generated using an older version of the format.
//      Index needs to be rebuilt. It was generated using a different indexing method.
      return 0;
    }

    if (verbose == true) {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, true, FormatString("Index needs to be rebuilt. It was generated using an older version.\n"), "LoadOrGenerate");
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("ret_load_from_file = %d\n", ret_load_from_file), "LoadOrGenerate");
      fflush(stderr);
    }
  }

  if (verbose == true) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, true, FormatString("Started generating new index from file '%s'...\n", reference_path.c_str()), "LoadOrGenerate");
    fflush(stderr);
  }

  GenerateFromFile(reference_path);

  if (verbose == true) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, true, FormatString("Storing new index to file '%s'...\n", out_index_path.c_str()), "LoadOrGenerate");
  }

  StoreToFile(out_index_path);

  if (verbose == true) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, true, FormatString("New index stored.\n"), "LoadOrGenerate");
  }

  return 0;
}

int Index::StoreToFile(std::string output_index_path) {
  clock_t time_start = clock();

  FILE *fp = fopen(output_index_path.c_str(), "w");

  if (fp == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "Path: '%s'", output_index_path.c_str()));
    return 1;
  }

  Serialize_(fp);

  fclose(fp);

  return 0;
}

int Index::InsertReverseSingleSequenceIntoData_(const SingleSequence *sequence) {
  SingleSequence *ascii_sequence = NULL;

  if (sequence->get_data_format() != kDataFormatAscii) {
    ascii_sequence = sequence->ConvertDataFormatAndReturn(kDataFormatAscii);
  }
  else {
    ascii_sequence = new SingleSequence();
    if (ascii_sequence == NULL) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: ascii_sequence."));
      return 1;
    }

    ascii_sequence->CopyFrom(*((SingleSequence *) (sequence)));
  }

  ascii_sequence->ReverseComplement();

  if (InsertSequence_((int8_t *) ascii_sequence->get_data(), ascii_sequence->get_sequence_length())) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Problem occurred during insertion of a new sequence to the index."));
    return 1;
  }

  if (ascii_sequence)
    delete ascii_sequence;
  ascii_sequence = NULL;

  return 1;
}

int Index::InsertHeaders_(const SequenceFile& sequence_file) {
  headers_.clear();
  for (SequenceVector::const_iterator sequence_iterator = sequence_file.get_sequences().begin();
      sequence_iterator != sequence_file.get_sequences().end(); sequence_iterator++) {
    InsertSingleHeader_((*sequence_iterator));
  }

  return 0;
}

int Index::InsertSequencesIntoData_(const SequenceFile& sequence_file) {
  if (data_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Data not initialized."));
    return 1;
  }

  for (SequenceVector::const_iterator sequence_iterator = sequence_file.get_sequences().begin(); sequence_iterator != sequence_file.get_sequences().end(); sequence_iterator++) {
    InsertSingleSequenceIntoData_((*sequence_iterator));
  }

  return 0;
}

int Index::InsertReverseSequencesIntoData_(const SequenceFile& sequence_file) {
  if (data_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Data not initialized."));
    return 1;
  }

  for (SequenceVector::const_iterator sequence_iterator = sequence_file.get_sequences().begin(); sequence_iterator != sequence_file.get_sequences().end(); sequence_iterator++) {
    InsertReverseSingleSequenceIntoData_((*sequence_iterator));
  }

  return 0;
}

int Index::InsertSingleHeader_(const SingleSequence *sequence) {
  return InsertHeader_(sequence->get_header(), sequence->get_header_length());
}

int Index::InsertHeader_(const char *header, uint64_t header_length) {
  std::string header_string(header);
  headers_.push_back(header_string);
  return 0;
}

int Index::Serialize_(FILE* fp_out) {
  uint64_t vector_length = 0;

  int64_t version_number = INDEX_VERSION;
//  IndexChunk::WriteChunk(fp_out, "VERSION", sizeof(version_number), (int8_t *) &version_number);
  char out_designator[] = "VERSION ";
  int64_t out_dlen = 1 * sizeof(int64_t);

//  std::string out_designator = "VERSION ";
//  int64_t out_dlen = 1;
//  if (out_designator.size() > 8) {
//    out_designator = out_designator.substr(0, 8);
//  } else if (out_designator.size() < 8) {
//    int num_empty_chars = 8 - out_designator.size();
//    out_designator += std::string(num_empty_chars, ' ');
//  }
  fwrite(out_designator, sizeof(char), 8, fp_out);
  fwrite(&out_dlen, sizeof(int64_t), 1, fp_out);
  fwrite(&version_number, sizeof(int64_t), 1, fp_out);



  fwrite(&num_sequences_, sizeof(num_sequences_), 1, fp_out);
  fwrite(&num_sequences_forward_, sizeof(num_sequences_forward_), 1, fp_out);

  for (uint64_t i=0; i<headers_.size(); i++) {
    uint64_t string_length = headers_[i].size();
    fwrite(&string_length, sizeof(string_length), 1, fp_out);
    fwrite(headers_[i].c_str(), sizeof(char), string_length, fp_out);
  }

  vector_length = reference_starting_pos_.size();
  fwrite(&vector_length, sizeof(vector_length), 1, fp_out);
  fwrite((reference_starting_pos_.data()), sizeof(uint64_t), vector_length, fp_out);

  vector_length = reference_lengths_.size();
  fwrite(&vector_length, sizeof(vector_length), 1, fp_out);
  fwrite((reference_lengths_.data()), sizeof(uint64_t), vector_length, fp_out);

  int ret_serialize_index = SerializeIndex_(fp_out);
  if (ret_serialize_index)
    return ret_serialize_index;

  fwrite(&data_length_, sizeof(data_length_), 1, fp_out);
  fwrite(&data_length_forward_, sizeof(data_length_forward_), 1, fp_out);
  fwrite(data_, sizeof(int8_t), data_length_, fp_out);

  return 0;
}

//void Index::Serialize_(FILE* fp_out) {
//  uint64_t vector_length = 0;
//
//  IndexChunk::WriteChunk(fp_out, "NUMSEQ", sizeof(num_sequences_), (int8_t *) &num_sequences_);
//  IndexChunk::WriteChunk(fp_out, "NUMSEQFW", sizeof(num_sequences_forward_), (int8_t *) &num_sequences_forward_);
//
////  fwrite(&num_sequences_, sizeof(num_sequences_), 1, fp_out);
////  fwrite(&num_sequences_forward_, sizeof(num_sequences_forward_), 1, fp_out);
//
//  for (uint64_t i=0; i<headers_.size(); i++) {
//    uint64_t string_length = headers_[i].size();
//    IndexChunk::WriteChunk(fp_out, "HEADER", sizeof(string_length), (int8_t *) headers_[i].c_str());
//
////    fwrite(&string_length, sizeof(string_length), 1, fp_out);
////    fwrite(headers_[i].c_str(), sizeof(char), string_length, fp_out);
//  }
//
//  IndexChunk::WriteChunk(fp_out, "REFST", sizeof(reference_starting_pos_.size()), (int8_t *) &reference_starting_pos_.size());
//  IndexChunk::WriteChunk(fp_out, "REFSTPOS", sizeof(uint64_t) * (int8_t *) &reference_starting_pos_.size(), (int8_t *) &reference_starting_pos_.data());
////  IndexChunk::WriteChunk(fp_out, "REFST", sizeof(reference_starting_pos_.size()), (int8_t *) &reference_starting_pos_.size());
////  IndexChunk::WriteChunk(fp_out, "REFSTPOS", sizeof(uint64_t) * (int8_t *) &reference_starting_pos_.size(), (int8_t *) &reference_starting_pos_.data());
//
////  vector_length = reference_starting_pos_.size();
////  fwrite(&vector_length, sizeof(vector_length), 1, fp_out);
////  fwrite((reference_starting_pos_.data()), sizeof(uint64_t), vector_length, fp_out);
//
//  vector_length = reference_lengths_.size();
//  fwrite(&vector_length, sizeof(vector_length), 1, fp_out);
//  fwrite((reference_lengths_.data()), sizeof(uint64_t), vector_length, fp_out);
//
//  SerializeIndex_(fp_out);
//
//  fwrite(&data_length_, sizeof(data_length_), 1, fp_out);
//  fwrite(&data_length_forward_, sizeof(data_length_forward_), 1, fp_out);
//  fwrite(data_, sizeof(int8_t), data_length_, fp_out);
//}

int Index::Deserialize_(FILE* fp_in) {
  uint64_t vector_length = 0;

//  int ret_value = 0;
//  ret_value = fseek(fp_in, 0, SEEK_END);
//  if (ret_value) {
//    return 1;
//  }
//
//  int64_t file_length = ftell(fp_in);
//  if (file_length == -1) {
//    return 1;
//  }
//
//  ret_value = fseek(fp_in, 0, SEEK_SET);
//  if (ret_value) {
//    return 1;
//  }
//
//  char *whole_file = (char *) malloc (sizeof(char) * file_length);



  char file_header[9];
//  IndexChunk::WriteChunk(fp_out, "VERSION", sizeof(version_number), (int8_t *) &version_number);
  if (fread(file_header, sizeof(char), 8, fp_in) != 8) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable file_header."));
    return 1;
  }
  file_header[8] = '\0';

//  printf ("file_header = '%s'\n", file_header);
//  fflush(stdout);

  if (std::string(file_header) != std::string("VERSION "))
    return -1;

//  Provjeriti da li je tekst stvarno jednak "VERSION", i ako nije, onda treba pokrenuti izgradnju indeksa ponovo!!

  int64_t temp_int = 0;
  if (fread(&temp_int, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable temp_int."));
    return 2;
  }
  int64_t version_number = 0;
  if (fread(&version_number, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable version_number."));
    return 3;
  }

//  printf ("temp_int = %ld\n", temp_int);
//  printf ("version_number = %ld\n", version_number);
//  std::cout << temp_int << "\n";
//  std::cout << version_number << "\n";
//  fflush(stdout);

  if (version_number != INDEX_VERSION) {
//    printf ("version_number = %ld\n", version_number);
//    fflush(stdout);
    return -3;
  }



  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Deserializing...\n"), "Deserialize_");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- num_sequences and num_sequences_forward_..."), "Deserialize_");

  if (fread(&num_sequences_, sizeof(num_sequences_), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable num_sequences_."));
    return 4;
  }

  if (fread(&num_sequences_forward_, sizeof(num_sequences_forward_), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable num_sequences_forward_."));
    return 5;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- headers_..."), "Deserialize_");
  headers_.clear();
  for (uint64_t i=0; i<num_sequences_forward_; i++) {
    uint64_t string_length = 0;
    if (fread(&string_length, sizeof(string_length), 1, fp_in) != 1) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable string_length."));
      return 6;
    }

    char *new_header = new char[string_length + 1];
    if (fread(new_header, sizeof(char), string_length, fp_in) != string_length) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable new_header."));
      return 7;
    }
    new_header[string_length] = '\0';
    headers_.push_back(std::string(new_header));
    delete[] new_header;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- reference_starting_pos_..."), "Deserialize_");
  if (fread(&vector_length, sizeof(vector_length), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_starting_pos_.size()."));
    return 8;
  }
  reference_starting_pos_.resize(vector_length);
  if (fread(reference_starting_pos_.data(), sizeof(uint64_t), vector_length, fp_in) != vector_length) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_starting_pos_."));
    return 9;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- reference_lengths_..."), "Deserialize_");
  if (fread(&vector_length, sizeof(vector_length), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_lengths_."));
    return 10;
  }
  reference_lengths_.resize(vector_length);
  if (fread(reference_lengths_.data(), sizeof(uint64_t), vector_length, fp_in) != vector_length) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_lengths_."));
    return 11;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  int ret_deserialize_index = DeserializeIndex_(fp_in);
  if (ret_deserialize_index)
    return (20 + ret_deserialize_index);

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- data_length_ and data_length_forward_..."), "Deserialize_");
  if (fread(&data_length_, sizeof(data_length_), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable data_length_."));
    return 12;
  }
  if (fread(&data_length_forward_, sizeof(data_length_forward_), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable data_length_forward_."));
    return 13;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("data_length_ = %ld, data_length_forward_ = %ld...", data_length_, data_length_forward_), "[]");
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- data_...\n"), "Deserialize_");
  if (data_)
    delete[] data_;
  data_ = new int8_t[data_length_ + 1];
  if (fread(data_, sizeof(int8_t), data_length_, fp_in) != data_length_) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable data_."));
    return 14;
  }
  data_[data_length_] = ((int8_t) '\0');
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "Deserialize_");

  return 0;
}

int Index::DeprecatedDeserialize_(FILE* fp_in) {
  uint64_t vector_length = 0;

  char file_header[9];
//  IndexChunk::WriteChunk(fp_out, "VERSION", sizeof(version_number), (int8_t *) &version_number);
  if (fread(file_header, sizeof(char), 8, fp_in) != 8) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable file_header."));
    return 1;
  }
  file_header[8] = '\0';

//  printf ("file_header = '%s'\n", file_header);
//  fflush(stdout);

  if (std::string(file_header) != std::string("VERSION "))
    return -1;

//  Provjeriti da li je tekst stvarno jednak "VERSION", i ako nije, onda treba pokrenuti izgradnju indeksa ponovo!!

  int64_t temp_int = 0;
  if (fread(&temp_int, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable temp_int."));
    return 2;
  }
  int64_t version_number = 0;
  if (fread(&version_number, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable version_number."));
    return 3;
  }

//  printf ("temp_int = %ld\n", temp_int);
//  printf ("version_number = %ld\n", version_number);
//  std::cout << temp_int << "\n";
//  std::cout << version_number << "\n";
//  fflush(stdout);

  if (version_number != INDEX_VERSION) {
//    printf ("version_number = %ld\n", version_number);
//    fflush(stdout);
    return -3;
  }



  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Deserializing...\n"), "Deserialize_");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- num_sequences and num_sequences_forward_..."), "Deserialize_");

  if (fread(&num_sequences_, sizeof(num_sequences_), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable num_sequences_."));
    return 4;
  }

  if (fread(&num_sequences_forward_, sizeof(num_sequences_forward_), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable num_sequences_forward_."));
    return 5;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- headers_..."), "Deserialize_");
  headers_.clear();
  for (uint64_t i=0; i<num_sequences_forward_; i++) {
    uint64_t string_length = 0;
    if (fread(&string_length, sizeof(string_length), 1, fp_in) != 1) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable string_length."));
      return 6;
    }

    char *new_header = new char[string_length + 1];
    if (fread(new_header, sizeof(char), string_length, fp_in) != string_length) {
      LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable new_header."));
      return 7;
    }
    new_header[string_length] = '\0';
    headers_.push_back(std::string(new_header));
    delete[] new_header;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- reference_starting_pos_..."), "Deserialize_");
  if (fread(&vector_length, sizeof(vector_length), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_starting_pos_.size()."));
    return 8;
  }
  reference_starting_pos_.resize(vector_length);
  if (fread(reference_starting_pos_.data(), sizeof(uint64_t), vector_length, fp_in) != vector_length) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_starting_pos_."));
    return 9;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- reference_lengths_..."), "Deserialize_");
  if (fread(&vector_length, sizeof(vector_length), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_lengths_."));
    return 10;
  }
  reference_lengths_.resize(vector_length);
  if (fread(reference_lengths_.data(), sizeof(uint64_t), vector_length, fp_in) != vector_length) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_lengths_."));
    return 11;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  int ret_deserialize_index = DeserializeIndex_(fp_in);
  if (ret_deserialize_index)
    return (20 + ret_deserialize_index);

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- data_length_ and data_length_forward_..."), "Deserialize_");
  if (fread(&data_length_, sizeof(data_length_), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable data_length_."));
    return 12;
  }
  if (fread(&data_length_forward_, sizeof(data_length_forward_), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable data_length_forward_."));
    return 13;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("data_length_ = %ld, data_length_forward_ = %ld...", data_length_, data_length_forward_), "[]");
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- data_...\n"), "Deserialize_");
  if (data_)
    delete[] data_;
  data_ = new int8_t[data_length_ + 1];
  if (fread(data_, sizeof(int8_t), data_length_, fp_in) != data_length_) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable data_."));
    return 14;
  }
  data_[data_length_] = ((int8_t) '\0');
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "Deserialize_");

  return 0;
}

/*
Index format specification v1.0:
Data is organized in chunks. Each chunk begins with a designator, followed by the data size in bytes, and then the data itself.
Chunks can be ordered in any way, except for the 'VERSION' chunk which needs to come first, always.

Chunk specification:
  - Designator (designator_) - 8 byte char array.
  - Data length (dlen_) - int64_t (8 byte signed) value indicating the data length in bytes.
  - Data (data_) - array of length dlen_ bytes.

Types of chunks currently supported:
  - "VERSION " - Contains the version of the index format used in the file. data_ of the chunk is a single int64_t number specifying the version number.
  - "TYPE    " - Specifies the type of the index that the file stores. data_ is a char array specifying the name of the index type in a C-style string manner. Valid types: "sa", "hash", "spacedhash".
  - "INFO    " - Additional info that can be stored in the file. data_ of the chunk contains a char array (a C-style string).
  - "NUMSEQ  " - Number of all sequences (forward and reverse strands) that the index encodes. data_ is an int64_t value.
  - "NUMSEQFW" - Number of forward sequences that the index encodes. The number specifies only the number of forward sequences. data_ is an int64_t value.
  - ""

Index needs to be rebuilt. It was generated using an older version of the format.
Index needs to be rebuilt. It was generated using a different indexing method.
 */

//void Index::Deserialize_(FILE* fp_in) {
//  uint64_t vector_length = 0;
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Deserializing...\n"), "Deserialize_");
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- num_sequences and num_sequences_forward_..."), "Deserialize_");
//
////  char version[9];
////  if (fread(version, sizeof(char), 8, fp_in) != 8) {
////    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable version."));
////    return;
////  }
////  if (fread(&num_sequences_, sizeof(num_sequences_), 1, fp_in) != 1) {
////    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable num_sequences_."));
////    return;
////  }
//
//
//
//  if (fread(&num_sequences_, sizeof(num_sequences_), 1, fp_in) != 1) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable num_sequences_."));
//    return;
//  }
//
//  if (fread(&num_sequences_forward_, sizeof(num_sequences_forward_), 1, fp_in) != 1) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable num_sequences_forward_."));
//    return;
//  }
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- headers_..."), "Deserialize_");
//  headers_.clear();
//  for (uint64_t i=0; i<num_sequences_forward_; i++) {
//    uint64_t string_length = 0;
//    if (fread(&string_length, sizeof(string_length), 1, fp_in) != 1) {
//      ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable string_length."));
//      return;
//    }
//
//    char *new_header = new char[string_length + 1];
//    if (fread(new_header, sizeof(char), string_length, fp_in) != string_length) {
//      ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable new_header."));
//      return;
//    }
//    new_header[string_length] = '\0';
//    headers_.push_back(std::string(new_header));
//    delete[] new_header;
//  }
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- reference_starting_pos_..."), "Deserialize_");
//  if (fread(&vector_length, sizeof(vector_length), 1, fp_in) != 1) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_starting_pos_.size()."));
//    return;
//  }
//  reference_starting_pos_.resize(vector_length);
//  if (fread(reference_starting_pos_.data(), sizeof(uint64_t), vector_length, fp_in) != vector_length) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_starting_pos_."));
//    return;
//  }
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- reference_lengths_..."), "Deserialize_");
//  if (fread(&vector_length, sizeof(vector_length), 1, fp_in) != 1) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_lengths_."));
//    return;
//  }
//  reference_lengths_.resize(vector_length);
//  if (fread(reference_lengths_.data(), sizeof(uint64_t), vector_length, fp_in) != vector_length) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable reference_lengths_."));
//    return;
//  }
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");
//
//  DeserializeIndex_(fp_in);
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- data_length_ and data_length_forward_..."), "Deserialize_");
//  if (fread(&data_length_, sizeof(data_length_), 1, fp_in) != 1) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable data_length_."));
//    return;
//  }
//  if (fread(&data_length_forward_, sizeof(data_length_forward_), 1, fp_in) != 1) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable data_length_forward_."));
//    return;
//  }
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("data_length_ = %ld, data_length_forward_ = %ld...", data_length_, data_length_forward_), "[]");
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "[]");
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- data_..."), "Deserialize_");
//  if (data_)
//    delete[] data_;
//  data_ = new int8_t[data_length_ + 1];
//  if (fread(data_, sizeof(int8_t), data_length_, fp_in) != data_length_) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable data_."));
//    return;
//  }
//  data_[data_length_] = ((int8_t) '\0');
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("done.\n"), "Deserialize_");
//}

int Index::InsertSequence_(const int8_t *sequence_data, uint64_t sequence_length) {
  reference_starting_pos_.push_back(data_ptr_);
  reference_lengths_.push_back(sequence_length);

  if (data_length_ < (data_ptr_ + sequence_length + 1)) {
    // Not enough preallocated memory.
    int8_t *data = new int8_t[data_ptr_ + sequence_length + 1];
    // memcpy is not used because it was updated recently and requires new GLIBC symbols,
    // thus binaries won't run on older versions of systems.
    memmove(data, data_, data_length_);
    if (data_)
      delete[] data_;

    data_ = data;
    data_length_ = data_ptr_ + sequence_length + 1;
  }

  num_sequences_ += 1;
  // Insert sequence in data + separator '!'
  for (uint64_t i = 0; i < sequence_length; ++i) {
    if ((sequence_data[i] == 'A' || sequence_data[i] == 'C' || sequence_data[i] == 'G' || sequence_data[i] == 'T')) {
      data_[data_ptr_++] = sequence_data[i];
    }
    else if ((sequence_data[i] == 'a' || sequence_data[i] == 'c' || sequence_data[i] == 'g' || sequence_data[i] == 't')) {
      data_[data_ptr_++] = toupper(sequence_data[i]);
    }
    else {
      data_[data_ptr_++] = 'N';
    }
  }

  data_[data_ptr_++] = '!';

  return 0;
}

int Index::InsertSingleSequenceIntoData_(const SingleSequence *sequence) {
  SingleSequence *ascii_sequence = NULL;
  if (sequence->get_data_format() != kDataFormatAscii) {
    ascii_sequence = sequence->ConvertDataFormatAndReturn(kDataFormatAscii);
  }
  else {
    ascii_sequence = ((SingleSequence *) (sequence));
  }

  if (InsertSequence_((int8_t *) ascii_sequence->get_data(), ascii_sequence->get_sequence_length())) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Problem occurred during insertion of a new sequence to the index."));
    return 1;
  }

  if (sequence->get_data_format() != kDataFormatAscii) {
    if (ascii_sequence)
      delete ascii_sequence;
    ascii_sequence = NULL;
  }

  return 0;
}


uint64_t Index::get_data_length_forward() const {
  return data_length_forward_;
}

const std::vector<std::string>& Index::get_headers() const {
  return headers_;
}

const int8_t* Index::get_data() const {
  return data_;
}

uint64_t Index::get_data_length() const {
  return data_length_;
}

uint64_t Index::get_num_sequences_forward() const {
  return num_sequences_forward_;
}

uint64_t Index::get_num_sequences() const {
  return num_sequences_forward_;
}

const std::vector<uint64_t>& Index::get_reference_lengths() const {
  return reference_lengths_;
}

const std::vector<uint64_t>& Index::get_reference_starting_pos() const {
  return reference_starting_pos_;
}

int64_t Index::RawPositionConverter(int64_t raw_position, int64_t query_length, int64_t *ret_absolute_position, int64_t *ret_relative_position, SeqOrientation *ret_orientation, int64_t *ret_reference_index_with_reverse) const {
  if (raw_position < 0 || raw_position >= data_length_)
    return -2;

  int64_t reference_index = RawPositionToReferenceIndexWithReverse(raw_position);
  if (reference_index < 0) {
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index. Values: reference_index = %ld, raw_position = %ld, data_length = %ld.", reference_index, raw_position, data_length_));
    return reference_index;
  }

  int64_t relative_pos = (raw_position - reference_starting_pos_[(uint64_t) reference_index]);
  SeqOrientation orientation = kForward;

  if (ret_reference_index_with_reverse != NULL)
    *ret_reference_index_with_reverse = reference_index;

  if (((uint64_t) reference_index) >= num_sequences_forward_) {
    // Relative position has to be changed, because, from the outside, it is expected that we have reverse-complemented the seed and not the reference sequence.
//    relative_pos = reference_lengths_[(uint64_t) reference_index] - relative_pos - query_length - 1 - reference_index;      // The '-1' is to compensate for the '!' character added at the end of every sequence in the data array.
    relative_pos = reference_lengths_[(uint64_t) reference_index] - relative_pos - query_length - 1;

    // Unlike BWA, we haven't reversed the order of sequences when their reverse complements were added to the index. That is why we only need to subtract the number of forward sequences, and not do (2*num_forward_sequences - gene_idx - 1).
    reference_index = reference_index - ((int64_t) num_sequences_forward_);
    orientation = kReverse;
  }

  if (ret_absolute_position != NULL)
    *ret_absolute_position = raw_position;

  if (ret_relative_position != NULL)
    *ret_relative_position = relative_pos;

  if (ret_orientation != NULL)
    *ret_orientation = orientation;

  return reference_index;
}

int64_t Index::RawPositionToReferenceIndexWithReverse(int64_t raw_position) const {
  int64_t low = 0;
  int64_t high = reference_starting_pos_.size() - 1;
  int64_t mid = 0;

  // Binary search.
  while (low <= high && low>=0 && high>=0) {
    if (low == high) {
      if (raw_position >= ((int64_t) reference_starting_pos_[low]) && raw_position < ((int64_t) (reference_starting_pos_[low] + reference_lengths_[low] + 1))) {
        return low;
      }

      break;
    }

    mid = (low + high) / 2;

    if (raw_position >= ((int64_t) reference_starting_pos_[mid]) &&
        raw_position < (((int64_t) (reference_starting_pos_[mid] + reference_lengths_[mid] + 1)))) {
      return mid;
    }
    else if (raw_position >= ((int64_t) (reference_starting_pos_[mid] + reference_lengths_[mid] + 1))) {
      // This is the case where the mid reference is entirely on the left of the raw position. We have to go right.
      low = mid + 1;

    }
    else if (raw_position < ((int64_t) reference_starting_pos_[mid])) {
      // This is the case where the mid reference is to the right of the raw position, we need to go to left.
      high = mid - 1;
    }
  }

  return -1;
}
