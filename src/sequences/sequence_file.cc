/*
 * Copyright 2014, Ivan Sovic.
 * All rights reserved.
 *
 * SequenceFile.cc
 *
 *  Created on: 15 May, 2014
 *      Author: Ivan Sovic
 */

#include "sequences/sequence_file.h"
#include "log_system/log_system.h"

//#include "bwa/kseq.h"
//KSEQ_DECLARE(gzFile)



SequenceFile::SequenceFile() {
  bwa_seq_ = NULL;
  Clear();
}

SequenceFile::SequenceFile(std::string file_path) {
  bwa_seq_ = NULL;
  Clear();
  LoadAllFromFastaOrFastq(file_path);
}

SequenceFile::SequenceFile(std::string file_path, uint64_t num_seqs_to_load) {
  bwa_seq_ = NULL;
  Clear();
  OpenFileForBatchLoading(file_path);
  LoadNextBatchNSequences(num_seqs_to_load);
}

SequenceFile::~SequenceFile() {
  Clear();
}

void SequenceFile::Clear() {
  ClearOnlyData();

  current_batch_id_ = 0;
  current_batch_starting_sequence_id_ = 0;
  open_file_path_ = "";

  if (bwa_seq_ != NULL) {
    kseq_destroy(bwa_seq_);
    bwa_seq_ = NULL;
  }
}

void SequenceFile::ClearOnlyData() {
  for (SequenceVector::iterator sequence_iterator = sequences_.begin();
      sequence_iterator != sequences_.end(); sequence_iterator++) {
    delete (*sequence_iterator);
  }

  sequences_.clear();
  current_data_size_ = 0;
}

void SequenceFile::AddSequence(SingleSequence *sequence) {
  sequences_.push_back(sequence);
  current_data_size_ += sequence->CalculateTotalSize(MEMORY_UNIT_BYTE);
}

const SequenceVector& SequenceFile::get_sequences() const {
  return sequences_;
}

uint64_t SequenceFile::GetNumberOfBases() {
  uint64_t ret = 0;
  for (uint64_t i=0; i<sequences_.size(); i++) {
    ret += sequences_[i]->get_sequence_length();
  }
  return ret;
}

void SequenceFile::set_sequences(const SequenceVector& sequences) {
  sequences_ = sequences;
}

int SequenceFile::LoadAllFromFastaOrFastq(std::string file_path, bool randomize_non_acgt_bases) {
  Clear();

  if (OpenFileForBatchLoading(file_path))
    return 1;

  int32_t l;
  uint64_t id = 0;
  SingleSequence *sequence = NULL;

  while ((l = kseq_read(bwa_seq_)) >= 0) {
    sequence = new SingleSequence();

    // BWA's parsing functions split the headers in two parts, name
    // and comment. We join them here again, but we must check that
    // these are not null pointers or empty strings.
    std::string header("");
    if (bwa_seq_->name.l > 0)
      header += std::string(bwa_seq_->name.s);
    if (bwa_seq_->comment.l > 0)
      header += std::string(" ") + std::string(bwa_seq_->comment.s);

    // If there are no quality values, do not initialize them.
    if (!bwa_seq_->qual.l) {
      sequence->InitHeaderAndDataFromAscii((char *) header.c_str(),
                                                header.length(),
                                                (int8_t *) bwa_seq_->seq.s,
                                                bwa_seq_->seq.l, id);
    } else {
      // If there are quality values, but for some reason the length of
      // the quality string is different from the length of the sequences
      // then the file is corrupted. Reporting the error here.
      if ((bwa_seq_->qual.l > 0 && bwa_seq_->seq.l != bwa_seq_->qual.l) || l == -2) {
        LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file '%s'. Sequence ID: %ld. Skipping rest of the file.", file_path.c_str(), id));
        id += 1;
        break;
      }

      sequence->InitAllFromAscii((char *) header.c_str(), header.length(),
                                      (int8_t *) bwa_seq_->seq.s,
                                      (int8_t *) bwa_seq_->qual.s, bwa_seq_->seq.l, id, id);
    }

    if (randomize_non_acgt_bases == true)
  	  sequence->RandomizeNonACGTBases();

    AddSequence(sequence);
    id += 1;

  }

  // If there are quality values, but for some reason the length of
  // the quality string is different from the length of the sequences
  // then the file is corrupted. Reporting the error here.
  // This check is repeated here because of BWA's kseq_read function
  // which can break on such sequences, thus to report the error non the less.
  // Additionally, if the functionality of kseq_read is changed, the error will
  // be reported anyway.
  if ((bwa_seq_->qual.l > 0 && bwa_seq_->seq.l != bwa_seq_->qual.l) || l == -2) {
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file '%s'. Sequence ID: %ld. Skipping rest of the file.", file_path.c_str(), id));
  }

  return CloseFileAfterBatchLoading();
}

int SequenceFile::LoadAllFromFastaOrFastqAsBatch(bool randomize_non_acgt_bases) {
  ClearOnlyData();

  current_batch_id_ = 0;
  current_batch_starting_sequence_id_ = 0;



  if (bwa_fp_ == NULL || bwa_seq_ == NULL)
  {
    if (bwa_fp_ == NULL)
      LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "Offending variable: bwt_fp_."));
    if (bwa_seq_ == NULL)
      LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Offending variable: bwt_seq_."));
    return 1;
  }

  int32_t l;
  uint64_t id = 0;
  SingleSequence *sequence = NULL;

  while ((l = kseq_read(bwa_seq_)) >= 0) {
    sequence = new SingleSequence();

    // BWA's parsing functions split the headers in two parts, name
    // and comment. We join them here again, but we must check that
    // these are not null pointers or empty strings.
    std::string header("");
    if (bwa_seq_->name.l > 0)
      header += std::string(bwa_seq_->name.s);
    if (bwa_seq_->comment.l > 0)
      header += std::string(" ") + std::string(bwa_seq_->comment.s);

    // If there are no quality values, do not initialize them.
    if (!bwa_seq_->qual.l) {
      sequence->InitHeaderAndDataFromAscii((char *) header.c_str(),
                                                header.length(),
                                                (int8_t *) bwa_seq_->seq.s,
                                                bwa_seq_->seq.l, id, id);
    } else {
      // If there are quality values, but for some reason the length of
      // the quality string is different from the length of the sequences
      // then the file is corrupted. Reporting the error here.
      if ((bwa_seq_->qual.l > 0 && bwa_seq_->seq.l != bwa_seq_->qual.l) || l == -2) {
        LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file '%s'. Sequence ID: %ld. Skipping rest of the file.", open_file_path_.c_str(), id));
        id += 1;
        break;
      }

      sequence->InitAllFromAscii((char *) header.c_str(), header.length(),
                                      (int8_t *) bwa_seq_->seq.s,
                                      (int8_t *) bwa_seq_->qual.s, bwa_seq_->seq.l, id, id);
    }

    if (randomize_non_acgt_bases == true)
  	  sequence->RandomizeNonACGTBases();

    AddSequence(sequence);
    id += 1;

  }

  // If there are quality values, but for some reason the length of
  // the quality string is different from the length of the sequences
  // then the file is corrupted. Reporting the error here.
  // This check is repeated here because of BWA's kseq_read function
  // which can break on such sequences, thus to report the error non the less.
  // Additionally, if the functionality of kseq_read is changed, the error will
  // be reported anyway.
  if ((bwa_seq_->qual.l > 0 && bwa_seq_->seq.l != bwa_seq_->qual.l) || l == -2) {
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file '%s'. Sequence ID: %ld. Skipping rest of the file.", open_file_path_.c_str(), id));
  }

  // This happens when EOF is reached. Not technically an error, but not as if the batch has been loaded.
  if (sequences_.size() == 0)
    return -1;

  return 0;
}

int SequenceFile::OpenFileForBatchLoading(std::string file_path) {
  bwa_fp_ = gzopen(file_path.c_str(), "r");

  if (bwa_fp_ == NULL) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "File path: '%s'.", file_path.c_str()));
    return 1;
  }

  bwa_seq_ = kseq_init(bwa_fp_);

  open_file_path_ = file_path;

  return 0;
}

int SequenceFile::CloseFileAfterBatchLoading() {
  if (bwa_seq_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_CLOSING_FILE, "Offending variable: bwt_seq_."));
    return 1;
  }

  kseq_destroy(bwa_seq_);

  bwa_seq_ = NULL;
  open_file_path_ = "";

  gzclose(bwa_fp_);

  return 0;
}

int SequenceFile::LoadNextBatchNSequences(uint64_t num_seqs_to_load, bool randomize_non_acgt_bases) {
  current_batch_starting_sequence_id_ += sequences_.size();

  ClearOnlyData();

  if (bwa_fp_ == NULL || bwa_seq_ == NULL)
  {
    if (bwa_fp_ == NULL)
      LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "Offending variable: bwt_fp_."));
    if (bwa_seq_ == NULL)
      LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Offending variable: bwt_seq_."));
    return 1;
  }

  int32_t l;
  uint64_t id = 0;
  uint64_t id_absolute = current_batch_starting_sequence_id_;
  SingleSequence *sequence = NULL;

  while ((l = kseq_read(bwa_seq_)) >= 0) {
    sequence = new SingleSequence();

    // BWA's parsing functions split the headers in two parts, name
    // and comment. We join them here again, but we must check that
    // these are not null pointers or empty strings.
    std::string header("");
    if (bwa_seq_->name.l > 0)
      header += std::string(bwa_seq_->name.s);
    if (bwa_seq_->comment.l > 0)
      header += std::string(" ") + std::string(bwa_seq_->comment.s);

    if (!bwa_seq_->qual.l) {  // If bwa_seq_->qual.l is equal to 0, then we are loading a FASTA file. In this case, only header and data need to be initialized.
      sequence->InitHeaderAndDataFromAscii((char *) header.c_str(),
                                                header.length(),
                                                (int8_t *) bwa_seq_->seq.s,
                                                bwa_seq_->seq.l, id, id_absolute);
    } else {  // If we got here, we are loading a FASTQ file. All three components (header, data and quality scores) need to be initialized.
      if (bwa_seq_->seq.l != bwa_seq_->qual.l) {
        LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file! Batch ID: %ld. Absolute sequence ID: %ld. Relative sequence ID: %ld. Skipping rest of the file.", current_batch_id_, (current_batch_starting_sequence_id_ + id), id));
        CloseFileAfterBatchLoading();
        return 1;
      }

      sequence->InitAllFromAscii((char *) header.c_str(), header.length(),
                                      (int8_t *) bwa_seq_->seq.s,
                                      (int8_t *) bwa_seq_->qual.s, bwa_seq_->seq.l, id, id_absolute);
    }

    if (randomize_non_acgt_bases == true)
  	  sequence->RandomizeNonACGTBases();

    AddSequence(sequence);
    id += 1;  // Increment the relative sequence id counter.
    id_absolute += 1;

    if (id >= num_seqs_to_load)  // Batch loading stopping condition.
      break;
  }

  if ((bwa_seq_->qual.l > 0 && bwa_seq_->seq.l != bwa_seq_->qual.l) || l == -2) {
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file! Batch ID: %ld. Absolute sequence ID: %ld. Relative sequence ID: %ld. Skipping rest of the file.", current_batch_id_, (current_batch_starting_sequence_id_ + id), id));
    CloseFileAfterBatchLoading();
    return 1;
  }

  // This happens when EOF is reached. Not technically an error, but not as if the batch has been loaded.
  if (sequences_.size() == 0)
    return -1;

  return 0;
}

int SequenceFile::LoadNextBatchInMegabytes(uint64_t megabytes_to_load, bool randomize_non_acgt_bases) {
  current_batch_starting_sequence_id_ += sequences_.size();

  ClearOnlyData();

  if (bwa_fp_ == NULL || bwa_seq_ == NULL)
  {
    if (bwa_fp_ == NULL)
      LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "Offending variable: bwt_fp_."));
    if (bwa_seq_ == NULL)
      LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Offending variable: bwt_seq_."));
    return 1;
  }

  int32_t l;
  uint64_t id = 0;
  uint64_t id_absolute = current_batch_starting_sequence_id_;
  SingleSequence *sequence = NULL;

  while ((l = kseq_read(bwa_seq_)) >= 0) {
    sequence = new SingleSequence();

    // BWA's parsing functions split the headers in two parts, name
    // and comment. We join them here again, but we must check that
    // these are not null pointers or empty strings.
    std::string header("");
    if (bwa_seq_->name.l > 0)
      header += std::string(bwa_seq_->name.s);
    if (bwa_seq_->comment.l > 0)
      header += std::string(" ") + std::string(bwa_seq_->comment.s);

    if (!bwa_seq_->qual.l) {  // If bwa_seq_->qual.l is equal to 0, then we are loading a FASTA file. In this case, only header and data need to be initialized.
      sequence->InitHeaderAndDataFromAscii((char *) header.c_str(),
                                                header.length(),
                                                (int8_t *) bwa_seq_->seq.s,
                                                bwa_seq_->seq.l, id, id_absolute);
    } else {  // If we got here, we are loading a FASTQ file. All three components (header, data and quality scores) need to be initialized.
      if (bwa_seq_->seq.l != bwa_seq_->qual.l) {
        LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file! Batch ID: %ld. Absolute sequence ID: %ld. Relative sequence ID: %ld. Skipping rest of the file.", current_batch_id_, (current_batch_starting_sequence_id_ + id), id));
        CloseFileAfterBatchLoading();
        return 1;
      }

      sequence->InitAllFromAscii((char *) header.c_str(), header.length(),
                                      (int8_t *) bwa_seq_->seq.s,
                                      (int8_t *) bwa_seq_->qual.s, bwa_seq_->seq.l, id, id_absolute);
    }

    if (randomize_non_acgt_bases == true)
  	  sequence->RandomizeNonACGTBases();

    AddSequence(sequence);
    id += 1;  // Increment the relative sequence id counter.
    id_absolute += 1;

    if (ConvertFromBytes(MEMORY_UNIT_MEGABYTE, current_data_size_) >= megabytes_to_load)  // Batch loading stopping condition.
      break;
  }

  if ((bwa_seq_->qual.l > 0 && bwa_seq_->seq.l != bwa_seq_->qual.l) || l == -2) {
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file! Batch ID: %ld. Absolute sequence ID: %ld. Relative sequence ID: %ld. Skipping rest of the file.", current_batch_id_, (current_batch_starting_sequence_id_ + id), id));
    CloseFileAfterBatchLoading();
    return 1;
  }

  // This happens when EOF is reached. Not technically an error, but not as if the batch has been loaded.
  if (sequences_.size() == 0)
    return -1;

  return 0;
}

uint64_t SequenceFile::CalculateTotalSize(int32_t memory_unit) {
  uint64_t total_size = 0;

  for (SequenceVector::iterator sequence_iterator = sequences_.begin();
      sequence_iterator != sequences_.end(); sequence_iterator++) {
    total_size += (*sequence_iterator)->CalculateTotalSize(MEMORY_UNIT_BYTE);
  }

  if (memory_unit == MEMORY_UNIT_BYTE)
    total_size = total_size / ((uint64_t) 1);
  else if (memory_unit == MEMORY_UNIT_KILOBYTE)
    total_size = total_size / ((uint64_t) 1024);
  else if (memory_unit == MEMORY_UNIT_MEGABYTE)
    total_size = total_size / (((uint64_t) 1024) * ((uint64_t) 1024));
  else if (memory_unit == MEMORY_UNIT_GIGABYTE)
    total_size = total_size / (((uint64_t) 1024) * ((uint64_t) 1024) * ((uint64_t) 1024));
  else
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Memory unit not recognized! Returning value in bytes."));

  return total_size;
}

void SequenceFile::Verbose(FILE *fp) {
  fprintf (fp, "Num sequences: %ld\n", sequences_.size());
  fprintf (fp, "Currently open file (only when batch loading): '%s'\n", open_file_path_.c_str());
  fprintf (fp, "Current batch ID: %ld\n", current_batch_id_);
  fprintf (fp, "Current batch starting sequence ID: %ld\n", current_batch_starting_sequence_id_);

  fprintf (fp, "\n");

  for (SequenceVector::iterator sequence_iterator = sequences_.begin();
      sequence_iterator != sequences_.end(); sequence_iterator++) {
    (*sequence_iterator)->Verbose(fp);
    fprintf (fp, "\n");
  }

  fflush(fp);
}
