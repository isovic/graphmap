/*
 * path_graph_entry.cc
 *
 *  Created on: Oct 16, 2014
 *      Author: ivan
 */

#include "containers/path_graph_entry.h"



//PathGraphEntry::PathGraphEntry() {
////  local_scores_id = -1;
////  region_id = -1;
////  source_node = -1;
////  sink_node = -1;
////  lcs_length = -1;
////  covered_bases = -1;
////  region_kmer_count = -1;
////  query.start = query.end = 0;
////  reference.start = reference.end = 0;
////  distance.start = distance.end = 0;
////  ratio = 0.0f;
////  ratio_suppress = 100.0f;
////  deviation = 0.0f;
//////  lcsl1_entries.Clear();
////
////  l1_l = 0;
////  l1_k = 0.0;
////  l1_lmin = 0;
////  l1_lmax = 0;
////  l1_confidence_abs = 0;
////  l1_std = 0;
////  l1_rough_start = 0;
////  l1_rough_end = 0;
////
////  edit_distance = -1;
////  edit_distance_end_pos = -1;
////
////  fpfilter = 0.0f;
////  fpfilter_score_std = 0.0f;
////  fpfilter_score_cov_bases = 0.0f;
////  fpfilter_score_read_len = 0.0f;
////  fpfilter_score_query_len = 0.0f;
////  num_covering_kmers = 0;
////
////  num_eq_ops = 0;
////  num_x_ops = 0;
////  num_i_ops = 0;
////  num_d_ops = 0;
//  read_ = NULL;
//  index_ = NULL;
//  parameters_ = NULL;
//}

PathGraphEntry::~PathGraphEntry() {
  alignments_secondary_.clear();
}

//PathGraphEntry::PathGraphEntry(Index *index, SingleSequence *read, ProgramParameters *parameters,
//                               int64_t lcs_length, int64_t cov_bases_query, int64_t cov_bases_ref, float deviation,
//                               int64_t query_start, int64_t query_end, int64_t reference_start, int64_t reference_end) {

//  Set(local_scores_id, region_id, source_node, sink_node, lcs_length, covered_bases, region_kmer_count, deviation, query_start, query_end, reference_start, reference_end, actual_query_start,
//      actual_query_end, actual_reference_start, actual_reference_end);
PathGraphEntry::PathGraphEntry(const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const Region &region,
                               MappingResults *mapping_data, L1Results *l1_data, AlignmentResults *alignment_data) {
  fpfilter_ = 0.0f;
  fpfilter_cov_bases_ = 0.0f;
  fpfilter_query_len_ = 0.0f;
  fpfilter_std_ = 0.0f;
  fpfilter_read_len_ = 0.0f;

  Set(index, read, parameters, region, mapping_data, l1_data, alignment_data);
}

void PathGraphEntry::Set(const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const Region &region,
                         MappingResults *mapping_data, L1Results *l1_data, AlignmentResults *alignment_data) {
  index_ = index;
  read_ = read;
  parameters_ = parameters;
  region_info_ = region;

  if (mapping_data != NULL)
    set_mapping_data((*mapping_data));
  if (l1_data != NULL)
    set_l1_data((*l1_data));
  if (alignment_data != NULL)
    AddSecondaryAlignmentData((*alignment_data));
}

void PathGraphEntry::Verbose(FILE *fp) const {
  fprintf (fp, "%s", VerboseToString().c_str());
}

std::string PathGraphEntry::VerboseInfoToString(std::string delimiter) const {
  std::stringstream ss;

  ss << "L1 info:\n";
  ss << "\t° l1_l = " << l1_info_.l1_l << ", l1_k = " << l1_info_.l1_k << "\n";
  ss << "\t° l1_lmin = " << l1_info_.l1_lmin << ", l1_lmax = " << l1_info_.l1_lmax << "\n";
  ss << "\t° l1_confidence_abs = " << l1_info_.l1_confidence_abs << ", l1_std = " << l1_info_.l1_std << "\n";
  ss << "\t° l1_rough_start = " << l1_info_.l1_rough_start << ", l1_rough_end = " << l1_info_.l1_rough_end;

  return ss.str();
}

std::string PathGraphEntry::VerboseToString(std::string delimiter) const {
  std::stringstream ss;

  std::string actual_delimiter = delimiter;
  if (delimiter == "\r")
    actual_delimiter = " ";

  int64_t query_start = mapping_info_.query_coords.start;
  int64_t query_end = (mapping_info_.query_coords.end);
  int64_t read_length = read_->get_sequence_length();
  int64_t db_size = index_->get_data_length_forward();
  int64_t cov_bases = mapping_info_.cov_bases_max;
  float error_rate = parameters_->error_rate;
  float std = mapping_info_.deviation;
  float ratio = CalcDistanceRatio();
  float ratio_suppress = CalcDistanceRatioSuppress();
  int64_t reference_start = index_->get_reference_starting_pos()[region_info_.reference_id];

  ss << "lcs_length=" << mapping_info_.lcs_length << actual_delimiter <<
        "cov_bases=" << mapping_info_.cov_bases_max << actual_delimiter <<
        "num_kmers=" << mapping_info_.num_covering_kmers << actual_delimiter <<
        "act_q[" << mapping_info_.query_coords.start << "," << mapping_info_.query_coords.end << "]" << actual_delimiter <<
        "act_r_l[" << (mapping_info_.ref_coords.start - reference_start) << "," << (mapping_info_.ref_coords.end - reference_start) << "]" << actual_delimiter <<
        "act_r[" << mapping_info_.ref_coords.start << "," << mapping_info_.ref_coords.end << "]" << actual_delimiter <<
        "act_d[" << (mapping_info_.query_coords.end - mapping_info_.query_coords.start) << "," << (mapping_info_.ref_coords.end - mapping_info_.ref_coords.start) << "]" << actual_delimiter;
  if (delimiter == "\r")
    ss << "\n        ";
  ss << "p[" << actual_delimiter << ratio << "]" << actual_delimiter <<
        "supp[" << ratio_suppress << "]" << actual_delimiter <<
        "std[" << mapping_info_.deviation << "]" << actual_delimiter <<
        "AS[" << FormatString("%.4f", CalcFPFilter()) << "]" << actual_delimiter <<
        "AS_std[" << FormatString("%.4f", CalcFPFactorStd_(std, read_length, error_rate)) << "]" << actual_delimiter <<
        "AS_cov_bases[" << FormatString("%.4f", CalcFPFactorCovBases_(cov_bases, query_start, query_end, error_rate)) << "]" << actual_delimiter <<
        "AS_read_len[" << FormatString("%.4f", CalcFPFactorReadLength_(read_length, db_size)) << "]" << actual_delimiter <<
        "AS_query_len[" << FormatString("%.4f", CalcFPFactorQueryLength_(query_start, query_end, read_length)) << "]" << actual_delimiter <<
        "is_mapped=" << ((mapping_info_.is_mapped == true) ? "true" : "false") << actual_delimiter <<
        "is_reverse=" << ((mapping_info_.is_reverse == true) ? "true" : "false");

  return ss.str();
}


float PathGraphEntry::CalcDistanceRatio() const {
  int64_t distance_query = mapping_info_.query_coords.end - mapping_info_.query_coords.start;
  int64_t distance_ref = mapping_info_.ref_coords.end - mapping_info_.ref_coords.start;
  float distance_ratio = 1.0f;
  if (distance_query != 0)
    distance_ratio = ((float) std::min(distance_query, distance_ref)) / ((float) std::max(distance_query, distance_ref));
  return distance_ratio;
}

float PathGraphEntry::CalcDistanceRatioSuppress() const {
  int64_t distance_query = mapping_info_.query_coords.end - mapping_info_.query_coords.start;
  int64_t distance_ref = mapping_info_.ref_coords.end - mapping_info_.ref_coords.start;
  float distance_ratio = 1.0f;
  if (distance_query != 0)
    distance_ratio = ((float) std::min(distance_query, distance_ref)) / ((float) std::max(distance_query, distance_ref));

  return (distance_ratio < 1.0f) ? (1.0f - distance_ratio) : (distance_ratio - 1.0f);
}

float PathGraphEntry::CalcFPFilter() const {
  int64_t query_start = mapping_info_.query_coords.start;
  int64_t query_end = (mapping_info_.query_coords.end);
  int64_t read_length = read_->get_sequence_length();
  int64_t db_size = index_->get_data_length_forward();
  int64_t cov_bases = mapping_info_.cov_bases_max;
  float error_rate = parameters_->error_rate;
  float std = mapping_info_.deviation;

  float fpfilter =  CalcFPFactorStd_(std, read_length, error_rate) *
                    CalcFPFactorQueryLength_(query_start, query_end, read_length) *
                    CalcFPFactorCovBases_(cov_bases, query_start, query_end, error_rate)  *
                    CalcFPFactorReadLength_(read_length, db_size);
  fpfilter = std::min(((float) 1.0f), fpfilter);
  fpfilter = std::max(((float) 0.0f), fpfilter);
  return fpfilter;
}

float PathGraphEntry::CalcFPFilterStatic(const MappingResults &mapping_info, const Index *index, const SingleSequence *read, const ProgramParameters *parameters) {
  int64_t query_start = mapping_info.query_coords.start;
  int64_t query_end = (mapping_info.query_coords.end + parameters->k_graph);
  int64_t read_length = read->get_sequence_length();
  int64_t db_size = index->get_data_length_forward();
  int64_t cov_bases = mapping_info.cov_bases_max;
  float error_rate = parameters->error_rate;
  float std = mapping_info.deviation;

  float fpfilter =  CalcFPFactorStd_(std, read_length, error_rate) *
                    CalcFPFactorQueryLength_(query_start, query_end, read_length) *
                    CalcFPFactorCovBases_(cov_bases, query_start, query_end, error_rate)  *
                    CalcFPFactorReadLength_(read_length, db_size);
  fpfilter = std::min(((float) 1.0f), fpfilter);
  fpfilter = std::max(((float) 0.0f), fpfilter);
  return fpfilter;
}

float PathGraphEntry::CalcFPFactorReadLength_(int64_t read_length, int64_t db_size) {
  float factor_read_length = 0.0f;
  if (read_length > 0 && db_size > 0)
    factor_read_length = log(read_length) / log(db_size);
  return factor_read_length;
}

float PathGraphEntry::CalcFPFactorQueryLength_(int64_t query_start, int64_t query_end, int64_t read_length) {
  float factor_query_length = 1.0f;
  if (read_length > 0)
    factor_query_length = ((float) (query_end - query_start)) / ((float) read_length);
  factor_query_length = std::min(((float) 1.0f), factor_query_length);
  factor_query_length = std::max(((float) 0.0f), factor_query_length);
  return factor_query_length;
}

float PathGraphEntry::CalcFPFactorCovBases_(int64_t covered_bases, int64_t query_start, int64_t query_end, float error_rate) {
  float factor_covered_bases = 1.0f;
  if (query_end != query_start && error_rate < 1.0f)
    factor_covered_bases = (((float) covered_bases) / ((float) (query_end - query_start) * (1.0f - error_rate)));
  factor_covered_bases = std::min(((float) 1.0f), factor_covered_bases);
  factor_covered_bases = std::max(((float) 0.0f), factor_covered_bases);
  return factor_covered_bases;
}

float PathGraphEntry::get_fpfilter() const {
  return fpfilter_;
}

void PathGraphEntry::set_fpfilter(float fpfilter) {
  fpfilter_ = fpfilter;
}

const L1Results& PathGraphEntry::get_l1_data() const {
  return l1_info_;
}

void PathGraphEntry::set_l1_data(L1Results& l1Data) {
  l1_info_ = l1Data;
}

const MappingResults& PathGraphEntry::get_mapping_data() const {
  return mapping_info_;
}

void PathGraphEntry::set_mapping_data(MappingResults& mappingData) {
  mapping_info_ = mappingData;

  int64_t query_start = mapping_info_.query_coords.start;
  int64_t query_end = (mapping_info_.query_coords.end);
  int64_t read_length = read_->get_sequence_length();
  int64_t db_size = index_->get_data_length_forward();
  int64_t cov_bases = mapping_info_.cov_bases_max;
  float error_rate = parameters_->error_rate;
  float std = mapping_info_.deviation;

  fpfilter_cov_bases_ = CalcFPFactorCovBases_(cov_bases, query_start, query_end, error_rate);
  fpfilter_query_len_ = CalcFPFactorQueryLength_(query_start, query_end, read_length);
  fpfilter_std_ = CalcFPFactorStd_(std, read_length, error_rate);
  fpfilter_read_len_ = CalcFPFactorReadLength_(read_length, db_size);
  fpfilter_ = CalcFPFilter();
}

const Region& PathGraphEntry::get_region_data() const {
  return region_info_;
}

void PathGraphEntry::set_region_data(Region& regionData) {
  region_info_ = regionData;
}

void PathGraphEntry::AddSecondaryAlignmentData(AlignmentResults alignment_info) {
  alignments_secondary_.push_back(alignment_info);
}

std::string PathGraphEntry::GenerateSAMFromInfoAlignment_(const AlignmentResults &alignment_info, const MappingMetadata &mapping_metadata, bool is_primary, int64_t verbose_sam_output) const {
  std::stringstream ss;

  std::string qname = ((std::string) (read_->get_header()));
  std::string qname_for_out = (verbose_sam_output < 4) ? (TrimToFirstSpace(qname)) : (qname);
  std::string rname = region_info_.rname;
  std::string rname_for_out = (verbose_sam_output < 4) ? (TrimToFirstSpace(rname)) : (rname);

  uint32_t reverse = (alignment_info.is_reverse == false) ? 0 : SAM_THIS_SEG_REVERSED;
  uint32_t mapped = (alignment_info.is_aligned) ? 0 : SAM_THIS_SEG_UNMAPPED;  // This means that the segment is mapped.
  uint32_t secondary_alignment = ((is_primary == true) ? (0) : (SAM_SECONDARY_ALIGNMENT));
  uint32_t flag = reverse | mapped | secondary_alignment;

  std::string rnext = "*";
  int64_t pnext = 0;
  int64_t tlen = 0;

  ss << qname_for_out << "\t";
  ss << flag << "\t";
  ss << rname_for_out << "\t";
  ss << alignment_info.pos_start << "\t";
  ss << ((int64_t) alignment_info.mapping_quality) << "\t";      // To avoid confusion with the definition of mapping quality, we will use the value of 255, and report the actual quality as AS optional parameter.
  ss << alignment_info.cigar << "\t";
  ss << rnext << "\t";
  ss << pnext << "\t";
  ss << tlen << "\t";

  if (alignment_info.is_reverse == false) {
    ss << read_->get_data() << "\t";

    if (verbose_sam_output < 5 && read_->get_quality_length() > 0 && read_->get_quality() != NULL) {
      ss << read_->get_quality() << "\t";
    } else {
      ss << "*" << "\t";
    }

  } else {
    ss << read_->GetReverseComplementAsString() << "\t";

    if (verbose_sam_output < 5 && read_->get_quality_length() > 0 && read_->get_quality() != NULL) {
      ss << read_->GetReverseQualityAsString() << "\t";
    } else {
      ss << "*" << "\t";
    }
  }

  std::stringstream ss_optional1;
  ss_optional1 << "NM:i:" << alignment_info.edit_distance << "\t"; // Specified by SAM format.
  ss_optional1 << "AS:i:" << alignment_info.alignment_score << "\t";
  ss_optional1 << "H0:i:" << alignment_info.num_secondary_alns << "\t"; // Specified by SAM format.
  ss_optional1 << "ZE:f:" << alignment_info.evalue << "\t";
  ss_optional1 << "ZF:f:" << fpfilter_ << "\t";
  ss_optional1 << "ZQ:i:" << read_->get_sequence_length() << "\t";
  ss_optional1 << "ZR:i:" << index_->get_reference_lengths()[region_info_.reference_id];

  ss << ss_optional1.str();

  if (verbose_sam_output >= 3) {
    float mismatch_rate = (((float) (alignment_info.num_x_ops + alignment_info.num_i_ops + alignment_info.num_d_ops)) / ((float) (alignment_info.num_eq_ops + alignment_info.num_x_ops + alignment_info.num_d_ops + alignment_info.num_i_ops)));
    float match_rate = (((float) alignment_info.num_eq_ops) / ((float) alignment_info.nonclipped_length)); // ((float) read_->get_sequence_length()));

    std::stringstream X3_ss;
    X3_ss << VerboseToString("_");
    X3_ss << "__region_votes=" << get_region_data().region_votes; // << "__max_region_votes=" << get_region_data()..bin_value << "__num_region_iterations=" << num_region_iterations_;
    X3_ss << "__num_eq_ops=" << alignment_info.num_eq_ops << "__num_x_ops=" << alignment_info.num_x_ops << "__num_i_ops=" << alignment_info.num_i_ops << "__num_d_ops=" << alignment_info.num_d_ops << "__nonclippedlen=" << alignment_info.nonclipped_length << "__match_rate=" << FormatString("%.2f", match_rate) << "__mismatch_rate=" << FormatString("%.2f", mismatch_rate);

    ss << "\tX3:Z:" << X3_ss.str();

    std::stringstream X4_ss;
    X4_ss << "Timings(sec)__regionselection=" << mapping_metadata.stats_time_region_selection << "__mapping=" << mapping_metadata.stats_time_mapping << "__alignment=" << mapping_metadata.stats_time_alignment;
    ss << "\tX4:Z:" << X4_ss.str();
  }

  return ss.str();
}

std::string PathGraphEntry::GenerateSAM(bool is_primary, int64_t verbose_sam_output) const {
  std::stringstream ss;

  ss << GenerateSAMFromInfoAlignment_(alignment_primary_, mapping_metadata_, is_primary, verbose_sam_output);

  for (int64_t i=0; i<alignments_secondary_.size(); i++) {
    ss << "\n";
    ss << GenerateSAMFromInfoAlignment_(alignments_secondary_[i], mapping_metadata_, false, verbose_sam_output);
  }

//  uint32_t reverse = (orientation == kForward) ? 0 : SAM_THIS_SEG_REVERSED;
//  uint32_t mapped = 0;  // This means that the segment is mapped.
//  uint32_t secondary_alignment = ((i == 0) ? (0) : (SAM_SECONDARY_ALIGNMENT));
//  uint32_t flag = reverse | mapped | secondary_alignment;
//
//  std::stringstream X1_ss;
//  X1_ss << metagen_alignment_score_;
//  std::stringstream X2_ss;
//  X2_ss << read_->get_sequence_length();
//  std::stringstream X3_ss;
//  X3_ss << final_mapping_ptrs_.at(i)->VerboseToString("_");
//  X3_ss << "__region_votes=" << final_mapping_ptrs_.at(i)->get_region_data().region_votes << "__max_region_votes=" << bins_.front().bin_value << "__num_region_iterations=" << num_region_iterations_;
////    X3_ss << "__num_eq_ops=" << num_eq_ops << "__num_x_ops=" << num_x_ops << "__num_i_ops=" << num_i_ops << "__num_d_ops=" << num_d_ops << "__match_rate=" << FormatString("%.2f", (((float) num_eq_ops) / ((float) read_->get_sequence_length()))) << "__mismatch_rate=" << FormatString("%.2f", (((float) (num_x_ops + num_i_ops + num_d_ops)) / ((float) (num_eq_ops + num_x_ops + num_d_ops))));
//  X3_ss << "__num_eq_ops=" << num_eq_ops << "__num_x_ops=" << num_x_ops << "__num_i_ops=" << num_i_ops << "__num_d_ops=" << num_d_ops << "__match_rate=" << FormatString("%.2f", match_rate) << "__mismatch_rate=" << FormatString("%.2f", mismatch_rate);
//  std::stringstream X5_ss;
////    X5_ss << fp_filter;
//  X5_ss << final_mapping_ptrs_.at(i)->fpfilter;
//
//  final_mapping_ptrs_.at(i)->sam.qname = std::string(read_->get_header());
//  final_mapping_ptrs_.at(i)->sam.flag = flag;
//  final_mapping_ptrs_.at(i)->sam.rname = index_->get_headers()[reference_id];
//  final_mapping_ptrs_.at(i)->sam.pos = relative_position_left_part + 1;
//  final_mapping_ptrs_.at(i)->sam.mapq = mapping_quality_;
//  final_mapping_ptrs_.at(i)->sam.cigar = (orientation == kForward) ? (cigar_left_part) : (ReverseCigarString(cigar_left_part));
//  final_mapping_ptrs_.at(i)->sam.seq = (orientation == kForward) ? (std::string((char *) read_->get_data())) : (read_->GetReverseComplementAsString());
//
//  final_mapping_ptrs_.at(i)->sam.NM = edit_distance;
//  final_mapping_ptrs_.at(i)->sam.H0 = num_same_mappings_;
//  final_mapping_ptrs_.at(i)->sam.AS = AS_left_part;
//  final_mapping_ptrs_.at(i)->sam.X1 = X1_ss.str();
//  final_mapping_ptrs_.at(i)->sam.X2 = X2_ss.str();
//  final_mapping_ptrs_.at(i)->sam.X3 = X3_ss.str();
//  final_mapping_ptrs_.at(i)->sam.evalue = evalue_left_part;
//
//  final_mapping_ptrs_.at(i)->num_eq_ops = num_eq_ops;
//  final_mapping_ptrs_.at(i)->num_x_ops = num_x_ops;
//  final_mapping_ptrs_.at(i)->num_i_ops = num_i_ops;
//  final_mapping_ptrs_.at(i)->num_d_ops = num_d_ops;
//  final_mapping_ptrs_.at(i)->sam.cigar_split_part = "";
//  final_mapping_ptrs_.at(i)->sam.pos_split_part = 0;
//  final_mapping_ptrs_.at(i)->sam.AS_split_part = 0;
//
////    final_mapping_ptrs_.at(i)->sam.ZF = X5_ss.str();
//  final_mapping_ptrs_.at(i)->sam.fpfilter = final_mapping_ptrs_.at(i)->fpfilter;
//
//  if (cigar_right_part.size() > 0) {
//    final_mapping_ptrs_.at(i)->sam.cigar_split_part = (orientation == kForward) ? (cigar_right_part) : (ReverseCigarString(cigar_right_part));
//    final_mapping_ptrs_.at(i)->sam.pos_split_part = relative_position_right_part + 1;
//    final_mapping_ptrs_.at(i)->sam.AS_split_part = AS_right_part;
//    final_mapping_ptrs_.at(i)->sam.evalue_split_part = evalue_right_part;
//    final_mapping_ptrs_.at(i)->sam.fpfilter_split_part = final_mapping_ptrs_.at(i)->fpfilter;
//  }
//
//  if (read_->get_quality_length() > 0) {
//    if (orientation == kForward) {
////        printf ("Forward!\n");
////        fflush(stdout);
//      if (read_->get_quality() != NULL)
//        final_mapping_ptrs_.at(i)->sam.qual = std::string((char *) read_->get_quality());  // GetSubstring((char *) read_->get_quality(), read_->get_quality_length());
//      else
//        final_mapping_ptrs_.at(i)->sam.qual = std::string("*");
//    } else {
////        printf ("Reverse!\n");
////        fflush(stdout);
//      std::string reverse_quality = read_->GetReverseQualityAsString();
//      if (reverse_quality.size() > 0)
//        final_mapping_ptrs_.at(i)->sam.qual = (reverse_quality);
//      else
//        final_mapping_ptrs_.at(i)->sam.qual = std::string("*");
//    }
//
////      final_mapping_ptrs_.at(i)->sam.qual = (orientation == kForward) ? (std::string((char *) read_->get_quality())) : (read_->GetReverseQualityAsString());
//  }

  return ss.str();
}

std::string PathGraphEntry::GenerateAMOSFromInfoMappping(const MappingResults &mapping_info) const {
  std::stringstream ss;

  int64_t l = l1_info_.l1_l - ((int64_t) index_->get_reference_starting_pos()[region_info_.reference_id]);
  int64_t alignment_end = l1_info_.l1_k * ((float) read_->get_sequence_length()) + l;

  std::string adj = (mapping_info.is_reverse == false) ? OVERLAP_NORMAL : OVERLAP_INNIE;
  int64_t read1_id = read_->get_sequence_id() + 1;
  int64_t read2_id = (region_info_.reference_id % index_->get_num_sequences_forward()) + 1;
  int64_t score = (mapping_info.query_coords.end - mapping_info.query_coords.start);
  // The L1 line is the projection of the median diagonal to the reference. In this case, the reference is another read.
  // If the L1 line starts before 0, then read1 has an overhang. The overhang is positive if the prefix of the read is not
  // overlapping, but positive if the suffix is not overlapping, so the l value needs to be inverted.
  int64_t ahang = 0 - l;
  int64_t bhang = index_->get_reference_lengths()[region_info_.reference_id] - alignment_end;

//  printf ("%s\n", VerboseToString().c_str());
//  printf ("ref_length = %ld\n", index_->get_reference_lengths()[region_info_.reference_id]);
//  printf ("\n");
//  fflush(stdout);

  ss << "{OVL" << "\n";
  ss << "adj:" << adj << "\n";
  ss << "rds:" << read1_id << "," << read2_id << "\n";
  ss << "scr:" << score << "\n";
  ss << "ahg:" << ahang << "\n";
  ss << "bhg:" << bhang << "\n";
  ss << "}";

  return ss.str();
}

std::string PathGraphEntry::GenerateAMOS() const {
  std::stringstream ss;

  ss << GenerateAMOSFromInfoAlignment_(alignment_primary_);

  for (int64_t i=0; i<alignments_secondary_.size(); i++) {
    ss << "\n";
    ss << GenerateAMOSFromInfoAlignment_(alignments_secondary_[i]);
  }

  return ss.str();
}

std::string PathGraphEntry::GenerateAMOSFromInfoAlignment_(const AlignmentResults &alignment_info) const {
  std::stringstream ss;

//  int64_t l = l1_info_.l1_l - ((int64_t) index_->get_reference_starting_pos()[region_info_.reference_id]);
//  int64_t alignment_end = l1_info_.l1_k * ((float) read_->get_sequence_length()) + l;

  int64_t alignment_start = alignment_info.pos_start - 1; /// -1 because the variable holds the 1-offset position for alignment as is supposed to be in the SAM format.
  int64_t alignment_end = alignment_start + alignment_info.nonclipped_length;

  char clip_op_front=0, clip_op_back=0;
  int64_t clip_count_front=0, clip_count_back=0;
  if (GetClippingOpsFromCigar(alignment_info.cigar, &clip_op_front, &clip_count_front, &clip_op_back, &clip_count_back)) {
    return std::string("");
  }

  std::string adj = (alignment_info.is_reverse == false) ? OVERLAP_NORMAL : OVERLAP_INNIE;
  int64_t read1_id = read_->get_sequence_id() + 1;
  int64_t read2_id = (region_info_.reference_id % index_->get_num_sequences_forward()) + 1;
  int64_t score = alignment_info.alignment_score;

  // The L1 line is the projection of the median diagonal to the reference. In this case, the reference is another read.
  // If the L1 line starts before 0, then read1 has an overhang. The overhang is positive if the prefix of the read is not
  // overlapping, but positive if the suffix is not overlapping, so the l value needs to be inverted.
//  int64_t ahang = 0 - l;
//  int64_t bhang = index_->get_reference_lengths()[region_info_.reference_id] - alignment_end;

  //ahg - Ahang. Length of the non-overlapping portion of the first read.
  //bhg - Bhang. Length of the non-overlapping portion of the second read.
  /// The position (alignment_start - clip_count_front) would be roughly where alignment of one reads starts on another.
  /// If this value is > 0, the first read starts within read2, and thus ahang needs to be negative (hence the '-' sign).
  int64_t ahang = - (alignment_start - clip_count_front);
  int64_t bhang = index_->get_reference_lengths()[region_info_.reference_id] - (alignment_end + clip_count_back);

  ss << "{OVL" << "\n";
  ss << "adj:" << adj << "\n";
  ss << "rds:" << read1_id << "," << read2_id << "\n";
  ss << "scr:" << score << "\n";
  ss << "ahg:" << ahang << "\n";
  ss << "bhg:" << bhang << "\n";
  ss << "}";

  return ss.str();
}

float PathGraphEntry::get_fpfilter_cov_bases() {
  return fpfilter_cov_bases_;
}

void PathGraphEntry::set_fpfilter_cov_bases(float fpfilterCovBases) {
  fpfilter_cov_bases_ = fpfilterCovBases;
}

float PathGraphEntry::get_fpfilter_query_len() {
  return fpfilter_query_len_;
}

void PathGraphEntry::set_fpfilter_query_len(float fpfilterQueryLen) {
  fpfilter_query_len_ = fpfilterQueryLen;
}

float PathGraphEntry::get_fpfilter_read_len() {
  return fpfilter_read_len_;
}

void PathGraphEntry::set_fpfilter_read_len(float fpfilterReadLen) {
  fpfilter_read_len_ = fpfilterReadLen;
}

float PathGraphEntry::get_fpfilter_std() {
  return fpfilter_std_;
}

void PathGraphEntry::set_fpfilter_std(float fpfilterStd) {
  fpfilter_std_ = fpfilterStd;
}

AlignmentResults& PathGraphEntry::get_alignment_primary() {
  return alignment_primary_;
}

void PathGraphEntry::set_alignment_primary(AlignmentResults& alignmentPrimary) {
  alignment_primary_ = alignmentPrimary;
}

std::vector<AlignmentResults>& PathGraphEntry::get_alignments_secondary() {
  return alignments_secondary_;
}

void PathGraphEntry::set_alignments_secondary(std::vector<AlignmentResults>& alignmentsSecondary) {
  alignments_secondary_ = alignmentsSecondary;
}

MappingMetadata& PathGraphEntry::get_mapping_metadata() {
  return mapping_metadata_;
}

void PathGraphEntry::set_mapping_metadata(MappingMetadata& mappingMetadata) {
  mapping_metadata_ = mappingMetadata;
}

float PathGraphEntry::CalcFPFactorStd_(float std, int64_t read_length, float error_rate) {
  float factor_std = 1.0f;
  if (error_rate > 0.0f && read_length > 0)
    factor_std = std::max(0.0f, (float) (1.0f - (std) / (error_rate / sqrt(2.0f) * read_length)));
  return factor_std;
}

//void PathGraphEntry::GenerateFPFilter(Index *index, SingleSequence *read_, ProgramParameters &parameters_) {
//  int64_t read_length = read->get_sequence_length();
//
//    float factor_std = 1.0f;
//    if (parameters.error_rate > 0.0f && read_length > 0)
//      factor_std = std::max(0.0f, (float) (1.0f - (deviation) / (parameters.error_rate / sqrt(2.0f) * read_length)));
//
//    float factor_query_coverage = 1.0f;
//    if (read_length > 0)
//      factor_query_coverage = ((float) (actual_query.end - actual_query.start)) / ((float) read_length);
//    factor_query_coverage = std::min(((float) 1.0f), factor_query_coverage);
//    factor_query_coverage = std::max(((float) 0.0f), factor_query_coverage);
//
//    float factor_covered_bases = 1.0f;
//    if (((float) (actual_query.end - actual_query.start + parameters.k_graph) * (1.0f - parameters.error_rate)) != 0.0f)
//      factor_covered_bases = (((float) covered_bases) / ((float) (actual_query.end - actual_query.start + parameters.k_graph) * (1.0f - parameters.error_rate)));
//    factor_covered_bases = std::min(((float) 1.0f), factor_covered_bases);
//    factor_covered_bases = std::max(((float) 0.0f), factor_covered_bases);
//
//    float factor_read_length = 0.0f;
//    if (read->get_sequence_length() > 0 && index->get_data_length_forward() > 0)
//      factor_read_length = log(read->get_sequence_length()) / log(index->get_data_length_forward());
//
//    this->fpfilter = factor_std * factor_query_coverage * factor_covered_bases* factor_read_length; // * (1.0f - ratio_suppress);
//    this->fpfilter_score_std = factor_std;
//    this->fpfilter_score_cov_bases = factor_covered_bases;
//    this->fpfilter_score_read_len = factor_read_length;
//    this->fpfilter_score_query_len = factor_query_coverage;
//
//    this->fpfilter = std::min(((float) 1.0f), this->fpfilter);
//    this->fpfilter = std::max(((float) 0.0f), this->fpfilter);
//}
//
//bool PathGraphEntry::IsLessThan(PathGraphEntry &op1)
//{
//  if (this->covered_bases < op1.covered_bases)
//    return true;
//  else if (this->covered_bases > op1.covered_bases)
//    return false;
//  else
//    return (this->ratio_suppress > op1.ratio_suppress);
//}
