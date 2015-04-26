/*
 * path_graph_entry.h
 *
 *  Created on: Oct 16, 2014
 *      Author: ivan
 */

#ifndef PATH_GRAPH_ENTRY_H_
#define PATH_GRAPH_ENTRY_H_

#include <stdio.h>
#include <string>
#include <sstream>
#include "containers/range.h"
#include "containers/region.h"
#include "containers/score_registry.h"
#include "utility/utility_general.h"
#include "utility/program_parameters.h"
//#include "sam/sam_entry.h"
#include "containers/vertices.h"

// N(ormal) would indicate both sequences are in the same orientation.
#define OVERLAP_NORMAL          "N"
// I(nnie) indicates the second sequence is in the opposite orientation than the first
#define OVERLAP_INNIE           "I"
#define OVERLAP_OUTIE           "O"
#define OVERLAP_ANTI_NORMAL     "A"



typedef struct InfoMapping {
  int64_t lcs_length = 0;
  int64_t cov_bases_max = 0;
  int64_t cov_bases_query = 0;
  int64_t cov_bases_ref = 0;
  int64_t num_covering_kmers = 0;
  float deviation = 0.0f;
  Range query_coords;
  Range ref_coords;
  bool is_mapped = false;
  bool is_reverse = false;
  int64_t local_score_id = 0;
//  Range distance;
//  float ratio;
//  float ratio_suppress;

//  Range actual_query;
//  Range actual_reference;
//  int64_t query_id;
//  int64_t query_length;
//  std::string qname;
//  std::string rname;      // Ovo dodati u Region!

//  float fpfilter;
//  float fpfilter_score_std;
//  float fpfilter_score_cov_bases;
//  float fpfilter_score_read_len;
//  float fpfilter_score_query_len;

} InfoMapping;

typedef struct InfoL1 {
  int64_t l1_l = 0;
  float l1_k = 1.0f;
  int64_t l1_lmin = 0;
  int64_t l1_lmax = 0;
  float l1_confidence_abs = 0;
  float l1_std = 0;
  int64_t l1_rough_start = 0;
  int64_t l1_rough_end = 0;
} InfoL1;

typedef struct InfoAlignment {
//  bool is_mapped = false;
  bool is_aligned = false;
  bool is_reverse = false;
  int64_t pos_start = 0;
  int64_t pos_end = 0;
  std::string cigar = "*";
  int64_t edit_distance = 0;
  int64_t alignment_score = 0;
  int64_t mapping_quality = 0;
  double evalue = 0.0f;
  int64_t num_same_mappings = 0;      // How many mapping positions have exactly the same score.

  int64_t num_eq_ops = 0;
  int64_t num_x_ops = 0;
  int64_t num_i_ops = 0;
  int64_t num_d_ops = 0;

  std::string unmapped_reason = "Not processed.";

} InfoAlignment;

class PathGraphEntry {
 public:
  PathGraphEntry(const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const Region &region, InfoMapping *mapping_data=NULL, InfoL1 *l1_data=NULL, InfoAlignment *alignment_data=NULL);
  ~PathGraphEntry();

  void Set(const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const Region &region, InfoMapping *mapping_data=NULL, InfoL1 *l1_data=NULL, InfoAlignment *alignment_data=NULL);
  void AddSecondaryAlignmentData(InfoAlignment alignment_info);

  std::string GenerateSAM(bool is_primary, int64_t verbose_sam_output) const;
  std::string GenerateAMOS() const;

  float CalcDistanceRatio() const;
  float CalcDistanceRatioSuppress() const;
  float CalcFPFilter() const;
  static float CalcFPFilterStatic(const InfoMapping &mapping_info, const Index *index, const SingleSequence *read, const ProgramParameters *parameters);

  void Verbose(FILE *fp) const;

  std::string VerboseToString(std::string delimiter="\r") const;
  float get_fpfilter() const;
  void set_fpfilter(float fpfilter);
  const InfoL1& get_l1_data() const;
  void set_l1_data(InfoL1& l1Data);
  const InfoMapping& get_mapping_data() const;
  void set_mapping_data(InfoMapping& mappingData);
  const Region& get_region_data() const;
  void set_region_data(Region& regionData);
  float get_fpfilter_cov_bases();
  void set_fpfilter_cov_bases(float fpfilterCovBases);
  float get_fpfilter_query_len();
  void set_fpfilter_query_len(float fpfilterQueryLen);
  float get_fpfilter_read_len();
  void set_fpfilter_read_len(float fpfilterReadLen);
  float get_fpfilter_std();
  void set_fpfilter_std(float fpfilterStd);
  InfoAlignment& get_alignment_primary();
  void set_alignment_primary(InfoAlignment& alignmentPrimary);
  std::vector<InfoAlignment>& get_alignments_secondary();
  void set_alignments_secondary(std::vector<InfoAlignment>& alignmentsSecondary);

//  inline bool IsLessThan(const PathGraphEntry &op1);
//
//  int64_t lcs_length;
//  int64_t covered_bases;
//  int64_t region_kmer_count;
//  Range query;
//  Range reference;
//  Range distance;
//  float ratio;
//  float ratio_suppress;
////  int64_t local_scores_id;
////  int64_t region_id;
////  int64_t source_node;
////  int64_t sink_node;
//  Range actual_query;
//  Range actual_reference;
//  float deviation;
//  Region region;
//  int64_t query_id;
//  int64_t query_length;
//  std::string qname;
//  std::string rname;

//  int64_t num_covering_kmers;



 private:
  const SingleSequence *read_;
  const Index *index_;
  const ProgramParameters *parameters_;

  InfoMapping mapping_info_;
  // Parameters holding the L1 results
  InfoL1 l1_info_;
//  std::vector<InfoAlignment> alignments_;
  InfoAlignment alignment_primary_;
  std::vector<InfoAlignment> alignments_secondary_;

  Region region_info_;

  float fpfilter_;
  float fpfilter_cov_bases_;
  float fpfilter_query_len_;
  float fpfilter_std_;
  float fpfilter_read_len_;

  std::string GenerateSAMFromInfoAlignment_(const InfoAlignment &alignment_info, bool is_primary, int64_t verbose_sam_output) const;



//  int64_t l1_l;
//  float l1_k;
//  int64_t l1_lmin;
//  int64_t l1_lmax;
//  int64_t l1_confidence_abs;
//  int64_t l1_std;
//  int64_t l1_rough_start;
//  int64_t l1_rough_end;
//
//  int64_t edit_distance;
//  int64_t edit_distance_end_pos;
//  float fpfilter;
//  float fpfilter_score_std;
//  float fpfilter_score_cov_bases;
//  float fpfilter_score_read_len;
//  float fpfilter_score_query_len;
//
//  int64_t num_eq_ops;
//  int64_t num_x_ops;
//  int64_t num_i_ops;
//  int64_t num_d_ops;
//
//  SamEntry sam;



//  OverlapEntry ovl;

// private:
  static float CalcFPFactorReadLength_(int64_t read_length, int64_t db_size);
  // Parameter query_end must include the k_graph length.
  static float CalcFPFactorQueryLength_(int64_t query_start, int64_t query_end, int64_t read_length);
  // Parameter query_end must include the k_graph length.
  static float CalcFPFactorCovBases_(int64_t covered_bases, int64_t query_start, int64_t query_end, float error_rate);
  static float CalcFPFactorStd_(float std, int64_t read_length, float error_rate);
};

// Used until 27.01.2015.
//struct less_than_key
//{
//    inline bool operator() (const PathGraphEntry& op1, const PathGraphEntry& op2) {
////      if (op1.lcs_length < op2.lcs_length)
////        return true;
////      else if (op1.lcs_length > op2.lcs_length)
////        return false;
////      else {
////        if (op1.covered_bases != op2.covered_bases)
////          return (op1.covered_bases < op2.covered_bases);
////        else
////          return (op1.region_kmer_count < op2.region_kmer_count);
////      }
//
//      if (op1.covered_bases < op2.covered_bases)
//        return true;
//      else if (op1.covered_bases > op2.covered_bases)
//        return false;
//      else {
//        if (op1.lcs_length != op2.lcs_length)
//          return (op1.lcs_length < op2.lcs_length);
//        else
//          return (op1.region_kmer_count < op2.region_kmer_count);
//      }
//
////      return (op1.depth < op2.depth);
////        return (op1.ratio_suppress > op2.ratio_suppress);
//      return false;
//    }
//};

///// Used until 11.03.2015.
//struct path_graph_entry_greater_than_key
//{
//    inline bool operator() (const PathGraphEntry* op1, const PathGraphEntry* op2) {
//      if (op1->num_covering_kmers > op2->num_covering_kmers)
//        return true;
//      else if (op1->num_covering_kmers < op2->num_covering_kmers)
//        return false;
//      else {
//        if (op1->lcs_length != op2->lcs_length) {
//          return (op1->lcs_length > op2->lcs_length);
//        } else {
//            return (op1->fpfilter > op2->fpfilter);
//        }
//      }
//
//      return false;
//    }
//};

struct path_graph_entry_greater_than_key
{
    inline bool operator() (PathGraphEntry* op1, PathGraphEntry* op2) {
      if (op1->get_fpfilter() != op2->get_fpfilter())
        return op1->get_fpfilter() > op2->get_fpfilter();
      else {
        if (op1->get_mapping_data().num_covering_kmers != op2->get_mapping_data().num_covering_kmers) {
          return (op1->get_mapping_data().num_covering_kmers > op2->get_mapping_data().num_covering_kmers);
        } else {
            return (op1->get_mapping_data().lcs_length > op2->get_mapping_data().lcs_length);
        }
      }

      return false;
    }
};

#endif /* PATH_GRAPH_ENTRY_H_ */



//The overlaps.afg file is a text file that contains the overlap information in an XML-like format that looks like this:
//
//{OVL
//adj:N
//rds:1,2
//scr:0
//ahg:139
//bhg:1685
//}
//
//This means that reads 1 and 2 overlap in "normal" orientation, meaning the end of read 1 overlaps the beginning of read 2. The read ids are internal numeric ids assigned by AMOS, the file reads.bnk/RED.0.map has the mapping from read number to the read name in the original fasta file. The "ahg" (ahang) is the number of bases of the first read in the overlap (read 1) that occur before the overlap and the "bhg" (bhang) are the number of bases of the second read (read 2) that extend past the overlap:
//
//1:  =============================================>        <-- 1685bp -->
//2:   <-- 139bp -->  =========================================================>
//
//Because the reads may be from either strand of the DNA, the other valid orientation is "I" for innie:
//
//1:  =============================================>        <-- 1685bp -->
//2:   <-- 139bp -->  <=========================================================
//
//Note the ahand and/or bhang may take negative values if the second read overlaps the beginning of the first read (negative ahang), or the first read extends past the second read (negative bhang)
//
//A:   <-- -ahang -->  =========================================================>
//B:  =============================================>        <-- -bhang -->



//Output format
//: for each pair of overlapping sequences you will have to report their relative
//orientation (assume first sequence is always “forward”) as well as the amount each sequence “hangs”
//(extends) outside of the overlap region.
//AMOS format for overlaps:
//{OVL
//adj:N
//rds:0,279
//scr:0
//ahg:-32
//bhg:-32
//}
//Details:
//rds – The overlapping reads, referenced by their IIDs.
//adj - Read adjacency. [N]ormal, [A]nti-normal, [I]nnie, [O]utie which are EB, BE, EE, BB overlaps
//respectively. W
//here I(nnie) indicates the second sequence is in the opposite orientation than the first
//one. N(ormal) would indicate both sequences are in the same orientation.
//ahg - Ahang. Length of the non-overlapping portion of the first read.
//bhg - Bhang. Length of the non-overlapping portion of the second read.
//Note
//: Positive Ahang values indicate seqA starts before seqB, negative indicate seqB starts before
//seqA.
//scr – An unsigned integer overlap score.
//flg – Universal_t flags plus one additional flag (A/B/C), default to zero if unspecified.
