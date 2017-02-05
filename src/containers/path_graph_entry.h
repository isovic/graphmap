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
#include "program_parameters.h"
//#include "sam/sam_entry.h"
#include "containers/vertices.h"
#include "containers/results.h"

// N(ormal) would indicate both sequences are in the same orientation.
#define OVERLAP_NORMAL          "N"
// I(nnie) indicates the second sequence is in the opposite orientation than the first
#define OVERLAP_INNIE           "I"
#define OVERLAP_OUTIE           "O"
#define OVERLAP_ANTI_NORMAL     "A"




class PathGraphEntry {
 public:
  PathGraphEntry(const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const Region &region, MappingResults *mapping_data=NULL, L1Results *l1_data=NULL, AlignmentResults *alignment_data=NULL);
  ~PathGraphEntry();

  void Set(const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const Region &region, MappingResults *mapping_data=NULL, L1Results *l1_data=NULL, AlignmentResults *alignment_data=NULL);
//  void AddSecondaryAlignmentData(AlignmentResults alignment_info);

  std::string GenerateSAM(bool is_primary, int64_t verbose_sam_output) const;
  std::string GenerateAFG() const;
  std::string GenerateM5(bool is_primary, int64_t verbose_sam_output) const;
  std::string GenerateMHAP(bool is_primary, int64_t verbose_sam_output) const;

  float CalcDistanceRatio() const;
  float CalcDistanceRatioSuppress() const;
  float CalcFPFilter() const;
  static float CalcFPFilterStatic(const MappingResults &mapping_info, const Index *index, const SingleSequence *read, const ProgramParameters *parameters);

  void Verbose(FILE *fp) const;

  std::string VerboseToString(std::string delimiter="\r") const;
  std::string VerboseInfoToString(std::string delimiter="\t") const;

  void SetMapped(bool is_mapped);
  bool IsMapped();
  void SetAligned(bool is_aligned);
  bool IsAligned();

  void AddAlignmentData(AlignmentResults alignment_info);

  float get_fpfilter() const;
  void set_fpfilter(float fpfilter);
  const L1Results& get_l1_data() const;
  void set_l1_data(L1Results& l1Data);
  MappingResults& mapping_data();
   const MappingResults& get_mapping_data() const;
  void set_mapping_data(MappingResults& mappingData);
  Region& region_data();
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
  std::vector<AlignmentResults>& get_alignments();
  void set_alignments(std::vector<AlignmentResults>& alignments);
//  AlignmentResults& get_alignment_primary();
//  void set_alignment_primary(AlignmentResults& alignmentPrimary);
//  std::vector<AlignmentResults>& get_alignments_secondary();
//  void set_alignments_secondary(std::vector<AlignmentResults>& alignmentsSecondary);
  MappingMetadata& get_mapping_metadata();
  void set_mapping_metadata(MappingMetadata& mappingMetadata);

 private:
  const SingleSequence *read_;
  const Index *index_;
  const ProgramParameters *parameters_;

  MappingResults mapping_info_;
  // Parameters holding the L1 results
  L1Results l1_info_;
//  std::vector<InfoAlignment> alignments_;
  std::vector<AlignmentResults> alignments_;
//  AlignmentResults alignment_primary_;
//  std::vector<AlignmentResults> alignments_secondary_;
  MappingMetadata mapping_metadata_;

  Region region_info_;

  float fpfilter_;
  float fpfilter_cov_bases_;
  float fpfilter_query_len_;
  float fpfilter_std_;
  float fpfilter_read_len_;

  std::string GenerateSAMFromInfoAlignment_(const AlignmentResults &alignment_info, const MappingMetadata &mapping_metadata, bool is_primary, int64_t verbose_sam_output) const;
  std::string GenerateAFGFromInfoAlignment_(const AlignmentResults &alignment_info) const;
  std::string GenerateAFGFromInfoMappping(const MappingResults &mapping_info) const;
  std::string GenerateM5FromInfoAlignment_(const AlignmentResults &alignment_info, const MappingMetadata &mapping_metadata, bool is_primary, int64_t verbose_sam_output) const;
  std::string GenerateMHAPFromInfoAlignment_(const AlignmentResults &alignment_info, const MappingMetadata &mapping_metadata, bool is_primary, int64_t verbose_sam_output) const;

// private:
  static float CalcFPFactorReadLength_(int64_t read_length, int64_t db_size);
  // Parameter query_end must include the k_graph length.
  static float CalcFPFactorQueryLength_(int64_t query_start, int64_t query_end, int64_t read_length);
  // Parameter query_end must include the k_graph length.
  static float CalcFPFactorCovBases_(int64_t covered_bases, int64_t query_start, int64_t query_end, float error_rate);
  static float CalcFPFactorStd_(float std, int64_t read_length, float error_rate);
};

// This one produced quite good results up until 01.04.2016. Testing with different sorting conditions.
//struct path_graph_entry_greater_than_key
//{
//    inline bool operator() (PathGraphEntry* op1, PathGraphEntry* op2) {
//      if (op1->get_fpfilter() != op2->get_fpfilter())
//        return op1->get_fpfilter() > op2->get_fpfilter();
//      else {
//        if (op1->get_mapping_data().num_covering_kmers != op2->get_mapping_data().num_covering_kmers) {
//          return (op1->get_mapping_data().num_covering_kmers > op2->get_mapping_data().num_covering_kmers);
//        } else {
//            return (op1->get_mapping_data().lcs_length > op2->get_mapping_data().lcs_length);
//        }
//      }
//
//      return false;
//    }
//};

//struct path_graph_entry_greater_than_key
//{
//    inline bool operator() (PathGraphEntry* op1, PathGraphEntry* op2) {
//      if (op1->get_mapping_data().cov_bases_max != op2->get_mapping_data().cov_bases_max)
//        return op1->get_mapping_data().cov_bases_max > op2->get_mapping_data().cov_bases_max;
//      else {
////        if (op1->get_fpfilter() != op2->get_fpfilter()) {
//          return (op1->get_fpfilter() > op2->get_fpfilter());
////        } else {
////          return (op1->get_mapping_data().query_coords.start);
////        }
//
////        if (op1->get_mapping_data().num_covering_kmers != op2->get_mapping_data().num_covering_kmers) {
////          return (op1->get_mapping_data().num_covering_kmers > op2->get_mapping_data().num_covering_kmers);
////        } else {
////        }
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
        return op1->get_mapping_data().cov_bases_max > op2->get_mapping_data().cov_bases_max;
      }

      return false;
    }
};

// Version which works (maybe) slightly worse:
//struct path_graph_entry_greater_than_key
//{
//    inline bool operator() (PathGraphEntry* op1, PathGraphEntry* op2) {
//      if (op1->get_fpfilter() != op2->get_fpfilter())
//        return op1->get_fpfilter() > op2->get_fpfilter();
//      else {
//        if (op1->get_mapping_data().cov_bases_max != op2->get_mapping_data().cov_bases_max) {
//          return op1->get_mapping_data().cov_bases_max > op2->get_mapping_data().cov_bases_max;
//        } else {
//          return (op1->get_mapping_data().num_covering_kmers > op2->get_mapping_data().num_covering_kmers);
//        }
//      }
//
//      return false;
//    }
//};

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
