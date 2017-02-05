/*
 * transcriptome_mod.cc
 *
 *  Created on: Jan 5, 2017
 *      Author: isovic
 */

#include "alignment/transcriptome_mod.h"
#include "alignment/cigargen.h"
#include "index/index_spaced_hash_fast.h"

int ConvertFromTranscriptomeToGenomeAln(const IndexSpacedHashFast *index, AlignmentResults *aln) {
  if (index->is_transcriptome() == false) { return 1; }

  // Get the info about the source chromosome.
  std::string tid = aln->ref_header;
  auto it_tid = index->get_trans_id_to_genome_id().find(tid);
  if (it_tid == index->get_trans_id_to_genome_id().end()) { return 1; }

  auto& chr_info = it_tid->second;
  const std::string& chr_name = chr_info.first;  // The header of the chromosome the tid corresponds to.
  char chr_orient = chr_info.second;    // Orientation of the transcriptome on the genome.

  // Get the split regions.
  auto& trans_id_to_regions = index->get_trans_id_to_regions();
  if (trans_id_to_regions.size() == 0) { return 2; }

  auto regions_it = trans_id_to_regions.find(aln->ref_header);
  if (regions_it == trans_id_to_regions.end()) { return 3; }
  auto& regions = regions_it->second;
  if (regions.size() == 0) { return 4; }

  if (chr_orient == '+') {
    std::vector<uint8_t> alignment;
    alignment.reserve(aln->alignment.size());

    int64_t pos_on_read = aln->query_start;
    int64_t pos_on_ref = aln->ref_start;

    int64_t reg_id = 0;
    int64_t region_genome_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
    int64_t region_genome_end = regions[reg_id].second - 1;      // They should be inclusive as well.
    int64_t region_genome_len = region_genome_end - region_genome_start + 1;      // +1 because it's inclusive.

    // This will keep the currently processed part of the transcriptome (start coordinate).
    int64_t processed_len = 0;

    // Find the exon where the alignment actually begins.
    for (reg_id=0; reg_id<regions.size() && pos_on_ref < processed_len; reg_id++) {
      processed_len += region_genome_len;

      region_genome_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
      region_genome_end = regions[reg_id].second - 1;      // They should be inclusive as well.
      region_genome_len = region_genome_end - region_genome_start;      // +1 because it's inclusive.
    }

    // Recalculate the new position of the alignment (on the genome).
    int64_t new_ref_start = (pos_on_ref - processed_len) + region_genome_start;
    int64_t new_ref_end = new_ref_start;

    // Process each alignment operation individually.
    for (int64_t i=0; i<aln->alignment.size(); i++) {
      int8_t op = aln->alignment[i];

      // The op needs to be copied regardless of whether there is a gap.
      alignment.push_back(op);

      // Check if there should be a gap following the current base. If so, add a bunch of 'N's.
      // First check if this is the last region. There can be I/S/H operations after that.
      if ((reg_id+1) < (regions.size()) && pos_on_ref >= (processed_len + region_genome_len - 1) && (op == EDLIB_M || op == EDLIB_EQUAL || op == EDLIB_X)) {
        int64_t next_region_start = regions[reg_id+1].first - 1;
        alignment.insert(alignment.end(), (next_region_start - region_genome_end - 1), EDLIB_N);

        processed_len += region_genome_len;
        reg_id += 1;
        if (reg_id > regions.size()) { break; }
        region_genome_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
        region_genome_end = regions[reg_id].second - 1;
        region_genome_len = region_genome_end - region_genome_start + 1;      // +1 because it's inclusive.
      }

      if (op == EDLIB_M || op == EDLIB_EQUAL || op == EDLIB_X) {
        // Matches, mismatches and deletions modify the position on the reference.
        new_ref_end = region_genome_start + (pos_on_ref - processed_len);
        pos_on_ref += 1;
        pos_on_read += 1;
      } else if (op == EDLIB_D) {
        // Matches, mismatches and deletions modify the position on the reference.
        new_ref_end = region_genome_start + (pos_on_ref - processed_len);
        pos_on_ref += 1;
      } else if (op == EDLIB_I || op == EDLIB_S) {
        pos_on_read += 1;
      }
    }

    // Modify the rname.
    aln->ref_header = chr_name; // Genome chromosome ID.

    // Copy the new alignment.
    aln->alignment = alignment;
    aln->ref_start = new_ref_start;
    aln->ref_end = new_ref_end;

  } else if (chr_orient == '-') {
    std::vector<uint8_t> alignment;
    alignment.reserve(aln->alignment.size());

    int64_t pos_on_ref = aln->ref_len - aln->ref_end - 1;   // Reverse the start position.

    // This will keep the currently processed part of the transcriptome (start coordinate).
    int64_t processed_len = 0;
    int64_t reg_id = 0;
    int64_t region_genome_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
    int64_t region_genome_end = regions[reg_id].second - 1;      // They should be inclusive as well.
    int64_t region_genome_len = region_genome_end - region_genome_start + 1;      // +1 because it's inclusive.

//    printf ("\n");

    // Find the exon where the alignment actually begins.
    for (reg_id=0; reg_id<regions.size() && pos_on_ref < processed_len && pos_on_ref ; reg_id++) {
//      printf ("(loop) reg_id = %ld, region_start = %ld, region_end = %ld, region_len = %ld, aln->ref_start = %ld, aln->ref_end = %ld, aln->ref_len = %ld, pos_on_ref = %ld, processed_len = %ld\n",
//              reg_id, region_genome_start, region_genome_end, region_genome_len, aln->ref_start, aln->ref_end, aln->ref_len, pos_on_ref, processed_len);
//      fflush(stdout);

      processed_len += region_genome_len;

      region_genome_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
      region_genome_end = regions[reg_id].second - 1;      // They should be inclusive as well.
      region_genome_len = region_genome_end - region_genome_start;      // +1 because it's inclusive.
    }

    // Recalculate the new position of the alignment (on the genome).
    int64_t new_ref_start = (pos_on_ref - processed_len) + region_genome_start;
    int64_t new_ref_end = new_ref_start;

//    printf ("(start) reg_id = %ld, region_start = %ld, region_end = %ld, region_len = %ld, aln->ref_start = %ld, aln->ref_end = %ld, aln->ref_len = %ld, pos_on_ref = %ld, processed_len = %ld\n",
//            reg_id, region_genome_start, region_genome_end, region_genome_len, aln->ref_start, aln->ref_end, aln->ref_len, pos_on_ref, processed_len);
//    fflush(stdout);

  //  for (int64_t i=0; i<aln->alignment.size(); i++) {
  //    printf ("%d", aln->alignment[i]);
  //  }
  //  printf ("\n");
  //  fflush(stdout);

    // Process each alignment operation individually.
    for (int64_t i=(aln->alignment.size()-1); i>=0; i--) {
      int8_t op = aln->alignment[i];

      // The op needs to be copied regardless of whether there is a gap.
      alignment.push_back(op);

      // Check if there should be a gap following the current base. If so, add a bunch of 'N's.
      // First check if this is the last region. There can be I/S/H operations after that.
      if ((reg_id+1) < (regions.size()) && pos_on_ref >= (processed_len + region_genome_len - 1) && (op == EDLIB_M || op == EDLIB_EQUAL || op == EDLIB_X)) {
  //      printf ("  reg_id = %ld, region_start = %ld, region_end = %ld, region_len = %ld, regions.size() = %ld\n",
  //              reg_id, region_genome_start, region_genome_end, region_genome_len, regions.size());
  //      fflush(stdout);
        int64_t next_region_start = regions[reg_id+1].first - 1;
        alignment.insert(alignment.end(), (next_region_start - region_genome_end - 1), EDLIB_N);

        processed_len += region_genome_len;
        reg_id += 1;
        if (reg_id > regions.size()) { break; }
        region_genome_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
        region_genome_end = regions[reg_id].second - 1;
        region_genome_len = region_genome_end - region_genome_start + 1;      // +1 because it's inclusive.
  //      printf ("  reg_id = %ld, region_start = %ld, region_end = %ld, region_len = %ld, aln->ref_start = %ld, aln->ref_end = %ld, processed_len = %ld\n",
  //              reg_id, region_genome_start, region_genome_end, region_genome_len, aln->ref_start, aln->ref_end, processed_len);
  //      fflush(stdout);
      }

      if (op == EDLIB_M || op == EDLIB_EQUAL || op == EDLIB_X) {
        // Matches, mismatches and deletions modify the position on the reference.
        new_ref_end = region_genome_start + (pos_on_ref - processed_len);
        // Position on read should change as well, but not important here.
        pos_on_ref += 1;
      } else if (op == EDLIB_D) {
        // Matches, mismatches and deletions modify the position on the reference.
        new_ref_end = region_genome_start + (pos_on_ref - processed_len);
        pos_on_ref += 1;
      } else if (op == EDLIB_I || op == EDLIB_S) {
        // Only position on read should change. We don't need this.
      }
    }

    // Modify the rname.
    aln->ref_header = chr_name; // Genome chromosome ID.

    // Copy the new alignment.
    aln->alignment = alignment;
    aln->ref_start = new_ref_start;
    aln->ref_end = new_ref_end;
    aln->is_reverse = !aln->is_reverse;
    aln->orientation = (aln->orientation == kForward) ? (kReverse) : (kForward);

//    printf ("new_ref_start = %ld\n", new_ref_start);
//    printf ("new_ref_end = %ld\n", new_ref_end);
//    fflush(stdout);

  } else {
    return 3;

  }

  return 0;
}
