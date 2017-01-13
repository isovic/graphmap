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

  auto& trans_id_to_regions = index->get_trans_id_to_regions();
  if (trans_id_to_regions.size() == 0) { return 2; }

  std::vector<uint8_t> alignment;
  alignment.reserve(aln->alignment.size());

  auto regions_it = trans_id_to_regions.find(aln->ref_header);
  if (regions_it == trans_id_to_regions.end()) { return 3; }
  auto& regions = regions_it->second;
  if (regions.size() == 0) { return 4; }

  int64_t pos_on_read = aln->query_start;
  int64_t pos_on_ref = aln->ref_start;

  int64_t reg_id = 0;
  int64_t region_genome_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
  int64_t region_genome_end = regions[reg_id].second - 1;      // They should be inclusive as well.
  int64_t region_genome_len = region_genome_end - region_genome_start + 1;      // +1 because it's inclusive.

//  int64_t ref_start_on_genome = region_start + aln->ref_start;
//  int64_t ref_end_on_genome = This is not as trivial

  int64_t processed_len = 0;

  for (reg_id=0; reg_id<regions.size() && pos_on_ref < processed_len; reg_id++) {
    processed_len += region_genome_len;

    region_genome_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
    region_genome_end = regions[reg_id].second - 1;      // They should be inclusive as well.
    region_genome_len = region_genome_end - region_genome_start;      // +1 because it's inclusive.
  }

  // Recalculate the new position of the alignment (on the genome).
  int64_t new_ref_start = (pos_on_ref - processed_len) + region_genome_start;

//  printf ("\nreg_id = %ld, region_start = %ld, region_end = %ld, region_len = %ld, aln->ref_start = %ld, aln->ref_end = %ld, processed_len = %ld\n",
//          reg_id, region_genome_start, region_genome_end, region_genome_len, aln->ref_start, aln->ref_end, processed_len);
//  fflush(stdout);

//  for (int64_t i=0; i<aln->alignment.size(); i++) {
//    printf ("%d", aln->alignment[i]);
//  }
//  printf ("\n");
//  fflush(stdout);

  for (int64_t i=0; i<aln->alignment.size(); i++) {
//    printf ("[i = %ld/%ld]\n", i, aln->alignment.size());
//    fflush(stdout);

    int8_t op = aln->alignment[i];

//    printf ("Tu sam 3!\n");
//    fflush(stdout);

    // The op needs to be copied regardless of whether there is a gap.
    alignment.push_back(op);

//    printf ("Tu sam 2!\n");
//    fflush(stdout);

    // TODO: Or >=, depends whether the coordinates are 0 or 1 based.
    // Check if there should be a gap following the current base.
    // If so, add a bunch of 'N's.
    // First check if this is the last region. There can be I/S/H operations after that.
    if ((reg_id+1) < (regions.size()) && pos_on_ref >= (processed_len + region_genome_len - 1) && (op == EDLIB_M || op == EDLIB_EQUAL || op == EDLIB_X)) {
//      printf ("  reg_id = %ld, region_start = %ld, region_end = %ld, region_len = %ld, regions.size() = %ld\n",
//              reg_id, region_genome_start, region_genome_end, region_genome_len, regions.size());
//      fflush(stdout);
      int64_t next_region_start = regions[reg_id+1].first - 1;
      alignment.insert(alignment.end(), (next_region_start - region_genome_end - 1), EDLIB_N);
//      printf ("\nAdding Ns!\n");
//      fflush(stdout);

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
      pos_on_ref += 1;
      pos_on_read += 1;
    } else if (op == EDLIB_D) {
      pos_on_ref += 1;
    } else if (op == EDLIB_I || op == EDLIB_S) {
      pos_on_read += 1;
    }
  }

  // Copy the new alignment.
  aln->alignment = alignment;

  std::string tid = aln->ref_header;
  auto it_tid = index->get_trans_id_to_genome_id().find(tid);
  if (it_tid == index->get_trans_id_to_genome_id().end()) { return 1; }
  aln->ref_header = it_tid->second;
  aln->ref_start = new_ref_start;

  return 0;
}
