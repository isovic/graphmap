/*
 * index_util.cc
 *
 *  Created on: Feb 6, 2017
 *      Author: isovic
 */

#include "index_util.h"

namespace is {

std::string GenerateSAMHeader(std::shared_ptr<is::MinimizerIndex> index,
                                  ProgramParameters& parameters) {
  // Output reference sequence information.
  std::stringstream ss_header;

  ss_header << "@HD\t" <<
               "VN:1.0\t" <<
               "SO:unknown" <<
               "\n";

  for (int64_t rid=0; rid<index->get_num_sequences_forward(); rid++) {
    std::string reference_header = index->get_headers()[rid];
    uint64_t rlen = (uint64_t) index->get_reference_lengths()[rid];

    ss_header << "@SQ\t" <<
                "SN:" << reference_header << "\t" <<
                "LN:" << rlen << "" <<
                "\n";
  }

  // If verbose_sam_output == 1, then print out a special version of the PG line. This was used for the web server
  // to omit paths from the output (not to share server sensitive information with users).
  if (parameters.verbose_sam_output == 1) {
    ss_header << "@PG\tID:graphmap\tPN:graphmap";
  } else {
    // Output the command line used to run the process to the file.
    ss_header << "@PG\t" <<
                 "ID:graphmap\t" <<
                 "PN:graphmap\t" <<
                 "CL:" << parameters.command_line << "\t" <<
                 "VN:" << std::string(GRAPHMAP_CURRENT_VERSION) << " compiled on " << std::string(GRAPHMAP_CURRENT_VERSION_RELEASE_DATE);
  }

  return ss_header.str();
}

std::string GenerateSAMHeader(std::shared_ptr<is::Transcriptome> transcriptome) {
  return transcriptome->GenerateSAMHeaders();
}

} /* namespace is */
