/*
 * anchor_aligner.cc
 *
 *  Created on: Aug 23, 2017
 *      Author: isovic
 */

#include "anchor_aligner.h"
#include "aligner_util.hpp"

namespace is {

std::shared_ptr<AnchorAligner> createAnchorAligner(std::shared_ptr<is::AlignerBase> aligner) {
  return std::shared_ptr<AnchorAligner>(new AnchorAligner(aligner));
}

AnchorAligner::AnchorAligner(std::shared_ptr<is::AlignerBase> aligner) : aligner_(aligner) {

}

AnchorAligner::~AnchorAligner() {

}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalEndToEnd(const char *ref, const char *query, const std::vector<AlignmentAnchor>& anchors) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);
  std::vector<AlignmentAnchor> final_anchors;
  if (anchors.size() > 0) {
    final_anchors.emplace_back(AlignmentAnchor(anchors.front().qstart, anchors.back().qend, anchors.front().rstart, anchors.back().rend));
  }
  return GlobalAnchored(ref, query, final_anchors);
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchored(const char *ref, const char *query, const std::vector<AlignmentAnchor>& anchors) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);

  if (anchors.size() == 0) {
    return result;
  }

  result->cigar.clear();

  // Align between anchors.
  for (int64_t i = 0; i < (anchors.size() - 1); i++) {
    aligner_->Global(query + anchors[i].qstart, anchors[i+1].qend - anchors[i].qstart,
                      ref + anchors[i].rstart, anchors[i+1].rend - anchors[i].rstart);

    auto aln_result = aligner_->getResults();

    auto left_part = ExtractCigarBetweenQueryCoords(aln_result->cigar,
                                                    0,
                                                    anchors[i+1].qstart - anchors[i].qstart); // Leave next anchor for the next alignment.

    result->cigar.insert(result->cigar.end(), left_part.begin(), left_part.end());
  }

  // Align the last anchor.
  aligner_->Global(query + anchors.back().qstart, anchors.back().qend - anchors.back().qstart,
                    ref + anchors.back().rstart, anchors.back().rend - anchors.back().rstart);

  auto aln_result = aligner_->getResults();

  result->cigar.insert(result->cigar.end(), aln_result->cigar.begin(), aln_result->cigar.end());

  // Fill the other alignment info.
  result->score = -1;
  result->edit_dist = EditDistFromExtCIGAR(result->cigar);

  result->position = is::AlignmentPosition(anchors.front().qstart, anchors.back().qend, anchors.front().rstart, anchors.back().rend);
  result->k = -1;
  result->rv = is::AlignmentReturnValue::OK;

  return result;
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchoredWithExtend(const char *ref, const char *query, const std::vector<AlignmentAnchor>& anchors) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);
  return result;
}

}
