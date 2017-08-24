/*
 * anchor_aligner.cc
 *
 *  Created on: Aug 23, 2017
 *      Author: isovic
 */

#include "anchor_aligner.h"

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

  printf ("anchors.front().qstart = %ld\n", anchors.front().qstart);
  printf ("anchors.front().qend = %ld\n", anchors.front().qend);
  printf ("anchors.front().rstart = %ld\n", anchors.front().rstart);
  printf ("anchors.front().rend = %ld\n\n", anchors.front().rend);

  printf ("anchors.back().qstart = %ld\n", anchors.back().qstart);
  printf ("anchors.back().qend = %ld\n", anchors.back().qend);
  printf ("anchors.back().rstart = %ld\n", anchors.back().rstart);
  printf ("anchors.back().rend = %ld\n\n", anchors.back().rend);

  aligner_->Global(query + anchors.front().qstart, anchors.back().qend - anchors.front().qstart,
                   ref + anchors.front().rstart, anchors.back().rend - anchors.front().rstart);

  return result;
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchored(const char *ref, const char *query, const std::vector<AlignmentAnchor>& anchors) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);
  assert(anchors.size() > 1);

  if (anchors.size() > 1) {
    for (int64_t i = 0; i < (anchors.size() - 1); i++) {
      printf ("anchors[i].qstart = %ld\n", anchors[i].qstart);
      printf ("anchors[i].qend = %ld\n", anchors[i].qend);
      printf ("anchors[i].rstart = %ld\n", anchors[i].rstart);
      printf ("anchors[i].rend = %ld\n\n", anchors[i].rend);

      printf ("anchors[i+1].qstart = %ld\n", anchors[i+1].qstart);
      printf ("anchors[i+1].qend = %ld\n", anchors[i+1].qend);
      printf ("anchors[i+1].rstart = %ld\n", anchors[i+1].rstart);
      printf ("anchors[i+1].rend = %ld\n\n", anchors[i+1].rend);

      aligner_->Global(query + anchors[i].qstart, anchors[i+1].qend - anchors[i].qstart,
                       ref + anchors[i].rstart, anchors[i+1].rend - anchors[i].rstart);

    }
  }

  return result;
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchoredWithExtend(const char *ref, const char *query, const std::vector<AlignmentAnchor>& anchors) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);
  return result;
}

}
