/*
 * anchor_aligner.h
 *
 *  Created on: Aug 23, 2017
 *      Author: isovic
 */

#ifndef SRC_ANCHOR_ALIGNER_H_
#define SRC_ANCHOR_ALIGNER_H_

#include <memory>
#include "aligner_base.h"
#include "containers/results.h"

namespace is {

class AnchorAligner;

std::shared_ptr<AnchorAligner> createAnchorAligner(std::shared_ptr<is::AlignerBase> aligner);

class AlignmentAnchor {
 public:
  AlignmentAnchor() : qstart(0), qend(0), rstart(0), rend(0) { }
  AlignmentAnchor(int64_t _qstart, int64_t _qend,
                  int64_t _rstart, int64_t _rend) :
                    qstart(_qstart), qend(_qend), rstart(_rstart), rend(_rend) { }

  int64_t qstart, qend;
  int64_t rstart, rend;
};

class AnchorAligner {
 public:
  friend std::shared_ptr<AnchorAligner> createAnchorAligner(std::shared_ptr<is::AlignerBase> aligner);

  ~AnchorAligner();

  /* Sorts anchors and then performs global alignment between the minimum and maximum anchor coordinates.
  */
  std::shared_ptr<AlignmentResult> GlobalEndToEnd(const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors);

  /* Sorts the anchors, and aligns every neighboring pair of anchors. It does not extend beyond
     the ends of the first and last anchor.
  */
  std::shared_ptr<AlignmentResult> GlobalAnchored(const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors);

  /* Sorts the anchors, and aligns every neighboring pair of anchors. This extends alignments beyond
     the ends of the first and last anchor in an attempt to produce end-to-end alignment.
  */
  std::shared_ptr<AlignmentResult> GlobalAnchoredWithExtend(const char *query, int64_t qlen, const char *ref, int64_t rlen,
                                                            const std::vector<AlignmentAnchor>& anchors, int32_t bandwidth, int32_t zdrop);

 private:
  AnchorAligner(const AnchorAligner&) = delete;
  AnchorAligner& operator=(const AnchorAligner&) = delete;

  AnchorAligner(std::shared_ptr<is::AlignerBase> aligner);

  const std::shared_ptr<is::AlignerBase> aligner_;
};

}

#endif
