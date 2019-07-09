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

#include <deque>
#include <stack>

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

  std::shared_ptr<AlignmentResult> CreateAlignmentResult(int qstart, int qend, int rstart, int rend, std::vector<is::CigarOp> rez);

  void AdjustEnds(int left_offset_ref, int right_offset_ref, const char *query, const char *ref, int64_t *start_position_ref, int64_t *start_position_read, int number_of_bases, std::stack<is::CigarOp> *cigar_stack, std::deque<is::CigarOp> *cigar_queue, bool type);

  /* Sorts anchors and then performs global alignment between the minimum and maximum anchor coordinates.
  */
  std::shared_ptr<AlignmentResult> GlobalEndToEnd(int64_t abs_ref_id, std::shared_ptr<is::MinimizerIndex> index, const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors);

  /* Sorts the anchors, and aligns every neighboring pair of anchors. It does not extend beyond
     the ends of the first and last anchor.
  */
  std::shared_ptr<AlignmentResult> GlobalAnchored(int64_t abs_ref_id, std::shared_ptr<is::MinimizerIndex> index, const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors, bool type);
  std::shared_ptr<AlignmentResult> GlobalAnchoredWithClipping(const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors);

  /* Sorts the anchors, and aligns every neighboring pair of anchors. This extends alignments beyond
     the ends of the first and last anchor in an attempt to produce end-to-end alignment.
  */
  std::shared_ptr<AlignmentResult> GlobalAnchoredWithExtend(int64_t abs_ref_id, std::shared_ptr<is::MinimizerIndex> index, const char *query, int64_t qlen, const char *ref, int64_t rlen,
                                                            const std::vector<AlignmentAnchor>& anchors, int32_t bandwidth, int32_t zdrop, bool type);

 private:
  AnchorAligner(const AnchorAligner&) = delete;
  AnchorAligner& operator=(const AnchorAligner&) = delete;

  AnchorAligner(std::shared_ptr<is::AlignerBase> aligner);

  const std::shared_ptr<is::AlignerBase> aligner_;
};

}

#endif
