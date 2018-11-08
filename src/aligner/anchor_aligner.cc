/*
 * anchor_aligner.cc
 *
 *  Created on: Aug 23, 2017
 *      Author: isovic
 */

#include "anchor_aligner.h"
#include "aligner_util.hpp"

#include <iostream>

namespace is {

std::shared_ptr<AnchorAligner> createAnchorAligner(std::shared_ptr<is::AlignerBase> aligner) {
  return std::shared_ptr<AnchorAligner>(new AnchorAligner(aligner));
}

AnchorAligner::AnchorAligner(std::shared_ptr<is::AlignerBase> aligner) : aligner_(aligner) {

}

AnchorAligner::~AnchorAligner() {

}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalEndToEnd(const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);
  std::vector<AlignmentAnchor> final_anchors;
  if (anchors.size() > 0) {
    final_anchors.emplace_back(AlignmentAnchor(anchors.front().qstart, anchors.back().qend, anchors.front().rstart, anchors.back().rend));
  }
  return GlobalAnchored(query, qlen, ref, rlen, final_anchors, true);
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchoredWithClipping(const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);

  if (anchors.size() == 0) {
    return result;
  }

  result->cigar.clear();
  result->score = 0;

  // Align between anchors.
  for (int64_t i = 0; i < (anchors.size() - 1); i++) {
    aligner_->Global(query + anchors[i].qstart, anchors[i+1].qend - anchors[i].qstart,
                      ref + anchors[i].rstart, anchors[i+1].rend - anchors[i].rstart, false);

    auto aln_result = aligner_->getResults();

    int64_t cigar_length = 0;
    auto left_part = ExtractCigarBetweenQueryCoords(aln_result->cigar,
                                                    0,
                                                    anchors[i+1].qstart - anchors[i].qstart, &cigar_length); // Leave next anchor for the next alignment.

    result->cigar.insert(result->cigar.end(), left_part.begin(), left_part.end());

    // TODO: This is wrong, because it takes the score for the next anchor twice. Need to rework this.
    result->score += aln_result->score;
  }

  // Align the last anchor.
  aligner_->Global(query + anchors.back().qstart, anchors.back().qend - anchors.back().qstart,
                    ref + anchors.back().rstart, anchors.back().rend - anchors.back().rstart, false);

  auto aln_result = aligner_->getResults();

  result->cigar.insert(result->cigar.end(), aln_result->cigar.begin(), aln_result->cigar.end());

  // Add the soft clippings at front and back.
  if ( anchors.front().qstart > 0) {
    result->cigar.insert(result->cigar.begin(), is::CigarOp('S', anchors.front().qstart));
  }
  if ((qlen - anchors.back().qend) > 0) {
    result->cigar.insert(result->cigar.end(), is::CigarOp('S', (qlen - anchors.back().qend)));
  }

  std::vector<is::CigarOp> cigarTmp = result->cigar;
  std::vector<is::CigarOp> cigarTmpContainer;

  std::reverse(cigarTmp.begin(), cigarTmp.end());

  auto result_After_rev = std::shared_ptr<AlignmentResult>(new AlignmentResult);

  std::vector<is::CigarOp> cigarInterContainer;

  bool foundLastGap = true;

  int sub_offset = 0;

  for (auto& c: cigarTmp) {
	  int numberOfBases = (int) c.count;

	  if (foundLastGap) {
		  cigarInterContainer.insert(cigarInterContainer.end(), c);
	  } else {
		  if ((c.op == 'D' || c.op == 'N') && numberOfBases >= 10) {
			  foundLastGap = true;
			  int numberOfMatch = 0;
			  int numberTotal = 0;
			  for (auto& cTmp: cigarTmpContainer) {
				  int numberOfBasesTmp = (int) cTmp.count;
				  if (cTmp.op == '=') {
					  numberOfMatch += numberOfBasesTmp;
				  }
				  if (cTmp.op != 'I') {
					  numberTotal += numberOfBasesTmp;
				  }
			  }
			  if (numberOfMatch > 20000) {
				  cigarInterContainer.insert(cigarInterContainer.end(), cigarTmpContainer.begin(), cigarTmpContainer.end());
			  } else {
				  sub_offset = numberTotal + numberOfBases;
			  }
		  } else {
			  cigarTmpContainer.insert(cigarTmpContainer.end(), c);
		  }
	  }
  }

  foundLastGap = true;
  std::vector<is::CigarOp> cigarTmpInterContainer;

  for (auto& c: cigarInterContainer) {
	  int numberOfBases = (int) c.count;

	  if (foundLastGap) {
		  result_After_rev->cigar.insert(result_After_rev->cigar.end(), c);
	  } else {
		  if ((c.op == 'D' || c.op == 'N') && numberOfBases >= 10) {
			  foundLastGap = true;
			  int numberOfMatch = 0;
			  int numberTotal = 0;
			  for (auto& cTmp: cigarTmpInterContainer) {
				  int numberOfBasesTmp = (int) cTmp.count;
				  if (cTmp.op == '=') {
					  if (numberOfBasesTmp > numberOfMatch) {
						  numberOfMatch = numberOfBasesTmp;
					}
				}
				  if (cTmp.op != 'I') {
					  numberTotal += numberOfBasesTmp;
				}
			}
			if (numberOfMatch >= 5 && numberTotal > 15) {
				result_After_rev->cigar.insert(result_After_rev->cigar.end(), cigarTmpInterContainer.begin(), cigarTmpInterContainer.end());
				result_After_rev->cigar.insert(result_After_rev->cigar.end(), c);
			} else {
				sub_offset += numberTotal + numberOfBases;
			}
		} else {
			cigarTmpInterContainer.insert(cigarTmpInterContainer.end(), c);
		}
	}
  }

  auto result_After = std::shared_ptr<AlignmentResult>(new AlignmentResult);

  std::vector<is::CigarOp> cigarTmp2 = result_After_rev->cigar;
  std::vector<is::CigarOp> cigarTmpContainer2;
  std::reverse(cigarTmp2.begin(), cigarTmp2.end());

  bool foundFirstGap = true;
  int pre_offset = 0;
  std::vector<is::CigarOp> cigarTmpInterContainer22;

  for (auto& c: cigarTmp2) {
	 int numberOfBases = (int) c.count;

	 if (foundFirstGap) {
		 cigarTmpInterContainer22.insert(cigarTmpInterContainer22.end(), c);
	 } else {
		 if ((c.op == 'D' || c.op == 'N') && numberOfBases >= 10) {
			 foundFirstGap = true;
			 int numberOfMatch = 0;
			 int numberTotal = 0;
			 for (auto& cTmp: cigarTmpContainer2) {
				 int numberOfBasesTmp = (int) cTmp.count;
				 if (cTmp.op == '=') {
					 numberOfMatch += numberOfBasesTmp;
				 }
				 if (cTmp.op != 'I') {
					 numberTotal += numberOfBasesTmp;
				 }
			 }
			 if (numberOfMatch > 20000) {
				 cigarTmpInterContainer22.insert(cigarTmpInterContainer22.end(), cigarTmpContainer2.begin(), cigarTmpContainer2.end());
			 } else {
				 pre_offset = numberTotal + numberOfBases;
			 }
		 } else {
			 cigarTmpContainer2.insert(cigarTmpContainer2.end(), c);
		 }
	 }
  }

  int i = 0;
  int wasInsertionAtBegining = false;
  for (auto& c: cigarTmpInterContainer22) {
	  if (i == 0) {
		  if ((c.op == 'D' || c.op == 'N') && c.count >= 10) {
			  int numberOfBases = (int) c.count;
			  pre_offset += numberOfBases;
		  } else if (c.op == 'I') {
			  wasInsertionAtBegining = true;
			  result_After->cigar.insert(result_After->cigar.end(), c);
		  } else {
			  result_After->cigar.insert(result_After->cigar.end(), c);
		  }
	  } else if (i == 1) {
		  if ((c.op == 'D' || c.op == 'N') && c.count >= 10 && wasInsertionAtBegining) {
			  int numberOfBases = (int) c.count;
			  pre_offset += numberOfBases;
		  } else {
			  result_After->cigar.insert(result_After->cigar.end(), c);
		  }
	  } else if (i == cigarTmpInterContainer22.size()-2) {
		  if ((c.op == 'D' || c.op == 'N') && c.count >= 10) {
			  if(cigarTmpInterContainer22[cigarTmpInterContainer22.size()-1].op == 'I') {
				  int numberOfBases = (int) c.count;
				  sub_offset += numberOfBases;
			  } else {
				  result_After->cigar.insert(result_After->cigar.end(), c);
			  }
		  } else {
			  result_After->cigar.insert(result_After->cigar.end(), c);
		  }
	  } else if (i == cigarTmpInterContainer22.size()-1) {
		  if ((c.op == 'D' || c.op == 'N') && c.count >= 10) {
			  int numberOfBases = (int) c.count;
			  sub_offset += numberOfBases;
		  } else {
			  result_After->cigar.insert(result_After->cigar.end(), c);
		  }
	  } else {
		  result_After->cigar.insert(result_After->cigar.end(), c);
	  }
	  i++;
  }

  result_After->position = is::AlignmentPosition(anchors.front().qstart, anchors.back().qend, anchors.front().rstart+(pre_offset-0), anchors.back().rend-(sub_offset-0));
  result_After->k = -1;
  result_After->rv = is::AlignmentReturnValue::OK;

  const int64_t MIN_INTRON_LEN = 10;
  for (auto& c: result_After->cigar) {
    if (c.op == 'D' && c.count >= MIN_INTRON_LEN) {
      c.op = 'N';
    }
  }

  if(result_After->cigar.size() < 2) {
	  return result_After;
  }

  int start_position = (result_After->position.tstart-anchors.front().qstart) + result_After->cigar[0].count;

  bool hasN = false;

  for(int i = 1; i < result_After->cigar.size()-1; i++) {
	int number_of_bases = result_After->cigar[i].count;
	if (result_After->cigar[i].op == 'N') {
		hasN = true;
		int left_offset = 0;
		bool found_left_offset = false;
		if(std::toupper(ref[start_position]) == 'C' && std::toupper(ref[start_position+1]) == 'T') {
			left_offset = 0;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position-1]) == 'C' && std::toupper(ref[start_position]) == 'T') {
			left_offset = -1;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position+1]) == 'C' && std::toupper(ref[start_position+2]) == 'T') {
			left_offset = 1;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position+2]) == 'C' && std::toupper(ref[start_position+3]) == 'T') {
			left_offset = 2;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position-2]) == 'C' && std::toupper(ref[start_position-1]) == 'T') {
			left_offset = -2;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position-3]) == 'C' && std::toupper(ref[start_position-2]) == 'T') {
			left_offset = -3;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position+3]) == 'C' && std::toupper(ref[start_position+4]) == 'T') {
			left_offset = 3;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position+4]) == 'C' && std::toupper(ref[start_position+5]) == 'T') {
			left_offset = 4;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position-4]) == 'C' && std::toupper(ref[start_position-3]) == 'T') {
			left_offset = -4;
			found_left_offset = true;
		}
		int right_offset = 0;
		bool found_right_offset = false;
		if(std::toupper(ref[start_position + number_of_bases - 2]) == 'A' && std::toupper(ref[start_position + number_of_bases - 1]) == 'C') {
			right_offset = 0;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 3]) == 'A' && std::toupper(ref[start_position + number_of_bases - 2]) == 'C') {
			right_offset = -1;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 1]) == 'A' && std::toupper(ref[start_position + number_of_bases]) == 'C') {
			right_offset = 1;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 4]) == 'A' && std::toupper(ref[start_position + number_of_bases - 3]) == 'C') {
			right_offset = -2;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases]) == 'A' && std::toupper(ref[start_position + number_of_bases + 1]) == 'C') {
			right_offset = 2;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 5]) == 'A' && std::toupper(ref[start_position + number_of_bases - 4]) == 'C') {
			right_offset = -3;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases + 1]) == 'A' && std::toupper(ref[start_position + number_of_bases + 2]) == 'C') {
			right_offset = 3;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 6]) == 'A' && std::toupper(ref[start_position + number_of_bases - 5]) == 'C') {
			right_offset = -4;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases + 2]) == 'A' && std::toupper(ref[start_position + number_of_bases + 3]) == 'C') {
			right_offset = 4;
			found_right_offset = true;
		}

		if(found_left_offset && found_right_offset) {
			result_After->cigar[i].count -= left_offset;
			result_After->cigar[i-1].count += left_offset;
			result_After->cigar[i].count += right_offset;
			result_After->cigar[i+1].count -= right_offset;
		} else {
			left_offset = 0;
			found_left_offset = false;
			if(std::toupper(ref[start_position]) == 'G' && std::toupper(ref[start_position+1]) == 'T') {
				left_offset = 0;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position-1]) == 'G' && std::toupper(ref[start_position]) == 'T') {
				left_offset = -1;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position+1]) == 'G' && std::toupper(ref[start_position+2]) == 'T') {
				left_offset = 1;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position+2]) == 'G' && std::toupper(ref[start_position+3]) == 'T') {
				left_offset = 2;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position-2]) == 'G' && std::toupper(ref[start_position-1]) == 'T') {
				left_offset = -2;
				found_left_offset = true;
			}else if(std::toupper(ref[start_position-3]) == 'G' && std::toupper(ref[start_position-2]) == 'T') {
				left_offset = -3;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position+3]) == 'G' && std::toupper(ref[start_position+4]) == 'T') {
				left_offset = 3;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position+4]) == 'G' && std::toupper(ref[start_position+5]) == 'T') {
				left_offset = 4;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position-4]) == 'G' && std::toupper(ref[start_position-3]) == 'T') {
				left_offset = -4;
				found_left_offset = true;
			}

			right_offset = 0;
			found_right_offset = false;
			if(std::toupper(ref[start_position + number_of_bases - 2]) == 'A' && std::toupper(ref[start_position + number_of_bases - 1]) == 'G') {
				right_offset = 0;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 3]) == 'A' && std::toupper(ref[start_position + number_of_bases - 2]) == 'G') {
				right_offset = -1;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 1]) == 'A' && std::toupper(ref[start_position + number_of_bases]) == 'G') {
				right_offset = 1;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 4]) == 'A' && std::toupper(ref[start_position + number_of_bases - 3]) == 'G') {
				right_offset = -2;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases]) == 'A' && std::toupper(ref[start_position + number_of_bases + 1]) == 'G') {
				right_offset = 2;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 5]) == 'A' && std::toupper(ref[start_position + number_of_bases - 4]) == 'G') {
				right_offset = -3;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases + 1]) == 'A' && std::toupper(ref[start_position + number_of_bases + 2]) == 'G') {
				right_offset = 3;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 6]) == 'A' && std::toupper(ref[start_position + number_of_bases - 5]) == 'G') {
				right_offset = -4;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases + 2]) == 'A' && std::toupper(ref[start_position + number_of_bases + 3]) == 'G') {
				right_offset = 4;
				found_right_offset = true;
			}

			if(found_left_offset && found_right_offset) {
				result_After->cigar[i].count -= left_offset;
				result_After->cigar[i-1].count += left_offset;
				result_After->cigar[i].count += right_offset;
				result_After->cigar[i+1].count -= right_offset;
			}
		}

	}
	if(result_After->cigar[i].op != 'I') {
		start_position += result_After->cigar[i].count;
	}
  }

  return result_After;
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchored(const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors, bool type) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);

  if (anchors.size() == 0) {
    return result;
  }

  result->cigar.clear();
  result->score = 0;
  int64_t offset = 0;
  int64_t firstQuery = anchors[0].rstart;
  int64_t refLen = 0;

  // Align between anchors.
  for (int64_t i = 0; i < (anchors.size() - 1); i++) {
	  int64_t start_ref = anchors[i].rstart + offset;
	  aligner_->Global(query + anchors[i].qstart, anchors[i+1].qend - anchors[i].qstart,
		                      ref + start_ref, anchors[i+1].rend - start_ref, type);

	  auto aln_result = aligner_->getResults();
	  int64_t cigar_length = 0;
	  auto left_part = ExtractCigarBetweenQueryCoords(aln_result->cigar,
		                                                    0,
		                                                    anchors[i+1].qstart - anchors[i].qstart, &cigar_length); // Leave next anchor for the next alignment.

	  refLen += cigar_length;
	  offset = refLen - (anchors[i+1].rstart - firstQuery);

	  result->cigar.insert(result->cigar.end(), left_part.begin(), left_part.end());
	  result->score += aln_result->score;
  }

  // Align the last anchor.
  aligner_->Global(query + anchors.back().qstart, anchors.back().qend - anchors.back().qstart,
                    ref + anchors.back().rstart + offset, anchors.back().rend - (anchors.back().rstart + offset), type);

  auto aln_result = aligner_->getResults();

  result->cigar.insert(result->cigar.end(), aln_result->cigar.begin(), aln_result->cigar.end());

  // Add the soft clippings at front and back.
  if ( anchors.front().qstart > 0) {
    result->cigar.insert(result->cigar.begin(), is::CigarOp('S', anchors.front().qstart));
  }
  if ((qlen - anchors.back().qend) > 0) {
    result->cigar.insert(result->cigar.end(), is::CigarOp('S', (qlen - anchors.back().qend)));
  }

  //  // Fill the other alignment info.
  result->edit_dist = EditDistFromExtCIGAR(result->cigar);
  result->position = is::AlignmentPosition(anchors.front().qstart, anchors.back().qend, anchors.front().rstart, anchors.back().rend);
  result->k = -1;
  result->rv = is::AlignmentReturnValue::OK;

  const int64_t MIN_INTRON_LEN = 10;
  for (auto& c: result->cigar) {
    if (c.op == 'D' && c.count >= MIN_INTRON_LEN) {
      c.op = 'N';
    }
  }

  if(result->cigar.size() < 2) {
	  return result;
  }

  int start_position = (result->position.tstart-anchors.front().qstart) + result->cigar[0].count;

  bool hasN = false;

  for(int i = 1; i < result->cigar.size()-1; i++) {
	int number_of_bases = result->cigar[i].count;
	if (result->cigar[i].op == 'N') {
		hasN = true;

		int left_offset = 0;
		bool found_left_offset = false;
		if(std::toupper(ref[start_position]) == 'C' && std::toupper(ref[start_position+1]) == 'T') {
			left_offset = 0;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position-1]) == 'C' && std::toupper(ref[start_position]) == 'T') {
			left_offset = -1;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position+1]) == 'C' && std::toupper(ref[start_position+2]) == 'T') {
			left_offset = 1;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position+2]) == 'C' && std::toupper(ref[start_position+3]) == 'T') {
			left_offset = 2;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position-2]) == 'C' && std::toupper(ref[start_position-1]) == 'T') {
			left_offset = -2;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position-3]) == 'C' && std::toupper(ref[start_position-2]) == 'T') {
			left_offset = -3;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position+3]) == 'C' && std::toupper(ref[start_position+4]) == 'T') {
			left_offset = 3;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position+4]) == 'C' && std::toupper(ref[start_position+5]) == 'T') {
			left_offset = 4;
			found_left_offset = true;
		} else if(std::toupper(ref[start_position-4]) == 'C' && std::toupper(ref[start_position-3]) == 'T') {
			left_offset = -4;
			found_left_offset = true;
		}
		int right_offset = 0;
		bool found_right_offset = false;
		if(std::toupper(ref[start_position + number_of_bases - 2]) == 'A' && std::toupper(ref[start_position + number_of_bases - 1]) == 'C') {
			right_offset = 0;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 3]) == 'A' && std::toupper(ref[start_position + number_of_bases - 2]) == 'C') {
			right_offset = -1;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 1]) == 'A' && std::toupper(ref[start_position + number_of_bases]) == 'C') {
			right_offset = 1;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 4]) == 'A' && std::toupper(ref[start_position + number_of_bases - 3]) == 'C') {
			right_offset = -2;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases]) == 'A' && std::toupper(ref[start_position + number_of_bases + 1]) == 'C') {
			right_offset = 2;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 5]) == 'A' && std::toupper(ref[start_position + number_of_bases - 4]) == 'C') {
			right_offset = -3;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases + 1]) == 'A' && std::toupper(ref[start_position + number_of_bases + 2]) == 'C') {
			right_offset = 3;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases - 6]) == 'A' && std::toupper(ref[start_position + number_of_bases - 5]) == 'C') {
			right_offset = -4;
			found_right_offset = true;
		} else if (std::toupper(ref[start_position + number_of_bases + 2]) == 'A' && std::toupper(ref[start_position + number_of_bases + 3]) == 'C') {
			right_offset = 4;
			found_right_offset = true;
		}

		if(found_left_offset && found_right_offset) {
			result->cigar[i].count -= left_offset;
			result->cigar[i-1].count += left_offset;
			result->cigar[i].count += right_offset;
			result->cigar[i+1].count -= right_offset;
		} else {
			left_offset = 0;
			found_left_offset = false;
			if(std::toupper(ref[start_position]) == 'G' && std::toupper(ref[start_position+1]) == 'T') {
				left_offset = 0;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position-1]) == 'G' && std::toupper(ref[start_position]) == 'T') {
				left_offset = -1;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position+1]) == 'G' && std::toupper(ref[start_position+2]) == 'T') {
				left_offset = 1;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position+2]) == 'G' && std::toupper(ref[start_position+3]) == 'T') {
				left_offset = 2;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position-2]) == 'G' && std::toupper(ref[start_position-1]) == 'T') {
				left_offset = -2;
				found_left_offset = true;
			}else if(std::toupper(ref[start_position-3]) == 'G' && std::toupper(ref[start_position-2]) == 'T') {
				left_offset = -3;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position+3]) == 'G' && std::toupper(ref[start_position+4]) == 'T') {
				left_offset = 3;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position+4]) == 'G' && std::toupper(ref[start_position+5]) == 'T') {
				left_offset = 4;
				found_left_offset = true;
			} else if(std::toupper(ref[start_position-4]) == 'G' && std::toupper(ref[start_position-3]) == 'T') {
				left_offset = -4;
				found_left_offset = true;
			}

			right_offset = 0;
			found_right_offset = false;
			if(std::toupper(ref[start_position + number_of_bases - 2]) == 'A' && std::toupper(ref[start_position + number_of_bases - 1]) == 'G') {
				right_offset = 0;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 3]) == 'A' && std::toupper(ref[start_position + number_of_bases - 2]) == 'G') {
				right_offset = -1;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 1]) == 'A' && std::toupper(ref[start_position + number_of_bases]) == 'G') {
				right_offset = 1;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 4]) == 'A' && std::toupper(ref[start_position + number_of_bases - 3]) == 'G') {
				right_offset = -2;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases]) == 'A' && std::toupper(ref[start_position + number_of_bases + 1]) == 'G') {
				right_offset = 2;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 5]) == 'A' && std::toupper(ref[start_position + number_of_bases - 4]) == 'G') {
				right_offset = -3;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases + 1]) == 'A' && std::toupper(ref[start_position + number_of_bases + 2]) == 'G') {
				right_offset = 3;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases - 6]) == 'A' && std::toupper(ref[start_position + number_of_bases - 5]) == 'G') {
				right_offset = -4;
				found_right_offset = true;
			} else if (std::toupper(ref[start_position + number_of_bases + 2]) == 'A' && std::toupper(ref[start_position + number_of_bases + 3]) == 'G') {
				right_offset = 4;
				found_right_offset = true;
			}
			if(found_left_offset && found_right_offset) {

				result->cigar[i].count -= left_offset;
				result->cigar[i-1].count += left_offset;
				result->cigar[i].count += right_offset;
				result->cigar[i+1].count -= right_offset;
			}
		}

	}
	if(result->cigar[i].op != 'I') {
		start_position += result->cigar[i].count;
	}
  }

  return result;
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchoredWithExtend(const char *query, int64_t qlen,
                                                                         const char *ref, int64_t rlen,
                                                                         const std::vector<AlignmentAnchor>& anchors,
                                                                         int32_t bandwidth, int32_t zdrop, bool type) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);

  if (anchors.size() == 0) {
    return result;
  }

  std::vector<AlignmentAnchor> updated_anchors = anchors;

  // Extend align last anchor forward to find max position.
  // Align the last anchor.
  aligner_->Extend(query + anchors.back().qstart, qlen - anchors.back().qstart,
                   ref + anchors.back().rstart, std::min(rlen - anchors.back().rstart, (qlen - anchors.back().qstart) * 2),
                   bandwidth, zdrop);
  auto ext_back = aligner_->getResults();
  // Get the correct coordinates from the alignment.
  int64_t max_q_pos_back = ext_back->max_q_pos + anchors.back().qstart + 1; // The "+1" because end coordinate is non-inclusive in GraphMap.
  int64_t max_t_pos_back = ext_back->max_t_pos + anchors.back().rstart + 1; // The max position is inclusive on the other hand.
  // If extend did not pan out (e.g. band is too narrow), do not extend.
  if (ext_back->max_q_pos >= 0 && ext_back->max_t_pos >= 0) {
    updated_anchors.back().qend = max_q_pos_back;
    updated_anchors.back().rend = max_t_pos_back;
  }

  // Extension of the front part.
  // Create a reverse copy (but not complemented) of the query and target.
  std::string rev_q_front, rev_t_front;
  std::reverse_copy(query , query + anchors.front().qend, std::back_inserter(rev_q_front));
  std::reverse_copy(ref + std::max((int64_t) 0, anchors.front().rend - 2 * anchors.front().qend),
                    ref + anchors.front().rend, std::back_inserter(rev_t_front));
  aligner_->Extend(rev_q_front.c_str(), rev_q_front.size(),
                   rev_t_front.c_str(), rev_t_front.size(),
                   bandwidth, zdrop);
  auto ext_front = aligner_->getResults();
  // Get the correct coordinates from the alignment.
  int64_t max_q_pos_front = anchors.front().qend - (ext_front->max_q_pos + 1); // The "+1" because end coordinate is non-inclusive in GraphMap.
  int64_t max_t_pos_front = anchors.front().rend - (ext_front->max_t_pos + 1); // The max position is inclusive on the other hand.
  // If extend did not pan out (e.g. band is too narrow), do not extend.
  // Added to fix crash when extending makes rstart bigger than rend
  if (ext_front->max_q_pos >= 0 && ext_front->max_t_pos >= 0 && max_t_pos_front < updated_anchors.front().rstart && max_q_pos_front < updated_anchors.front().qstart) {
    updated_anchors.front().qstart = max_q_pos_front;
    updated_anchors.front().rstart = max_t_pos_front;
  }

  // Align the updated coordinates.
  return GlobalAnchored(query, qlen, ref, rlen, updated_anchors, type);
}

}
