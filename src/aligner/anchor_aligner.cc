/*
 * anchor_aligner.cc
 *
 *  Created on: Aug 23, 2017
 *      Author: isovic
 */

#include "anchor_aligner.h"
#include "aligner_util.hpp"
#include <iostream>

//#include <chrono>

namespace is {

std::shared_ptr<AnchorAligner> createAnchorAligner(std::shared_ptr<is::AlignerBase> aligner) {
  return std::shared_ptr<AnchorAligner>(new AnchorAligner(aligner));
}

AnchorAligner::AnchorAligner(std::shared_ptr<is::AlignerBase> aligner) : aligner_(aligner) {

}

AnchorAligner::~AnchorAligner() {

}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalEndToEnd(int64_t abs_ref_id, std::shared_ptr<is::MinimizerIndex> index, const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);
  std::vector<AlignmentAnchor> final_anchors;
  if (anchors.size() > 0) {
    final_anchors.emplace_back(AlignmentAnchor(anchors.front().qstart, anchors.back().qend, anchors.front().rstart, anchors.back().rend));
  }
  return GlobalAnchored(abs_ref_id, index, query, qlen, ref, rlen, final_anchors, true);
}

bool FindRefOffsets(const char *ref, char first_base, char second_base, char third_base, char fourth_base, int64_t *left_offset, int64_t *right_offset, int64_t start_position, int number_of_bases) {
	bool found_left_offset = false;
	if(std::toupper(ref[start_position]) == first_base && std::toupper(ref[start_position+1]) == second_base) {
		*left_offset = 0;
		found_left_offset = true;
	} else if(std::toupper(ref[start_position-1]) == first_base && std::toupper(ref[start_position]) == second_base) {
		*left_offset = -1;
		found_left_offset = true;
	} else if(std::toupper(ref[start_position+1]) == first_base && std::toupper(ref[start_position+2]) == second_base) {
		*left_offset = 1;
		found_left_offset = true;
	} else if(std::toupper(ref[start_position+2]) == first_base && std::toupper(ref[start_position+3]) == second_base) {
		*left_offset = 2;
		found_left_offset = true;
	} else if(std::toupper(ref[start_position-2]) == first_base && std::toupper(ref[start_position-1]) == second_base) {
		*left_offset = -2;
		found_left_offset = true;
	} else if(std::toupper(ref[start_position-3]) == first_base && std::toupper(ref[start_position-2]) == second_base) {
		*left_offset = -3;
		found_left_offset = true;
	} else if(std::toupper(ref[start_position+3]) == first_base && std::toupper(ref[start_position+4]) == second_base) {
		*left_offset = 3;
		found_left_offset = true;
	} else if(std::toupper(ref[start_position+4]) == first_base && std::toupper(ref[start_position+5]) == second_base) {
		*left_offset = 4;
		found_left_offset = true;
	} else if(std::toupper(ref[start_position-4]) == first_base && std::toupper(ref[start_position-3]) == second_base) {
		*left_offset = -4;
		found_left_offset = true;
	}
//	else if(std::toupper(ref[start_position+5]) == first_base && std::toupper(ref[start_position+6]) == second_base) {
//		*left_offset = 5;
//		found_left_offset = true;
//	} else if(std::toupper(ref[start_position-5]) == first_base && std::toupper(ref[start_position-4]) == second_base) {
//		*left_offset = -5;
//		found_left_offset = true;
//	}
//	else if(std::toupper(ref[start_position+6]) == first_base && std::toupper(ref[start_position+7]) == second_base) {
//		*left_offset = 6;
//		found_left_offset = true;
//	} else if(std::toupper(ref[start_position-6]) == first_base && std::toupper(ref[start_position-5]) == second_base) {
//		*left_offset = -6;
//		found_left_offset = true;
//	}
//	else if(std::toupper(ref[start_position+7]) == first_base && std::toupper(ref[start_position+8]) == second_base) {
//		*left_offset = 7;
//		found_left_offset = true;
//	} else if(std::toupper(ref[start_position-7]) == first_base && std::toupper(ref[start_position-6]) == second_base) {
//		*left_offset = -7;
//		found_left_offset = true;
//	}

	bool found_right_offset = false;
	if(std::toupper(ref[start_position + number_of_bases - 2]) == third_base && std::toupper(ref[start_position + number_of_bases - 1]) == fourth_base) {
		*right_offset = 0;
		found_right_offset = true;
	} else if (std::toupper(ref[start_position + number_of_bases - 3]) == third_base && std::toupper(ref[start_position + number_of_bases - 2]) == fourth_base) {
		*right_offset = -1;
		found_right_offset = true;
	} else if (std::toupper(ref[start_position + number_of_bases - 1]) == third_base && std::toupper(ref[start_position + number_of_bases]) == fourth_base) {
		*right_offset = 1;
		found_right_offset = true;
	} else if (std::toupper(ref[start_position + number_of_bases - 4]) == third_base && std::toupper(ref[start_position + number_of_bases - 3]) == fourth_base) {
		*right_offset = -2;
		found_right_offset = true;
	} else if (std::toupper(ref[start_position + number_of_bases]) == third_base && std::toupper(ref[start_position + number_of_bases + 1]) == fourth_base) {
		*right_offset = 2;
		found_right_offset = true;
	} else if (std::toupper(ref[start_position + number_of_bases - 5]) == third_base && std::toupper(ref[start_position + number_of_bases - 4]) == fourth_base) {
		*right_offset = -3;
		found_right_offset = true;
	} else if (std::toupper(ref[start_position + number_of_bases + 1]) == third_base && std::toupper(ref[start_position + number_of_bases + 2]) == fourth_base) {
		*right_offset = 3;
		found_right_offset = true;
	} else if (std::toupper(ref[start_position + number_of_bases - 6]) == third_base && std::toupper(ref[start_position + number_of_bases - 5]) == fourth_base) {
		*right_offset = -4;
		found_right_offset = true;
	} else if (std::toupper(ref[start_position + number_of_bases + 2]) == third_base && std::toupper(ref[start_position + number_of_bases + 3]) == fourth_base) {
		*right_offset = 4;
		found_right_offset = true;
	}
//	else if (std::toupper(ref[start_position + number_of_bases - 7]) == third_base && std::toupper(ref[start_position + number_of_bases - 6]) == fourth_base) {
//		*right_offset = -5;
//		found_right_offset = true;
//	} else if (std::toupper(ref[start_position + number_of_bases + 3]) == third_base && std::toupper(ref[start_position + number_of_bases + 4]) == fourth_base) {
//		*right_offset = 5;
//		found_right_offset = true;
//	}
//	else if (std::toupper(ref[start_position + number_of_bases - 8]) == third_base && std::toupper(ref[start_position + number_of_bases - 7]) == fourth_base) {
//		*right_offset = -6;
//		found_right_offset = true;
//	} else if (std::toupper(ref[start_position + number_of_bases + 4]) == third_base && std::toupper(ref[start_position + number_of_bases + 5]) == fourth_base) {
//		*right_offset = 6;
//		found_right_offset = true;
//	}
//	else if (std::toupper(ref[start_position + number_of_bases - 9]) == third_base && std::toupper(ref[start_position + number_of_bases - 8]) == fourth_base) {
//		*right_offset = -7;
//		found_right_offset = true;
//	} else if (std::toupper(ref[start_position + number_of_bases + 5]) == third_base && std::toupper(ref[start_position + number_of_bases + 6]) == fourth_base) {
//		*right_offset = 7;
//		found_right_offset = true;
//	}

	return found_left_offset && found_right_offset;
}

int FindReadLeftOffset(const char *query, int left_offset_ref, int64_t start_position_read, std::stack<is::CigarOp> *cigar_stack) {
	int counter = -left_offset_ref;
	int read_offset = 0;

	while (!cigar_stack->empty()) {
		is::CigarOp c = cigar_stack->top();
		cigar_stack->pop();
		int count = c.count;

		if(c.op == 'N' || c.op == 'S') {
			cigar_stack->push(c);
			break;
		}

		if(c.op == '=' || c.op == 'X') {
			if(count > counter) {
				read_offset += counter;
				is:CigarOp new_c = is::CigarOp(c.op, count - counter);
				cigar_stack->push(new_c);
				break;
			} else if(count == counter) {
				read_offset += counter;
				break;
			} else {
				read_offset += count;
				counter -= count;
			}
		}

		if(c.op == 'D') {
			if(count > counter) {
				is::CigarOp new_c = is::CigarOp(c.op, count - counter);
				cigar_stack->push(new_c);
				break;
			} else if(count == counter) {
				break;
			} else {
				counter -= count;
			}
		}

		if(c.op == 'I') {
			read_offset += count;
		}
	}

	return read_offset;
}

int FindReadRightOffset(const char *query, int right_offset_ref, int64_t start_position_read, std::deque<is::CigarOp> *cigar_queue) {
	int counter = right_offset_ref;
	int read_offset = 0;

	while (!cigar_queue->empty()) {
		is::CigarOp c = cigar_queue->front();
		cigar_queue->pop_front();
		int count = c.count;

		if(c.op == 'N' || c.op == 'S') {
			cigar_queue->push_front(c);
			break;
		}

		if(c.op == '=' || c.op == 'X') {
			if(count > counter) {
				read_offset += counter;
				is:CigarOp new_c = is::CigarOp(c.op, count - counter);
				cigar_queue->push_front(new_c);
				break;
			} else if(count == counter) {
				read_offset += counter;
				break;
			} else {
				read_offset += count;
				counter -= count;
			}
		}

		if(c.op == 'D') {
			if(count > counter) {
				is::CigarOp new_c = is::CigarOp(c.op, count - counter);
				cigar_queue->push_front(new_c);
				break;
			} else if(count == counter) {
				break;
			} else {
				counter -= count;
			}
		}

		if(c.op == 'I') {
			read_offset += count;
		}
	}

	return read_offset;
}

void AnchorAligner::AdjustEnds(int left_offset_ref, int right_offset_ref, const char *query, const char *ref, int64_t *start_position_ref, int64_t *start_position_read, int number_of_bases, std::stack<is::CigarOp> *cigar_stack, std::deque<is::CigarOp> *cigar_queue, bool type) {
	int left_offset_read = left_offset_ref > 0 ? 0 : FindReadLeftOffset(query, left_offset_ref, *start_position_read, cigar_stack);
	int right_offset_read = right_offset_ref < 0 ? 0 : FindReadRightOffset(query, right_offset_ref, *start_position_read + number_of_bases, cigar_queue);

//	std::cout << "lr " << left_offset_read << " " << left_offset_ref << std::endl;
//	std::cout << "rr " << right_offset_read << " " << right_offset_ref << std::endl;

	if (left_offset_ref >= 0 && right_offset_ref <= 0) {
//		std::cout << "1" << std::endl;
		if(left_offset_ref > 0) {
			is::CigarOp c_left = is::CigarOp('D', left_offset_ref);
			cigar_stack->push(c_left);
		}
		is::CigarOp c_gap = is::CigarOp('N', number_of_bases + right_offset_ref + (-left_offset_ref));
		cigar_stack->push(c_gap);
		if(-right_offset_ref > 0) {
			is::CigarOp c_right = is::CigarOp('D', -right_offset_ref);
			cigar_stack->push(c_right);
		}

		*start_position_ref += number_of_bases;

	} else if (left_offset_ref <= 0 && right_offset_ref >= 0) {
//		std::cout << "2" << std::endl;
		if(left_offset_read > 0) {
			is::CigarOp c_left = is::CigarOp('I', left_offset_read);
			cigar_stack->push(c_left);
		}
		is::CigarOp c_gap = is::CigarOp('N', number_of_bases + right_offset_ref + (-left_offset_ref));
		cigar_stack->push(c_gap);
		if(right_offset_read > 0) {
			is::CigarOp c_right = is::CigarOp('I', right_offset_read);
			cigar_stack->push(c_right);
		}

		*start_position_ref += (number_of_bases + right_offset_ref);
		*start_position_read += (right_offset_read);

	} else if (left_offset_ref >= 0 && right_offset_ref >= 0) {
//		std::cout << "3" << std::endl;
		std::string ref_string;
		for(int i = 0; i < left_offset_ref; i++) {
			ref_string.push_back(ref[*start_position_ref + i]);
		}
		std::string read_string;
		for(int i = 0; i < right_offset_read; i++) {
			read_string.push_back(query[*start_position_read + i]);
		}
		const char* ref_sub = ref_string.c_str();
		const char* read_sub = read_string.c_str();

		aligner_->Global(read_sub, right_offset_read, ref_sub, left_offset_ref, type);
		auto aln_result_sub = aligner_->getResults();

		for (auto& c: aln_result_sub->cigar) {
			cigar_stack->push(c);
		}

		is::CigarOp c_gap = is::CigarOp('N', number_of_bases + right_offset_ref + (-left_offset_ref));
		cigar_stack->push(c_gap);

		*start_position_ref += (number_of_bases + right_offset_ref);
		*start_position_read += (right_offset_read);

	} else if (left_offset_ref <= 0 && right_offset_ref <= 0) {
//		std::cout << "4" << std::endl;
		std::string ref_string;
		for(int i = right_offset_ref; i < 0; i++) {
			ref_string.push_back(ref[*start_position_ref + number_of_bases - i]);
		}
		std::string read_string;
		for(int i = -left_offset_read; i < 0; i++) {
			read_string.push_back(query[*start_position_read + left_offset_read]);
		}
		const char* ref_sub = ref_string.c_str();
		const char* read_sub = read_string.c_str();

		aligner_->Global(read_sub, left_offset_read, ref_sub, -right_offset_ref, type);
		auto aln_result_sub = aligner_->getResults();

		is::CigarOp c_gap = is::CigarOp('N', number_of_bases + right_offset_ref + (-left_offset_ref));
		cigar_stack->push(c_gap);

		for (auto& c: aln_result_sub->cigar) {
			cigar_stack->push(c);
		}

		*start_position_ref += (number_of_bases);
	}
}

std::shared_ptr<AlignmentResult> AnchorAligner::CreateAlignmentResult(int rstart, int rend, int qstart, int qend, std::vector<is::CigarOp> rez) {
	  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);

	  for(is::CigarOp &c: rez) {
		  result->cigar.push_back(c);
	  }

	  result->edit_dist = EditDistFromExtCIGAR(result->cigar);
	  result->position = is::AlignmentPosition(qstart, qend, rstart, rend);
	  result->k = -1;
	  result->rv = is::AlignmentReturnValue::OK;

	  return result;
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchoredWithClipping(const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);

  if (anchors.size() == 0) {
    return result;
  }

  result->cigar.clear();
  result->score = 0;

  int64_t cigar_length = 0;
  int64_t cigar_length_q = 0;

  // Align between anchors.
  for (int64_t i = 0; i < (anchors.size() - 1); i++) {
    aligner_->Global(query + anchors[i].qstart, anchors[i+1].qend - anchors[i].qstart,
                      ref + anchors[i].rstart, anchors[i+1].rend - anchors[i].rstart, false);

    auto aln_result = aligner_->getResults();

    auto left_part = ExtractCigarBetweenQueryCoords(aln_result->cigar,
                                                    0,
                                                    anchors[i+1].qstart - anchors[i].qstart, &cigar_length, &cigar_length_q); // Leave next anchor for the next alignment.

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
//  if ( anchors.front().qstart > 0) {
//    result->cigar.insert(result->cigar.begin(), is::CigarOp('S', anchors.front().qstart));
//  }
//  if ((qlen - anchors.back().qend) > 0) {
//    result->cigar.insert(result->cigar.end(), is::CigarOp('S', (qlen - anchors.back().qend)));
//  }

  for (auto& c: result->cigar) {
	  if(c.op != 'N' && c.op != 'D') {
		  cigar_length_q += c.count;
	  }
  }

  result->cigar.push_back(is::CigarOp('S', (qlen- cigar_length_q)));

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
	  return result;
	}

	int64_t start_position_ref = result_After->position.tstart;
	int64_t start_position_read = result_After->position.qstart;

	std::stack<is::CigarOp> cigar_stack;
	std::deque<is::CigarOp> cigar_queue;

	for(int i = 0; i < result_After->cigar.size(); i++) {
	  cigar_queue.push_back(result_After->cigar[i]);
	}

	while (!cigar_queue.empty()) {
		is::CigarOp cigar_op = cigar_queue.front();
		cigar_queue.pop_front();
		int number_of_bases = cigar_op.count;

		if (cigar_op.op == 'S') {
			cigar_stack.push(cigar_op);
			continue;
		}

		if (cigar_op.op == 'N') {
			int64_t left_offset_ref = 0;
			int64_t right_offset_ref = 0;
			bool found_base_pairs = FindRefOffsets(ref, 'G', 'T', 'A', 'G', &left_offset_ref, &right_offset_ref, start_position_ref, number_of_bases);
			if(found_base_pairs && (left_offset_ref != 0 || right_offset_ref != 0)) {
				AdjustEnds(left_offset_ref, right_offset_ref, query, ref, &start_position_ref, &start_position_read, number_of_bases, &cigar_stack, &cigar_queue, 1);
			} else {
				int64_t left_offset_ref = 0;
				int64_t right_offset_ref = 0;
				bool found_base_pairs = FindRefOffsets(ref, 'C', 'T', 'A', 'C', &left_offset_ref, &right_offset_ref, start_position_ref, number_of_bases);
				if(found_base_pairs && (left_offset_ref != 0 || right_offset_ref != 0)) {
					AdjustEnds(left_offset_ref, right_offset_ref, query, ref, &start_position_ref, &start_position_read, number_of_bases, &cigar_stack, &cigar_queue, 1);
				} else {
					cigar_stack.push(cigar_op);

					if(cigar_op.op != 'I') {
						start_position_ref += cigar_op.count;
					}
					if(cigar_op.op != 'D' && cigar_op.op != 'N') {
						start_position_read += cigar_op.count;
					}
				}
			}
		} else {
			cigar_stack.push(cigar_op);

			if(cigar_op.op != 'I') {
				start_position_ref += cigar_op.count;
			}
			if(cigar_op.op != 'D' && cigar_op.op != 'N') {
				start_position_read += cigar_op.count;
			}
		}
	}

    std::stack<is::CigarOp> tmp_stack;

    is::CigarOp previous_op = cigar_stack.top();
    cigar_stack.pop();

    while(!cigar_stack.empty()) {
  	  is::CigarOp tmp_op = cigar_stack.top();
  	  cigar_stack.pop();
  	  if(tmp_op.op == previous_op.op) {
  		  previous_op = is::CigarOp(previous_op.op, previous_op.count + tmp_op.count);
  	  } else {
  		  tmp_stack.push(previous_op);
  		  previous_op = tmp_op;
  	  }
    }
    tmp_stack.push(previous_op);

    result_After->cigar.clear();

    while(!tmp_stack.empty()) {
  	  is::CigarOp c = tmp_stack.top();
  	  tmp_stack.pop();
  	  result_After->cigar.push_back(c);
    }

  return result_After;
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchored(int64_t abs_ref_id, std::shared_ptr<is::MinimizerIndex> index, const char *query, int64_t qlen, const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors, bool type) {
  auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);

  if (anchors.size() == 0) {
    return result;
  }

  result->cigar.clear();
  result->score = 0;

  int64_t offset = 0;
  int64_t firstQuery = anchors[0].rstart;
  int64_t refLen = 0;

  int64_t offset_q = 0;
  int64_t firstQuery_q = anchors[0].qstart;
  int64_t refLen_q = 0;

  // Align between anchors.
  for (int64_t i = 0; i < (anchors.size() - 1); i++) {
	  int64_t start_ref = anchors[i].rstart + offset;
	  int64_t start_ref_q = anchors[i].qstart + offset_q;

	  if(((anchors[i+1].qend - start_ref_q) > 85000 || (anchors[i+1].rend - start_ref) > 85000) && type == 0) {
		  return result;
	  }


	  aligner_->Global(query + start_ref_q, anchors[i+1].qend - start_ref_q,
		                      ref + start_ref, anchors[i+1].rend - start_ref, type);

	  auto aln_result = aligner_->getResults();
	  int64_t cigar_length = 0;
	  int64_t cigar_length_q = 0;

	  auto left_part = ExtractCigarBetweenQueryCoords(aln_result->cigar,
		                                                    0,
		                                                    anchors[i+1].qstart - anchors[i].qstart, &cigar_length, &cigar_length_q); // Leave next anchor for the next alignment.

	  refLen += cigar_length;
	  refLen_q += cigar_length_q;

	  offset = refLen - (anchors[i+1].rstart - firstQuery);
	  offset_q = refLen_q - (anchors[i+1].qstart - firstQuery_q);

	  result->cigar.insert(result->cigar.end(), left_part.begin(), left_part.end());
	  result->score += aln_result->score;
  }

  // Align the last anchor.
  aligner_->Global(query + anchors.back().qstart + offset_q, anchors.back().qend - (anchors.back().qstart + offset_q),
                    ref + anchors.back().rstart + offset, anchors.back().rend - (anchors.back().rstart + offset), type);

  auto aln_result = aligner_->getResults();

  result->cigar.insert(result->cigar.end(), aln_result->cigar.begin(), aln_result->cigar.end());

  const int64_t MIN_INTRON_LEN = 10;

  int64_t s_min_value = 13;
  int64_t exon_min_value = 15;
  int64_t window = 8000;
  int64_t windowBase = 8000;
  int64_t minimumWindow = 800;

  int backS = qlen - anchors.back().qend;

  std::vector<is::CigarOp> tmp_cigar;
  int64_t len_tmp = 0;
  int64_t current_len = 0;
  int64_t current_ref_len = 0;
  int64_t minimum_found_exon = 0;
  bool isFirst = true;

  int64_t lower_bound = index->get_reference_starting_pos()[abs_ref_id];
  int64_t ref_data_len = index->get_reference_lengths()[abs_ref_id];
  int64_t upper_bound = lower_bound + ref_data_len;

  double threshold = 0.85;

  window = std::min(windowBase, upper_bound - anchors.back().rend);

  if(backS > s_min_value && window > minimumWindow) {

	  bool isAdapter = false;
	  int countA = 0;
	  int countT = 0;

	  for(int i = 0; i < anchors.front().qstart; i++) {
		  if(query[i] == 'A') {
			  countA += 1;
			  countT = 0;
		  } else if(query[i] == 'T') {
			  countT += 1;
			  countA = 0;
		  } else {
			  countA = 0;
			  countT = 0;
		  }

		  if(countA >= 10 || countT >= 10) {
			  isAdapter = true;
		  }
	  }

	  if(!isAdapter && backS < 100) {
		  aligner_->Global(query + anchors.back().qend, backS, ref + anchors.back().rend, window, false);

		  auto aln_result_tmp = aligner_->getResults();

		  for (auto& c: aln_result_tmp->cigar) {
			  if ((c.op == 'N' || c.op == 'D') && c.count >= MIN_INTRON_LEN) {
				  if(isFirst) {
					  isFirst = false;
					  tmp_cigar.push_back(c);
				  } else {
					  double rez = (double) minimum_found_exon / (double) len_tmp;
					  if((minimum_found_exon < exon_min_value) && rez > threshold) {
						  break;
					  } else {
						  current_len += len_tmp;
						  len_tmp = 0;
						  minimum_found_exon = 0;
						  for (auto& c1: tmp_cigar) {
							  if(c1.op != 'I') {
								  current_ref_len += c1.count;
							  }
							  result->cigar.push_back(c1);
						  }
						  tmp_cigar.clear();
						  tmp_cigar.push_back(c);
					  }
				  }
			  } else {
				  isFirst = false;
				if(c.op != 'D') {
					len_tmp += c.count;
					if(c.op == '=') {
						minimum_found_exon += c.count;
					}
				}
				tmp_cigar.push_back(c);
			  }
		  }
		  double rez = (double) minimum_found_exon / (double) len_tmp;
		  if(minimum_found_exon >= exon_min_value  && rez > threshold && aln_result_tmp->cigar.size() > 0) {
			  current_len += len_tmp;
			  len_tmp = 0;
			  minimum_found_exon = 0;
			  for (auto& c1: tmp_cigar) {
				  if(c1.op != 'I') {
					  current_ref_len += c1.count;
				  }
				  result->cigar.push_back(c1);
			  }
			  tmp_cigar.clear();
		  }
	  }
  }

  window = std::min(windowBase, anchors.front().rstart - lower_bound);

  int64_t frontS = anchors.front().qstart;

  tmp_cigar.clear();
  len_tmp = 0;
  minimum_found_exon = 0;
  int64_t current_len2 = 0;
  int64_t current_ref_len2 = 0;
  isFirst = true;

  if(frontS > s_min_value && window > minimumWindow) {

	  std::string read_String;

	  bool isAdapter = false;
	  int countA = 0;
	  int countT = 0;

	  for(int i = 0; i < anchors.front().qstart; i++) {
		  read_String.insert(0, 1, query[i]);
		  if(query[i] == 'A') {
			  countA += 1;
			  countT = 0;
		  } else if(query[i] == 'T') {
			  countT += 1;
			  countA = 0;
		  } else {
			  countA = 0;
			  countT = 0;
		  }

		  if(countA >= 10 || countT >= 10) {
			  isAdapter = true;
		  }
	  }

	  if(!isAdapter && frontS < 100) {
		  std::string ref_String;

		  for(int i = anchors.front().rstart-window; i < anchors.front().rstart; i++) {
			  ref_String.insert(0, 1, ref[i]);
		  }

		  aligner_->Global(read_String.c_str(), read_String.size(), ref_String.c_str(), ref_String.size(), false);

		  auto aln_result_tmp = aligner_->getResults();

		  for (auto& c: aln_result_tmp->cigar) {
			  if ((c.op == 'N' || c.op == 'D') && c.count >= MIN_INTRON_LEN) {
				  if(isFirst) {
					  isFirst = false;
					  tmp_cigar.push_back(c);
				  } else {
					  double rez = (double) minimum_found_exon / (double) len_tmp;
					  if(minimum_found_exon < exon_min_value && rez > threshold) {
						  break;
					  } else {
						  current_len2 += len_tmp;
						  len_tmp = 0;
						  minimum_found_exon = 0;
						  for (auto& c1: tmp_cigar) {
							  if(c1.op != 'I') {
								  current_ref_len2 += c1.count;
							  }
							  result->cigar.insert(result->cigar.begin(), c1);
						  }
						  tmp_cigar.clear();
						  tmp_cigar.push_back(c);
					  }
				  }
			  } else {
				  isFirst = false;
				if(c.op != 'D') {
					len_tmp += c.count;
					if(c.op == '=') {
						minimum_found_exon += c.count;
					}
				}
				tmp_cigar.push_back(c);
			  }
		  }

		  double rez = (double) minimum_found_exon / (double) len_tmp;
		  if(minimum_found_exon >= exon_min_value  && rez > threshold && aln_result_tmp->cigar.size() > 0) {
			  current_len2 += len_tmp;
			  len_tmp = 0;
			  minimum_found_exon = 0;
			  for (auto& c1: tmp_cigar) {
				  if(c1.op != 'I') {
					  current_ref_len2 += c1.count;
				  }
				  result->cigar.insert(result->cigar.begin(), c1);
			  }
			  tmp_cigar.clear();
		  }

  	  } else {
  	  }
  }

  // Add the soft clippings at front and back.
  if ( (anchors.front().qstart-current_len2) > 0) {
    result->cigar.insert(result->cigar.begin(), is::CigarOp('S', anchors.front().qstart - current_len2));
  }
  if ((qlen - (anchors.back().qend+current_len)) > 0) {
    result->cigar.insert(result->cigar.end(), is::CigarOp('S', (qlen - (anchors.back().qend+current_len))));
  }

  //  // Fill the other alignment info.
  result->edit_dist = EditDistFromExtCIGAR(result->cigar);
  result->position = is::AlignmentPosition(anchors.front().qstart-current_len2, anchors.back().qend+current_len, anchors.front().rstart-current_ref_len2, anchors.back().rend + current_ref_len);
  result->k = -1;
  result->rv = is::AlignmentReturnValue::OK;

  for (auto& c: result->cigar) {
    if (c.op == 'D' && c.count >= MIN_INTRON_LEN) {
      c.op = 'N';
    }
  }

  if(result->cigar.size() < 2) {
	  return result;
  }

  int64_t start_position_ref = result->position.tstart;
  int64_t start_position_read = result->position.qstart;

  std::stack<is::CigarOp> cigar_stack;
  std::deque<is::CigarOp> cigar_queue;

  for(int i = 0; i < result->cigar.size(); i++) {
	  cigar_queue.push_back(result->cigar[i]);
  }

  while (!cigar_queue.empty()) {
	is::CigarOp cigar_op = cigar_queue.front();
	cigar_queue.pop_front();
	int number_of_bases = cigar_op.count;

	if (cigar_op.op == 'S') {
		cigar_stack.push(cigar_op);
		continue;
	}

	if (cigar_op.op == 'N') {
		int64_t left_offset_ref = 0;
		int64_t right_offset_ref = 0;
		bool found_base_pairs = FindRefOffsets(ref, 'G', 'T', 'A', 'G', &left_offset_ref, &right_offset_ref, start_position_ref, number_of_bases);
		if(found_base_pairs && (left_offset_ref != 0 || right_offset_ref != 0)) {
//			std::cout << "prvi" << std::endl;
//			std::cout << ref[start_position_ref - 3] << ref[start_position_ref - 2] << ref[start_position_ref - 1] << ref[start_position_ref] << ref[start_position_ref + 1] << ref[start_position_ref + 2] << ref[start_position_ref + 3] << ref[start_position_ref + 4] << std::endl;
//			std::cout << ref[start_position_ref + number_of_bases - 4] << ref[start_position_ref + number_of_bases - 3] << ref[start_position_ref + number_of_bases - 2] << ref[start_position_ref + number_of_bases - 1] << ref[start_position_ref + number_of_bases] << ref[start_position_ref + number_of_bases + 1] << ref[start_position_ref + number_of_bases + 2] << ref[start_position_ref + number_of_bases + 3] << ref[start_position_ref + number_of_bases + 4] << std::endl;
//			std::cout << left_offset_ref << " " << right_offset_ref << std::endl;
			AdjustEnds(left_offset_ref, right_offset_ref, query, ref, &start_position_ref, &start_position_read, number_of_bases, &cigar_stack, &cigar_queue, 1);
		} else {
			int64_t left_offset_ref = 0;
			int64_t right_offset_ref = 0;
			bool found_base_pairs = FindRefOffsets(ref, 'C', 'T', 'A', 'C', &left_offset_ref, &right_offset_ref, start_position_ref, number_of_bases);
			if(found_base_pairs && (left_offset_ref != 0 || right_offset_ref != 0)) {
//				std::cout << "drugi" << std::endl;
//				std::cout << ref[start_position_ref - 3] << ref[start_position_ref - 2] << ref[start_position_ref - 1] << ref[start_position_ref] << ref[start_position_ref + 1] << ref[start_position_ref + 2] << ref[start_position_ref + 3] << ref[start_position_ref + 4] << std::endl;
//				std::cout << ref[start_position_ref + number_of_bases - 4] << ref[start_position_ref + number_of_bases - 3] << ref[start_position_ref + number_of_bases - 2] << ref[start_position_ref + number_of_bases - 1] << ref[start_position_ref + number_of_bases] << ref[start_position_ref + number_of_bases + 1] << ref[start_position_ref + number_of_bases + 2] << ref[start_position_ref + number_of_bases + 3] << ref[start_position_ref + number_of_bases + 4] << std::endl;
//				std::cout << left_offset_ref << " " << right_offset_ref << std::endl;
				AdjustEnds(left_offset_ref, right_offset_ref, query, ref, &start_position_ref, &start_position_read, number_of_bases, &cigar_stack, &cigar_queue, 1);
			} else {
				cigar_stack.push(cigar_op);

				if(cigar_op.op != 'I') {
					start_position_ref += cigar_op.count;
				}
				if(cigar_op.op != 'D' && cigar_op.op != 'N') {
					start_position_read += cigar_op.count;
				}
			}
		}
	} else {
		cigar_stack.push(cigar_op);

		if(cigar_op.op != 'I') {
			start_position_ref += cigar_op.count;
		}
		if(cigar_op.op != 'D' && cigar_op.op != 'N') {
			start_position_read += cigar_op.count;
		}
	}
  }

  std::stack<is::CigarOp> tmp_stack;

  is::CigarOp previous_op = cigar_stack.top();
  cigar_stack.pop();

  while(!cigar_stack.empty()) {
	  is::CigarOp tmp_op = cigar_stack.top();
	  cigar_stack.pop();
	  if(tmp_op.op == previous_op.op) {
		  previous_op = is::CigarOp(previous_op.op, previous_op.count + tmp_op.count);
	  } else {
		  tmp_stack.push(previous_op);
		  previous_op = tmp_op;
	  }
  }
  tmp_stack.push(previous_op);

  result->cigar.clear();

  while(!tmp_stack.empty()) {
	  is::CigarOp c = tmp_stack.top();
	  tmp_stack.pop();
	  result->cigar.push_back(c);
  }

  return result;
}

std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchoredWithExtend(int64_t abs_ref_id, std::shared_ptr<is::MinimizerIndex> index,
																	   const char *query, int64_t qlen,
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
  return GlobalAnchored(abs_ref_id, index, query, qlen, ref, rlen, updated_anchors, type);
}

}
