/*
 * cigargen.cc
 *
 *  Created on: Aug 28, 2014
 *      Author: ivan
 */

#include "alignment/cigargen.h"
#include "utility/utility_general.h"



std::string AlignmentToCigar(unsigned char *alignment, int alignmentLength, bool extended_format) {
  if (alignment == NULL || alignmentLength == 0)
    return (std::string("*-"));

  std::stringstream ss;

  char* cigar = NULL;
  if (extended_format == true) {
    AlignmentToExtendedCigar(alignment, alignmentLength, &cigar);
  } else {
    AlignmentToBasicCigar(alignment, alignmentLength, &cigar);
  }
  ss << cigar;
  if (cigar)
    free(cigar);

  return ss.str();
}

std::string ReverseCigarString(std::string &cigar) {
  std::stringstream ret_cigar;

  int64_t num_cigarops = 0;
  for (int64_t i=0; i<((int64_t) cigar.size()); i++)
    if (cigar[i] < '0' || cigar[i] > '9')
      num_cigarops += 1;

  std::vector<std::string> cigarops(num_cigarops, "");

  std::string current_cigar = "";
  for (int64_t i=0; i<((int64_t) cigar.size()); i++) {
    current_cigar += cigar[i];
    if (cigar[i] < '0' || cigar[i] > '9') {
      cigarops.push_back(current_cigar);
      current_cigar = "";
    }
  }

  for (int64_t i=0; i<((int64_t) cigarops.size()); i++) {
    ret_cigar << cigarops[cigarops.size() - i - 1];
  }

  cigarops.clear();

  return ret_cigar.str();
}

int ConvertInsertionsToClipping(unsigned char* alignment, int alignmentLength) {
  for (int64_t i=0; i<alignmentLength; i++) {
    if (alignment[i] == EDLIB_I)
      alignment[i] = EDLIB_S;
    else break;
  }
  for (int64_t i=(((int64_t) alignmentLength) - 1); i>=0; i--) {
    if (alignment[i] == EDLIB_I)
      alignment[i] = EDLIB_S;
    else break;
  }

  return 0;
}

int CountClippedBases(unsigned char* alignment, int alignmentLength, int64_t *ret_num_clipped_front, int64_t *ret_num_clipped_back) {
  if (ret_num_clipped_front == NULL || ret_num_clipped_back == NULL)
    return 1;

  int64_t num_clipped_front = 0, num_clipped_back = 0;
  for (int64_t i=0; i<alignmentLength; i++) {
    if (alignment[i] == EDLIB_S || alignment[i] == EDLIB_I)
      num_clipped_front += 1;
    else break;
  }
  for (int64_t i=(((int64_t) alignmentLength) - 1); i>=0; i--) {
    if (alignment[i] == EDLIB_S || alignment[i] == EDLIB_I)
      num_clipped_back += 1;
    else break;
  }

  *ret_num_clipped_front = num_clipped_front;
  *ret_num_clipped_back = num_clipped_back;

  return 0;
}

std::vector<unsigned char> FixAlignment(unsigned char* alignment, int alignmentLength) {
  int move_type = -1, prev_move_type = -1;
  int64_t num_moves = 0, prev_num_moves = 0;
  std::vector<unsigned char> new_alignment;
  new_alignment.reserve(alignmentLength);

  for (int64_t i=0; i<alignmentLength; i++) {

  }

  return new_alignment;
}

std::vector<unsigned char> FixAlignment1(unsigned char* alignment, int alignmentLength) {
  int move_type = -1, prev_move_type = -1;
  int64_t num_moves = 0, prev_num_moves = 0;
  std::vector<unsigned char> new_alignment;
  new_alignment.reserve(alignmentLength);
//  return new_alignment;

  for (int64_t i=0; i<alignmentLength; i++) {
//    printf ("[%ld] prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", i, prev_num_moves, num_moves, prev_move_type, move_type);
//    fflush(stdout);

    if (i == 0) {
      // This is the first op.
      move_type = alignment[i];
      num_moves = 1;
    } else if (alignment[i] != alignment[i-1]) {
      // Operation changed. Check if there was another one before it, or was the current streak the first one in the alignment.
      if (prev_move_type != -1) {
        if ((move_type == EDLIB_D && prev_move_type == EDLIB_I) || (move_type == EDLIB_I && prev_move_type == EDLIB_D)) {
          if (prev_num_moves < num_moves) {
            std::vector<int> temp_op(prev_num_moves, EDLIB_X);
            new_alignment.insert(new_alignment.end(), temp_op.begin(), temp_op.end());
//            printf ("-> (1) [%ld] prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", i, prev_num_moves, num_moves, prev_move_type, move_type);
//            fflush(stdout);
            num_moves = num_moves - prev_num_moves;

          } else if (prev_num_moves > num_moves) {
            std::vector<int> temp_op(prev_num_moves - num_moves, prev_move_type);
            new_alignment.insert(new_alignment.end(), temp_op.begin(), temp_op.end());
//            printf ("-> (2) [%ld] prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", i, prev_num_moves, num_moves, prev_move_type, move_type);
//            fflush(stdout);
            move_type = EDLIB_X;

          } else {
            std::vector<int> temp_op(prev_num_moves, EDLIB_X);
            new_alignment.insert(new_alignment.end(), temp_op.begin(), temp_op.end());
//            printf ("-> (3) [%ld] prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", i, prev_num_moves, num_moves, prev_move_type, move_type);
//            fflush(stdout);
            num_moves = 0;
            move_type = -1;
          }
        } else {
          std::vector<int> temp_op(prev_num_moves, prev_move_type);
          new_alignment.insert(new_alignment.end(), temp_op.begin(), temp_op.end());
//          printf ("-> (4) [%ld] prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", i, prev_num_moves, num_moves, prev_move_type, move_type);
//          fflush(stdout);
        }
      } else {
//        std::vector<int> temp_op(prev_num_moves, prev_move_type);
//        new_alignment.insert(new_alignment.end(), temp_op.begin(), temp_op.end());
//        printf ("-> (4) [%ld] prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", i, prev_num_moves, num_moves, prev_move_type, move_type);
//        fflush(stdout);
      }
      prev_move_type = move_type;
      prev_num_moves = num_moves;
      move_type = alignment[i];
      num_moves = 1;
    } else {
      num_moves += 1;
    }
  }
  if (alignmentLength > 0) {
//    std::vector<int> temp_op(prev_num_moves, prev_move_type);
//    new_alignment.insert(new_alignment.end(), temp_op.begin(), temp_op.end());
    if ((move_type == EDLIB_D && prev_move_type == EDLIB_I) || (move_type == EDLIB_I && prev_move_type == EDLIB_D)) {
      if (prev_num_moves < num_moves) {
        std::vector<int> temp_op(prev_num_moves, EDLIB_X);
        new_alignment.insert(new_alignment.end(), temp_op.begin(), temp_op.end());
//        printf ("-> (1) prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", prev_num_moves, num_moves, prev_move_type, move_type);
//        fflush(stdout);
        num_moves = num_moves - prev_num_moves;
        std::vector<int> temp_op2(num_moves, move_type);
        new_alignment.insert(new_alignment.end(), temp_op2.begin(), temp_op2.end());

      } else if (prev_num_moves > num_moves) {
        std::vector<int> temp_op(prev_num_moves - num_moves, prev_move_type);
        new_alignment.insert(new_alignment.end(), temp_op.begin(), temp_op.end());
//        printf ("-> (2) prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", prev_num_moves, num_moves, prev_move_type, move_type);
//        fflush(stdout);
        move_type = EDLIB_X;
        std::vector<int> temp_op2(num_moves, move_type);
        new_alignment.insert(new_alignment.end(), temp_op2.begin(), temp_op2.end());

      } else {
        std::vector<int> temp_op(prev_num_moves, EDLIB_X);
        new_alignment.insert(new_alignment.end(), temp_op.begin(), temp_op.end());
//        printf ("-> (3) prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", prev_num_moves, num_moves, prev_move_type, move_type);
//        fflush(stdout);
        num_moves = 0;
        move_type = -1;
      }
    } else {
      std::vector<int> temp_op1(prev_num_moves, prev_move_type);
      new_alignment.insert(new_alignment.end(), temp_op1.begin(), temp_op1.end());

      std::vector<int> temp_op2(num_moves, move_type);
      new_alignment.insert(new_alignment.end(), temp_op2.begin(), temp_op2.end());

//      printf ("-> (4) [%ld] prev_num_moves = %ld, num_moves == %ld\tprev_move_type = %d\tmove_type = %d\n", prev_num_moves, num_moves, prev_move_type, move_type);
//      fflush(stdout);
    }

  }

  return new_alignment;
}

int AlignmentToBasicCigar(unsigned char* alignment, int alignmentLength,
                          char** cigar_) {

    unsigned char *alignment_with_clipping = alignment;

    if (alignment[0] == EDLIB_I || alignment[alignmentLength-1] == EDLIB_I) {
      alignment_with_clipping = (unsigned char *) malloc(sizeof(unsigned char) * (alignmentLength + 1));
      memcpy(alignment_with_clipping, alignment, alignmentLength);
      alignment_with_clipping[alignmentLength] = '\0';

      for (int64_t i=0; i<alignmentLength; i++) {
        if (alignment_with_clipping[i] == EDLIB_I)
          alignment_with_clipping[i] = EDLIB_S;
        else break;
      }
      for (int64_t i=(((int64_t) alignmentLength) - 1); i>=0; i--) {
        if (alignment_with_clipping[i] == EDLIB_I)
          alignment_with_clipping[i] = EDLIB_S;
        else break;
      }
    }

    std::vector<char>* cigar = new std::vector<char>();
    unsigned char lastMove = -1;  // Code of last move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
      char alignment_char = 0;
      if (i < alignmentLength) {
        alignment_char = alignment_with_clipping[i];
        /// Replace the 'X' operation with 'M'.
        if (alignment_char == 3)
          alignment_char = 0;
      }

        // if new sequence of same moves started
        if (i == alignmentLength || alignment_char != lastMove) {
            if (i > 0) {  // if previous sequence of same moves ended
                // Write number of moves to cigar string.
                int numDigits = 0;
                for (; numOfSameMoves; numOfSameMoves /= 10) {
                    cigar->push_back('0' + numOfSameMoves % 10);
                    numDigits++;
                }
                reverse(cigar->end() - numDigits, cigar->end());
                // Write code of move to cigar string.
                if (lastMove == 0) {
//                    cigar->push_back('=');
                    cigar->push_back('M');
                } else if (lastMove == 1) {
                    cigar->push_back('I');
                } else if (lastMove == 2) {
                    cigar->push_back('D');
                } else if (lastMove == 3) {
                    cigar->push_back('X');
//                    cigar->push_back('M');
                } else if (lastMove == EDLIB_N) {
                    cigar->push_back('N');
                } else if (lastMove == 4) {
                    cigar->push_back('S');
                } else {
                    delete cigar;
                    return EDLIB_STATUS_ERROR;
                }
            }
            if (i < alignmentLength) {
                numOfSameMoves = 0;
                lastMove = alignment_char;
            }
        }
        numOfSameMoves++;
    }
    cigar->push_back(0);  // Null character termination.
    *cigar_ = (char*) malloc(cigar->size() * sizeof(char));
    // memcpy is not used because it was updated recently and requires new GLIBC symbols,
    // thus binaries won't run on older versions of systems.
    memmove(*cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    if (alignment_with_clipping != alignment)
      free(alignment_with_clipping);

    return EDLIB_STATUS_OK;
}

int AlignmentToExtendedCigar(unsigned char* alignment, int alignmentLength, char** cigar_) {
  std::vector<char>* cigar = new std::vector<char>();

  unsigned char *alignment_with_clipping = alignment;

  if (alignment[0] == EDLIB_I || alignment[alignmentLength-1] == EDLIB_I) {
    alignment_with_clipping = (unsigned char *) malloc(sizeof(unsigned char) * (alignmentLength + 1));
    memcpy(alignment_with_clipping, alignment, alignmentLength);
    alignment_with_clipping[alignmentLength] = '\0';

    for (int64_t i=0; i<alignmentLength; i++) {
      if (alignment_with_clipping[i] == EDLIB_I)
        alignment_with_clipping[i] = EDLIB_S;
      else break;
    }
    for (int64_t i=(((int64_t) alignmentLength) - 1); i>=0; i--) {
      if (alignment_with_clipping[i] == EDLIB_I)
        alignment_with_clipping[i] = EDLIB_S;
      else break;
    }
  }

    unsigned char lastMove = -1;  // Code of last move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || alignment_with_clipping[i] != lastMove) {
            if (i > 0) {  // if previous sequence of same moves ended
                // Write number of moves to cigar string.
                int numDigits = 0;
                for (; numOfSameMoves; numOfSameMoves /= 10) {
                    cigar->push_back('0' + numOfSameMoves % 10);
                    numDigits++;
                }
                reverse(cigar->end() - numDigits, cigar->end());
                // Write code of move to cigar string.
                if (lastMove == EDLIB_M || lastMove == EDLIB_EQUAL) {
                    cigar->push_back('=');
//                    cigar->push_back('M');
               } else if (lastMove == EDLIB_I) {
                    cigar->push_back('I');
                } else if (lastMove == EDLIB_D) {
                    cigar->push_back('D');
                } else if (lastMove == EDLIB_X) {
                    cigar->push_back('X');
//                    cigar->push_back('M');
                } else if (lastMove == EDLIB_N) {
                    cigar->push_back('N');
                } else if (lastMove == EDLIB_S) {
                    cigar->push_back('S');
                } else {
                    delete cigar;
                    return EDLIB_STATUS_ERROR;
                }
            }
            if (i < alignmentLength) {
                numOfSameMoves = 0;
                lastMove = alignment_with_clipping[i];
            }
        }
        numOfSameMoves++;
    }
    cigar->push_back(0);  // Null character termination.
    *cigar_ = (char*) malloc(cigar->size() * sizeof(char));
    // memcpy is not used because it was updated recently and requires new GLIBC symbols,
    // thus binaries won't run on older versions of systems.
    memmove(*cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    if (alignment_with_clipping != alignment)
      free(alignment_with_clipping);

    return EDLIB_STATUS_OK;
}

// Converts the alignment array to a more compact CIGAR form, where each CIGAR operation is a std::pair of single-letter operation name and count.
int AlignmentToExtendedCigarArray(unsigned char* alignment, int alignmentLength, std::vector<CigarOp> &cigar) {
  cigar.clear();

  std::vector<unsigned char> alignment_with_clipping;
  alignment_with_clipping.assign(alignment, alignment + alignmentLength);
  alignment_with_clipping.push_back('\0');

  if (alignment[0] == EDLIB_I || alignment[alignmentLength-1] == EDLIB_I) {

    for (int64_t i=0; i<alignmentLength; i++) {
      if (alignment_with_clipping[i] == EDLIB_I)
        alignment_with_clipping[i] = EDLIB_S;
      else break;
    }
    for (int64_t i=(((int64_t) alignmentLength) - 1); i>=0; i--) {
      if (alignment_with_clipping[i] == EDLIB_I)
        alignment_with_clipping[i] = EDLIB_S;
      else break;
    }
  }

  unsigned char lastMove = -1;  // Code of last move.
  int numOfSameMoves = 0;
  for (int i = 0; i <= alignmentLength; i++) {
      // if new sequence of same moves started
      if (i == alignmentLength || alignment_with_clipping[i] != lastMove) {
          if (i > 0) {  // if previous sequence of same moves ended
              // Write number of moves to cigar string.
              char move_op = EdlibOpToCharExtended(lastMove);
              if (move_op == 0) {
                return EDLIB_STATUS_ERROR;
              }
              CigarOp cigar_op(move_op, numOfSameMoves, 0, 0);
//              cigar.push_back(std::make_pair(move_op, numOfSameMoves));
              cigar.push_back(cigar_op);
          }
          if (i < alignmentLength) {
              numOfSameMoves = 0;
              lastMove = alignment_with_clipping[i];
          }
      }
      numOfSameMoves++;
  }

  // Initialize reference and read positions.
  int64_t pos_on_ref = 0;
  int64_t pos_on_read = 0;
  for (int64_t i=0; i<cigar.size(); i++) {
    cigar[i].pos_query = pos_on_read;
    cigar[i].pos_ref = pos_on_ref;

    if (cigar[i].op == 'M' || cigar[i].op == '=' || cigar[i].op == 'X') {
      pos_on_read += cigar[i].count;
      pos_on_ref += cigar[i].count;
    } else if (cigar[i].op == 'I' || cigar[i].op == 'S') {
      pos_on_read += cigar[i].count;
    } else if (cigar[i].op == 'D') {
      pos_on_ref += cigar[i].count;
    }
  }

  return EDLIB_STATUS_OK;
}

int64_t CalculateReconstructedLength(unsigned char *alignment, int alignmentLength) {
  if (alignment == NULL || alignmentLength == 0)
      return 0;

  int64_t length = 0;
  std::stringstream ss;

  for (int i=0; i<alignmentLength; i++) {
    if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL || alignment[i] == EDLIB_X || alignment[i] == EDLIB_D)
      length += 1;
  }

  return length;
}

std::string PrintAlignmentToString(const unsigned char* query, const int queryLength,
                    const unsigned char* target, const int targetLength,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const int modeCode, int row_width) {
    std::stringstream ss;

    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
          if (alignment[i] != EDLIB_I && alignment[i] != EDLIB_S)
                tIdx--;
        }
    }

    for (int start = 0; start < alignmentLength; start += row_width) {
        // target
        ss << "T: ";
        int startTIdx = 0;
        for (int j = start; j < start + row_width && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_I || alignment[j] == EDLIB_S) {
                ss << "_";

            } else {
                ss << (char) target[++tIdx];
            }
            if (j == start)
                startTIdx = tIdx;
        }
        ss << " (" << std::max(startTIdx, 0) << " - " << tIdx << ")\n";

        // Mismatch
        ss << "   ";
        for (int j = start; j < start + row_width && j < alignmentLength; j++) {
            if (alignment[j] == 0)
              ss << "|";
            else if (alignment[j] == EDLIB_X)
              ss << "X";
            else if (alignment[j] == EDLIB_S)
              ss << "-";
            else
              ss << " ";
        }
        ss << "\n";

        // query
        ss << "Q: ";
        int startQIdx = qIdx;
        for (int j = start; j < start + row_width && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_D) {
              ss << "_";
            } else {
                ss << (char) query[++qIdx];
            }
            if (j == start)
                startQIdx = qIdx;
        }
        ss << " (" << std::max(startQIdx, 0) << " - " << qIdx << ")\n\n";
    }

    return ss.str();
}

int CountAlignmentOperations(std::vector<unsigned char>& alignment, const int8_t *read_data, const int8_t *ref_data, int64_t reference_hit_id, int64_t alignment_position_start, SeqOrientation orientation,
                             int64_t match, int64_t mismatch, int64_t gap_open, int64_t gap_extend,
                             bool skip_leading_and_trailing_insertions,
                             int64_t* ret_eq, int64_t* ret_x, int64_t* ret_i, int64_t* ret_d, int64_t *ret_alignment_score, int64_t *ret_edit_dist, int64_t *ret_nonclipped_length) {
  unsigned char last_move = -1;  // Code of last move.
  int64_t num_same_moves = 0;
  int64_t read_position = 0;
  int64_t ref_position = 0;

  int64_t num_eq = 0;
  int64_t num_x = 0;
  int64_t num_i = 0;
  int64_t num_d = 0;
  int64_t alignment_score = 0;

  int64_t nonclipped_length = 0;

  int64_t start_op = 0, end_op = alignment.size() - 1;
  if (skip_leading_and_trailing_insertions == true) {
    for (start_op = 0; start_op < alignment.size(); start_op++, read_position++) {
      if (alignment[start_op] != EDLIB_I && alignment[start_op] != EDLIB_S) { break; }
    }
    for (end_op = (alignment.size() - 1); end_op >= 0; end_op--) {
      if (alignment[end_op] != EDLIB_I && alignment[start_op] != EDLIB_S) { break; }
    }
  }

  for (int i = start_op; i <= end_op; i++) {
    char align_op = 255;
    align_op = alignment[i];

    if (align_op == EDLIB_M || align_op == EDLIB_EQUAL || align_op == EDLIB_X) {
      if (read_data[read_position] == ref_data[alignment_position_start + ref_position]) {
        num_eq += 1;
        alignment_score += match;

      } else {
        num_x += 1;
        alignment_score -= mismatch;
      }

      nonclipped_length += 1;

    } else if (align_op == EDLIB_I) {
      num_i += 1;
      /// This is for the (gap_open + (N) * gap_extend).
      ///      alignment_score -= ((i == 0 || (i > 0 && alignment[i-1] != align_op)) ? (gap_open + gap_extend) : gap_extend);
      /// This is for the (gap_open + (N - 1) * gap_extend).
      alignment_score -= ((i == 0 || (i > 0 && alignment[i-1] != align_op)) ? (gap_open) : gap_extend);
      nonclipped_length += 1;

    } else if (align_op == EDLIB_D) {
      num_d += 1;
      /// This is for the (gap_open + (N) * gap_extend).
      ///      alignment_score -= ((i == 0 || (i > 0 && alignment[i-1] != align_op)) ? (gap_open + gap_extend) : gap_extend);
      /// This is for the (gap_open + (N - 1) * gap_extend).
      alignment_score -= ((i == 0 || (i > 0 && alignment[i-1] != align_op)) ? (gap_open) : gap_extend);
    }

    // Increase coordinates.
    if (align_op == EDLIB_M || align_op == EDLIB_EQUAL || align_op == EDLIB_X || align_op == EDLIB_I || align_op == EDLIB_S)
      read_position += 1;
    if (align_op == EDLIB_M || align_op == EDLIB_EQUAL || align_op == EDLIB_X || align_op == EDLIB_D)
      ref_position += 1;
  }

  *ret_eq = num_eq;
  *ret_x = num_x;
  *ret_i = num_i;
  *ret_d = num_d;
  *ret_alignment_score = alignment_score;
  *ret_nonclipped_length = nonclipped_length;
  *ret_edit_dist = num_x + num_i + num_d;

  return 0;
}

/** MD is a string describing the alignment on the reference from the reference's point of view.
 * It is composed only of match, mismatch and deletion CIGAR operations, other operations are
 * ignored (including clipping).
 * One needs to be careful to add-up matches/mismatches normally separated by I operations, e.g.
 * for a CIGAR "5=2I3=" the MD would be "8".
 */
std::string AlignmentToMD(std::vector<unsigned char>& alignment, const int8_t *read_data, const int8_t *ref_data, int64_t reference_hit_id, int64_t alignment_position_start) {
  std::vector<CigarOp> cigar_array;
  AlignmentToExtendedCigarArray(&alignment[0], alignment.size(), cigar_array);
//  printf("\n");
//  for (int32_t i=0; i<alignment.size(); i++) {
//    printf ("%d", alignment[i]);
//  }
//  printf ("\n");
//  fflush(stdout);

  // This loop truncates the CIGAR operations to skip anything but M and D operations.
  // It keeps track of an offset which points to the last index to where an op was copied.
  // If the current operation is the same as a previous one, then the count are added-up.
  int32_t offset = 0;
  for (int32_t i=0; i<cigar_array.size(); i++) {
    if (cigar_array[i].op == 'I' || cigar_array[i].op == 'S') {     // I and S should not be counted. Truncate them.
      offset += 1;

    } else if (i > 0 && (i-offset-1) > 0 && cigar_array[i].op == cigar_array[i-offset-1].op) {  // Handle same operations
      // The "(i-offset-1) > 0" comes into play when there are leading insertions/soft clips. Prevents invalid reads of size 1.
      cigar_array[i-offset-1].count += cigar_array[i].count;
      offset += 1;

    } else {      // Just truncate the current op.
      cigar_array[i-offset] = cigar_array[i];
    }
  }

  // Offset is the last CIGAR operation that was updated.
  cigar_array.resize(cigar_array.size() - offset);

  std::stringstream md;

  for (int32_t i=0; i<cigar_array.size(); i++) {
    int64_t ref_position = cigar_array[i].pos_ref + alignment_position_start;

    if (cigar_array[i].op == '=') {
      md << cigar_array[i].count;
//      md << "_";
    } else if (cigar_array[i].op == 'X') {
      for (int32_t j=0; j<cigar_array[i].count; j++) {
        md << ref_data[ref_position + j];
        if ((j + 1) < cigar_array[i].count) {
          md << '0';
        }
      }
    } else if (cigar_array[i].op == 'D') {
      md << '^';
      for (int32_t j=0; j<cigar_array[i].count; j++) {
        md << ref_data[ref_position + j];
      }
    }

    if ((i + 1) < cigar_array.size() && cigar_array[i].op != '=' && cigar_array[i+1].op != '=') {
      md << '0';
    }
  }

  return md.str();
}

int GetAlignmentPatterns(const unsigned char* query, const int64_t queryLength,
                         const unsigned char* target, const int64_t targetLength,
                         const unsigned char* alignment, const int64_t alignmentLength,
                         std::string& ret_query, std::string& ret_target, std::string& ret_match_pattern) {
  std::stringstream ss_query, ss_target, ss_pattern;
  int64_t qpos = 0, tpos = 0;
  for (int64_t i=0; i<alignmentLength; i++) {
    if (alignment[i] == EDLIB_H) {
      continue;

    } else if (alignment[i] == EDLIB_S) {
      qpos += 1;
      continue;

    } else if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL || alignment[i] == EDLIB_X) {
      ss_query << query[qpos];
      ss_target << target[tpos];
      ss_pattern << ((alignment[i] == EDLIB_X) ? "*" : "|");
      qpos += 1;
      tpos += 1;

    } else if (alignment[i] == EDLIB_I) {
      ss_query << query[qpos];
      ss_target << "-";
      ss_pattern << "*";
      qpos += 1;

    } else if (alignment[i] == EDLIB_D) {
      ss_query << "-";
      ss_target << target[tpos];
      ss_pattern << "*";
      tpos += 1;

    } else {
      fprintf (stderr, "ERROR: Unknown EDLIB operation %d!\n", alignment[i]);
    }
  }

  ret_query = ss_query.str();
  ret_target = ss_target.str();
  ret_match_pattern = ss_pattern.str();

  return 0;
}

/* Sometimes the aligner (Edlib) can produce deletions right after/before clipping (insertion operations).
 * Normally, a combination of I and D would be an M, but if D's follow I's right at the beginning of the read
 * then they can be removed provided that the start and end reference coordinates of the alignment are moved.
 */
void FixAlignmentLeadingTrailingID(std::vector<unsigned char>& alignment, int64_t *ref_start, int64_t *ref_end) {
  // Find the D stretch at the front, if any.
  int64_t front_start = 0;
  int64_t front_end = 0;
  int32_t front_state = alignment[0];
  for (int64_t i=0; i<alignment.size(); i++) {
    if (front_state == EDLIB_I) {
      if (alignment[i] == EDLIB_I) {
        // Pass.
      } else if (alignment[i] == EDLIB_S) {
        // Pass.
      } else if (alignment[i] == EDLIB_D) {
        front_start = i;
        front_end = i + 1;
        front_state = EDLIB_D;
      } else {
        break;
      }
    } else if (front_state == EDLIB_D) {
      if (alignment[i] == EDLIB_D) {
        front_end = i + 1;
      } else {
        break;
      }
    } else {
      break;
    }
  }

  // Find the D stretch at the front, if any.
  int64_t back_start = alignment.size();
  int64_t back_end = alignment.size();
  int32_t back_state = alignment.back();
  for (int64_t i=(alignment.size()-1); i>=0; i--) {
    if (back_state == EDLIB_I) {
      if (alignment[i] == EDLIB_I) {
        // Pass.
      } else if (alignment[i] == EDLIB_S) {
        // Pass.
      } else if (alignment[i] == EDLIB_D) {
        back_end = i + 1;
        back_start = i;
        back_state = EDLIB_D;
      } else {
        break;
      }
    } else if (back_state == EDLIB_D) {
      if (alignment[i] == EDLIB_D) {
        back_start = i;
      } else {
        break;
      }
    } else {
      break;
    }
  }

  // Otherwise, no changes should be made.
  if (front_end > front_start || back_start < back_end) {
    std::vector<unsigned char> new_alignment;

    if (front_start > 0) {
      new_alignment.insert(new_alignment.end(), alignment.begin(), (alignment.begin()+front_start));
    }
    new_alignment.insert(new_alignment.end(), (alignment.begin() + front_end), (alignment.begin() + back_start));
    if (back_end < alignment.size()) {
      new_alignment.insert(new_alignment.end(), (alignment.begin() + back_end), alignment.end());
    }

    alignment = new_alignment;
    *ref_start += (front_end - front_start);
    *ref_end -= (back_end - back_start);
  }
}
