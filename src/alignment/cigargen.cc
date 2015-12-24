/*
 * cigargen.cc
 *
 *  Created on: Aug 28, 2014
 *      Author: ivan
 */

#include "alignment/cigargen.h"
#include "utility/utility_general.h"



std::string AlignmentToCigar(unsigned char *alignment, int alignmentLength) {
  if (alignment == NULL || alignmentLength == 0)
    return (std::string("*-"));

  std::stringstream ss;

  char* cigar = NULL;
#ifdef USE_EXTENDED_CIGAR_FORMAT
  edlibAlignmentToCigarForward2(alignment, alignmentLength, &cigar);
#else
  edlibAlignmentToCigarForward(alignment, alignmentLength, &cigar);
#endif
  ss << cigar;
  if (cigar)
    free(cigar);

  return ss.str();
}

std::string AlignmentToCigarReverse(unsigned char *alignment, int alignmentLength) {
  if (alignment == NULL || alignmentLength == 0)
    return (std::string("*-"));

  std::stringstream ss;

  char* cigar = NULL;
#ifdef USE_EXTENDED_CIGAR_FORMAT
  edlibAlignmentToCigarReverse2(alignment, alignmentLength, &cigar);
#else
  edlibAlignmentToCigarReverse(alignment, alignmentLength, &cigar);
#endif
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

int edlibAlignmentToCigarForward(unsigned char* alignment, int alignmentLength,
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
                } else if (lastMove == 4) {
                    cigar->push_back('S');
                } else {
                    delete cigar;
                    return MYERS_STATUS_ERROR;
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

    return MYERS_STATUS_OK;
}

int edlibAlignmentToCigarForward2(unsigned char* alignment, int alignmentLength, char** cigar_) {
  std::vector<char>* cigar = new std::vector<char>();
    unsigned char lastMove = -1;  // Code of last move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || alignment[i] != lastMove) {
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
                    cigar->push_back('=');
//                    cigar->push_back('M');
               } else if (lastMove == 1) {
                    cigar->push_back('I');
                } else if (lastMove == 2) {
                    cigar->push_back('D');
                } else if (lastMove == 3) {
                    cigar->push_back('X');
//                    cigar->push_back('M');
                } else if (lastMove == 4) {
                    cigar->push_back('S');
                } else {
                    delete cigar;
                    return MYERS_STATUS_ERROR;
                }
            }
            if (i < alignmentLength) {
                numOfSameMoves = 0;
                lastMove = alignment[i];
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

    return MYERS_STATUS_OK;
}

int edlibAlignmentToCigarReverse(unsigned char* alignment, int alignmentLength,
                          char** cigar_) {
    std::vector<char>* cigar = new std::vector<char>();
    unsigned char lastMove = -1;  // Code of last move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
      char alignment_char = 0;
      if (i < alignmentLength) {
        alignment_char = alignment[alignmentLength - i - 1];
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
                } else if (lastMove == 4) {
                    cigar->push_back('S');
                } else {
                    delete cigar;
                    return MYERS_STATUS_ERROR;
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

    return MYERS_STATUS_OK;
}

int edlibAlignmentToCigarReverse2(unsigned char* alignment, int alignmentLength,
                          char** cigar_) {
  std::vector<char>* cigar = new std::vector<char>();
    unsigned char lastMove = -1;  // Code of last move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || alignment[alignmentLength - i - 1] != lastMove) {
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
                    cigar->push_back('=');
//                    cigar->push_back('M');
               } else if (lastMove == 1) {
                    cigar->push_back('I');
                } else if (lastMove == 2) {
                    cigar->push_back('D');
                } else if (lastMove == 3) {
                    cigar->push_back('X');
//                    cigar->push_back('M');
                } else if (lastMove == 4) {
                    cigar->push_back('S');
                } else {
                    delete cigar;
                    return MYERS_STATUS_ERROR;
                }
            }
            if (i < alignmentLength) {
                numOfSameMoves = 0;
                lastMove = alignment[alignmentLength - i - 1];
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

    return MYERS_STATUS_OK;
}

int64_t CalculateReconstructedLength(unsigned char *alignment, int alignmentLength) {
  if (alignment == NULL || alignmentLength == 0)
      return 0;

  int64_t length = 0;
  std::stringstream ss;

  for (int i=0; i<alignmentLength; i++) {
//    if (alignment[i] != 1)
    if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL || alignment[i] == EDLIB_X || alignment[i] == EDLIB_D)
      length += 1;
  }

  return length;
}

std::string PrintAlignmentToString(const unsigned char* query, const int queryLength,
                    const unsigned char* target, const int targetLength,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const int modeCode) {
    std::stringstream ss;

    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == MYERS_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != 1)
                tIdx--;
        }
    }

    int row_width = 100;

    for (int start = 0; start < alignmentLength; start += row_width) {
        // target
        ss << "T: ";
        int startTIdx = 0;
        for (int j = start; j < start + row_width && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_I || alignment[j] == EDLIB_S)
                ss << "_";
            else
                ss << (char) target[++tIdx];
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
            if (alignment[j] == EDLIB_D)
              ss << "_";
            else
                ss << (char) query[++qIdx];
            if (j == start)
                startQIdx = qIdx;
        }
        ss << " (" << std::max(startQIdx, 0) << " - " << qIdx << ")\n\n";
    }

    return ss.str();
}
