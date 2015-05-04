/*
 * cigargen.h
 *
 *  Created on: Aug 28, 2014
 *      Author: ivan
 */

#ifndef CIGARGEN_H_
#define CIGARGEN_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <sstream>
#include <algorithm>

#include "alignment/myers.h"

#define EDLIB_M 0
#define EDLIB_EQUAL 0
#define EDLIB_X 3
#define EDLIB_I 1
#define EDLIB_D 2
#define EDLIB_S 4

std::string AlignmentToCigar(unsigned char *alignment, int alignmentLength);
int64_t CalculateReconstructedLength(unsigned char *alignment, int alignmentLength);
void PrintAlignment(const unsigned char* query, const int queryLength,
                    const unsigned char* target, const int targetLength,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const int modeCode, const char* idxToLetter);
std::string PrintAlignmentToString(const unsigned char* query, const int queryLength,
                    const unsigned char* target, const int targetLength,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const int modeCode);
std::string AlignmentToCigarReverse(unsigned char *alignment, int alignmentLength);
int edlibAlignmentToCigarForward(unsigned char* alignment, int alignmentLength, char** cigar_);
int edlibAlignmentToCigarForward2(unsigned char* alignment, int alignmentLength, char** cigar_);
int edlibAlignmentToCigarReverse(unsigned char* alignment, int alignmentLength, char** cigar_);
int edlibAlignmentToCigarReverse2(unsigned char* alignment, int alignmentLength, char** cigar_);
int ConvertInsertionsToClipping(unsigned char* alignment, int alignmentLength);

std::string ReverseCigarString(std::string &cigar);

#endif /* CIGARGEN_H_ */
