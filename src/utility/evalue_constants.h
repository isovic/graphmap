/*
swsharp - CUDA parallelized Smith Waterman with applying Hirschberg's and 
Ukkonen's algorithm and dynamic cell pruning.
Copyright (C) 2013 Matija Korpar, contributor Mile Šikić

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact the author by mkorpar@gmail.com.
*/
/**
@file

@brief Project constants definitions header.
*/

#ifndef __SW_SHARP_CONSTANTSH__
#define __SW_SHARP_CONSTANTSH__

#ifdef __cplusplus 
extern "C" {
#endif

/*!
@brief Constant used for defining the semiglobal alignment.
*/
#define HW_ALIGN    0

/*!
@brief Constant used for defining the global alignment.
*/
#define NW_ALIGN    1

/*!
@brief Constant used for defining the local alignment.
*/
#define SW_ALIGN    2

/*!
@brief Constant used for defining the overlap alignment.
*/
#define OV_ALIGN    3

/*!
@brief Constant used for defining an unknown score.
*/
#define NO_SCORE    -1000000001

/*!
@brief Constant used for defining score minimum, in other words -infinity.
*/
#define SCORE_MIN   -1000000000

//!@{
/*!
@brief Similarity matrix table scores.
*/
//extern int BLOSUM_45_TABLE[26 * 26];
//extern int BLOSUM_50_TABLE[26 * 26];
//extern int BLOSUM_62_TABLE[26 * 26];
//extern int BLOSUM_80_TABLE[26 * 26];
//extern int BLOSUM_90_TABLE[26 * 26];
//
//extern int PAM_30_TABLE[26 * 26];
//extern int PAM_70_TABLE[26 * 26];
//extern int PAM_250_TABLE[26 * 26];

extern int EDNA_FULL_TABLE_5_4[26 * 26];
extern int EDNA_FULL_TABLE_1_1[26 * 26];
extern int EDNA_FULL_TABLE_4_5[26 * 26];
//!@}

#ifdef __cplusplus 
}
#endif
#endif // __SW_SHARP_CONSTANTSH__
