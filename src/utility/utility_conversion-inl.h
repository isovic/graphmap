/*
 * Copyright 2014, Ivan Sovic.
 * All rights reserved.
 *
 * utility_conversion-inc.h
 *
 *  Created on: 15 May, 2014
 *      Author: Ivan Sovic
 */

#ifndef UTILITY_CONVERSION_INC_H_
#define UTILITY_CONVERSION_INC_H_

#include <stdint.h>
#include <math.h>

#define MEMORY_UNIT_BYTE 0
#define MEMORY_UNIT_KILOBYTE 1
#define MEMORY_UNIT_MEGABYTE 2
#define MEMORY_UNIT_GIGABYTE 3

// Each nucleotide (base) is represented with 2 bits, to maximally use the entropy of the data. Four nucleotides can be packed within a single byte
#define BASE_A  0x00
#define BASE_C  0x01
#define BASE_T  0x02
#define BASE_G  0x03

#define ComplementBase(base)    (base^0x02)  // Binary codes for bases are chosen in such a way that by changing the second LSB bit of the code we obtain the bases complement.
#define ComplementByte(bases)   (base^0xAA)  // In case the entire byte (containing 4 bases) needs to be complemented, complementByte(bases) can be used. In this case, every other bit is changed.

// #define ReverseComplement(bases)  ((((bases & 0x03)<<6) | ((bases & 0x0C)<<2) | ((bases & 0x30)>>2) | ((bases & 0xC0)>>6)) ^ 0xAA)

inline char ComplementBaseAscii(char base) {
  if    (base == 'A')  return 'T';
  else if (base == 'C')  return 'G';
  else if (base == 'T')  return 'A';
  else if (base == 'G')  return 'C';

  return 'N';
}

inline unsigned char BaseToChar(char base)
{
  if    (base == BASE_A)  return 'A';
  else if (base == BASE_C)  return 'C';
  else if (base == BASE_T)  return 'T';
  else if (base == BASE_G)  return 'G';

  return 'N';
}

inline unsigned char CharToBase(char charBase)
{
  if    (charBase == 'A')   return BASE_A;
  else if (charBase == 'C') return BASE_C;
  else if (charBase == 'T') return BASE_T;
  else if (charBase == 'G') return BASE_G;

  return 0x00;
}

inline unsigned char CharToBaseBWA(char charBase)
{
  if    (charBase == 'A')   return 0;
  else if (charBase == 'C') return 1;
  else if (charBase == 'T') return 3;
  else if (charBase == 'G') return 2;

  return 0x00;
}

inline unsigned char PackByte(char bases[4])
{
  return ((CharToBase(bases[0]) << 0) |
      (CharToBase(bases[1]) << 2) |
      (CharToBase(bases[2]) << 4) |
      (CharToBase(bases[3]) << 6));
}

inline unsigned char PackByteBWA(char bases[4])
{
  return ((CharToBaseBWA(bases[3]) << 0) |
      (CharToBaseBWA(bases[2]) << 2) |
      (CharToBaseBWA(bases[1]) << 4) |
      (CharToBaseBWA(bases[0]) << 6));
}

inline unsigned char BaseToCharBWA(char base)
{
  if    (base == 0) return 'A';
  else if (base == 1) return 'C';
  else if (base == 3) return 'T';
  else if (base == 2) return 'G';

  return 'N';
}



inline uint8_t ComplementBaseBWA(uint8_t base)
{
  if (base < 4)
    return (3 - base);

  return 4;
}

inline unsigned long long int CalcNumSeedsPerSeq(uint64_t seq_length, uint64_t seed_length, uint64_t seed_step)
{
  if (seed_step==0 || seed_length==0 || seq_length==0)
    return 0;

  return (ceil(((float) (seq_length - seed_length + 1)) / seed_step));
}

inline int64_t FastAbs(int64_t a)
{
  int64_t mask = (a >> (sizeof(uint64_t)*8 - 1));
  return (a + mask) ^ mask;
}

const int8_t kBwaToBase[256] = {
  //A, C,  G,  T,    N
  65, 67, 71, 84,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 0 - 15
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 16 - 31
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 32 - 47
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 48 - 63
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 64 - 79
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 80 - 95
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 96 - 111
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 112 - 127
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 128 - 143
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 144 - 159
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 160 - 176
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 176 - 191
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 192 - 208
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 208 - 223
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 224 - 239
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78   // 240 - 256
};

const int8_t kBaseToBwa[256] = {
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 0 - 15
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 16 - 31
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 32 - 47
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 48 - 63
  4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4,  // 64 - 79 (A, C, G)
  4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 80 - 95 (T)
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 96 - 111
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 112 - 127
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 128 - 143
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 144 - 159
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 160 - 176
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 176 - 191
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 192 - 208
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 208 - 223
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 224 - 239
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4  // 240 - 256
};

const uint8_t kBaseToBwaUnsigned[256] = {
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 0 - 15
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 16 - 31
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 32 - 47
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 48 - 63
  4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4,  // 64 - 79 (A, C, G)
  4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 80 - 95 (T)
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 96 - 111
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 112 - 127
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 128 - 143
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 144 - 159
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 160 - 176
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 176 - 191
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 192 - 208
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 208 - 223
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,  // 224 - 239
  4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4  // 240 - 256
};

const int8_t kIsBase[256] = {
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 0 - 15
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 16 - 31
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 32 - 47
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 48 - 63
  0, 1, 0, 1,   0, 0, 0, 1,   0, 0, 0, 0,   0, 0, 0, 0,  // 64 - 79 (A, C, G)
  0, 0, 0, 0,   1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 80 - 95 (T)
  0, 1, 0, 1,   0, 0, 0, 1,   0, 0, 0, 0,   0, 0, 0, 0,  // 96 - 111
  0, 0, 0, 0,   1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 112 - 127
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 128 - 143
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 144 - 159
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 160 - 176
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 176 - 191
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 192 - 208
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 208 - 223
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,  // 224 - 239
  0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0  // 240 - 256
};

const int8_t kBaseComplement[256] = {
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 0 - 15
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 16 - 31
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 32 - 47
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 48 - 63
  78, 84, 78, 71,   78, 78, 78, 67,   78, 78, 78, 78,   78, 78, 78, 78,  // 64 - 79 (A, C, G)
  78, 78, 78, 78,   65, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 80 - 95 (T)
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 96 - 111
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 112 - 127
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 128 - 143
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 144 - 159
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 160 - 176
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 176 - 191
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 192 - 208
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 208 - 223
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,  // 224 - 239
  78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78,   78, 78, 78, 78  // 240 - 256
};

inline uint64_t ConvertFromBytes (uint32_t to_memory_unit, uint64_t original_value) {
  uint64_t converted_value = 0;

  if (to_memory_unit == MEMORY_UNIT_BYTE)
    converted_value = original_value / ((uint64_t) 1);
  else if (to_memory_unit == MEMORY_UNIT_KILOBYTE)
    converted_value = original_value / ((uint64_t) 1024);
  else if (to_memory_unit == MEMORY_UNIT_MEGABYTE)
    converted_value = original_value / (((uint64_t) 1024) * ((uint64_t) 1024));
  else if (to_memory_unit == MEMORY_UNIT_GIGABYTE)
    converted_value = original_value / (((uint64_t) 1024) * ((uint64_t) 1024) * ((uint64_t) 1024));

  return converted_value;
}

inline uint64_t ConvertToBytes (uint32_t from_memory_unit, uint64_t original_value) {
  uint64_t converted_value = 0;

  if (from_memory_unit == MEMORY_UNIT_BYTE)
    converted_value = original_value * ((uint64_t) 1);
  else if (from_memory_unit == MEMORY_UNIT_KILOBYTE)
    converted_value = original_value * ((uint64_t) 1024);
  else if (from_memory_unit == MEMORY_UNIT_MEGABYTE)
    converted_value = original_value * (((uint64_t) 1024) * ((uint64_t) 1024));
  else if (from_memory_unit == MEMORY_UNIT_GIGABYTE)
    converted_value = original_value * (((uint64_t) 1024) * ((uint64_t) 1024) * ((uint64_t) 1024));

  return converted_value;
}


#endif /* UTILITY_CONVERSION_INC_H_ */
