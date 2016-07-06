/*
 * compiled_shape.h
 *
 *  Created on: Jun 24, 2016
 *      Author: isovic
 */

#ifndef SRC_OWLER2_COMPILED_SHAPE_H_
#define SRC_OWLER2_COMPILED_SHAPE_H_

#include <stdint.h>
#include <string>
#include <vector>

struct Mask {
  int32_t start = 0;
  int32_t len = 0;
  uint64_t bits = 0;
  int32_t shift = 0;
};

class CompiledShape {
 public:
  std::vector<Mask> masks;
  std::string shape = "";
  int32_t num_incl_bits = 0;

  CompiledShape() { };
  CompiledShape(const std::string new_shape) { Compile(new_shape); };
  void Compile(const std::string new_shape);

  /// Takes a buffer of bases (max 32 bases in 64-bits), 2bit packed, and extracts those inclusive ones (defined by a shape).
  /// Parameters:
  ///  @buffer an integer containing 2-bit packed values of the sequence, max. 32 bases from the starting position.
  ///  @shape a string specifying the shape to be extracted from the buffer. Specified with '1' as inclusive bases and '0' as don't care bases.
  uint64_t CreateSeedFromShape(uint64_t bases2bit) const;
};

std::vector<CompiledShape> CompileShapes(const std::vector<std::string> &shapes);

// For a given shape, for every don't cate ('0') base generate all three combinations for the same seed, containing at the position: 0 (deletion), 1 (match/mismatch) and 2 (insertion).
// E.g. for a shape "1110111", this function will generate three shapes: "11111", "1110111" and "11100111".
// The shapes vector is not cleared, so the function can be re-used to add multiple shapes.
int CreateLookupShapes(std::string index_shape, std::vector<std::string> &shapes);

// Converts a base 10 number x to a base N number with n digits.
void Base10ToBaseN(int32_t x, int32_t N, int32_t n, std::vector<int8_t> &digits);

#endif /* SRC_OWLER2_COMPILED_SHAPE_H_ */
