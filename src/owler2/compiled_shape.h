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

#endif /* SRC_OWLER2_COMPILED_SHAPE_H_ */
