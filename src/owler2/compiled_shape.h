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
};

#endif /* SRC_OWLER2_COMPILED_SHAPE_H_ */
