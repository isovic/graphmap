/*
 * compiled_shape.cc
 *
 *  Created on: Jun 24, 2016
 *      Author: isovic
 */

#include "compiled_shape.h"
#include <math.h>

void CompiledShape::Compile(const std::string new_shape) {
  shape = new_shape;
  masks.clear();
  num_incl_bits = 0;

  Mask mask;
  for (int32_t i=0; i<new_shape.size(); i++) {
//      mask.start = (new_shape[i] == "1" && (i == 0 || (i > 0 && new_shape[i-1] == "0"))) ? i : mask.start;

    if (shape[i] == '1' && (i == 0 || (i > 0 && shape[i-1] == '0'))) {  /// This detects the start of a new chain of inclusive bases.
      mask.start = i*2;
      mask.len = 2; /// Counts two bits per base.
    } else if (shape[i] == '0' && (i > 0 && shape[i-1] == '1')) {       /// This detects the end of the chain of inclusive bases.
      masks.push_back(mask);
      mask.start = 0;
      mask.len = 0;
    } else if (shape[i] == '1') { /// Count the number of inclusive bases.
      mask.len += 2; /// Counts two bits per base.
    }
  }

  /// Leftover mask part.
  if (mask.len > 0) {
    masks.push_back(mask);
    mask.start = 0;
    mask.len = 0;
  }

//    shape = "11110111101111";
//    Mask mask;
//    mask.start = 0; mask.len = 4*2; masks.push_back(mask);
//    mask.start = 5*2; mask.len = 4*2; masks.push_back(mask);
//    mask.start = 10*2; mask.len = 4*2; masks.push_back(mask);

  for (int32_t i=0; i<masks.size(); i++) {
    int64_t shift_mask_raw = 64 - masks[i].start - masks[i].len;
    masks[i].bits = ((uint64_t) (pow(2, masks[i].len) - 1)) << (shift_mask_raw);
    masks[i].shift = (i > 0) ? (masks[i].start - (masks[i-1].start + masks[i-1].len) + masks[i-1].shift) : 0;
    num_incl_bits += masks[i].len;
//      printf ("%s\t\t[%d] start = %d, len = %d, shift = %d, shape_incl_bases = %d\n", ConvertToBinary(masks[i].bits).c_str(), i, masks[i].start, masks[i].len, masks[i].shift, num_incl_bases);
  }
}

uint64_t CompiledShape::CreateSeedFromShape(uint64_t bases2bit) const {
  uint64_t seed = 0;
  for (int32_t i=0; i<masks.size(); i++) {
    uint64_t buffer_part = bases2bit & masks[i].bits;
    buffer_part = buffer_part << masks[i].shift;
    seed |= buffer_part;
  }
  seed = seed >> (64 - num_incl_bits);
  return seed;
}

std::vector<CompiledShape> CompileShapes(const std::vector<std::string> &shapes) {
  std::vector<CompiledShape> compiled_shapes;
  for (int32_t i=0; i<shapes.size(); i++) {
    CompiledShape compiled_shape(shapes[i]);
    compiled_shapes.push_back(compiled_shape);
  }
  return compiled_shapes;
}

