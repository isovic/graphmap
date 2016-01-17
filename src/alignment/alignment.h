/*
 * alignment.h
 *
 *  Created on: Jan 17, 2016
 *      Author: isovic
 */

#ifndef SRC_ALIGNMENT_ALIGNMENT_H_
#define SRC_ALIGNMENT_ALIGNMENT_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <sstream>
#include <vector>

#include "sequences/single_sequence.h"
#include "utility/utility_general.h"
#include "utility/program_parameters.h"
#include "utility/utility_conversion-inl.h"
#include "containers/path_graph_entry.h"
#include "libs/myers.h"
#include "alignment/cigargen.h"
#include "log_system/log_system.h"
#include "containers/region.h"
#include "seqan/basic.h"
#include "seqan/align.h"
#include "seqan/sequence.h"
#include "seqan/stream.h"

#include "alignment/local_realignment_wrappers.h"



typedef int (*AlignmentFunctionType)(const int8_t*, int64_t, const int8_t*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t*, int64_t*, int64_t*, std::vector<unsigned char> &);
typedef int (*EditDistanceFunctionType)(const int8_t*, int64_t, const int8_t*, int64_t, int64_t*, int64_t*, int);

//typedef int (*AlignmentFunctionTypeMex)(const int8_t*, int64_t, const int8_t*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t*, int64_t*, int64_t*, std::vector<unsigned char> &);




#endif /* SRC_ALIGNMENT_ALIGNMENT_H_ */
