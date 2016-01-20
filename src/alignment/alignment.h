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
#include "utility/evalue.h"

#include "alignment/local_realignment_wrappers.h"
#include "alignment/local_realignment.h"



int AlignRegion(const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const EValueParams *evalue_params, bool extend_to_end, PathGraphEntry *region_results);
int SemiglobalAlignment(AlignmentFunctionType AlignmentFunction, const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const EValueParams *evalue_params, PathGraphEntry *region_results);
void VerboseAlignment(const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const AlignmentResults *aln);



#endif /* SRC_ALIGNMENT_ALIGNMENT_H_ */
