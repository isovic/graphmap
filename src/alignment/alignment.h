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
int SemiglobalAlignment(AlignmentFunctionType AlignmentFunction,
                        const SingleSequence *read, const Index *index, const ProgramParameters *parameters,
                        const EValueParams *evalue_params, PathGraphEntry *region_results);
int AnchoredAlignmentNew(AlignmentFunctionType AlignmentFunctionNW, AlignmentFunctionType AlignmentFunctionSHW,
                         const SingleSequence *read, const Index *index, const ProgramParameters *parameters,
                         const EValueParams *evalue_params, PathGraphEntry *region_results, bool align_end_to_end, bool spliced_alignment);

void VerboseAlignment(const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const AlignmentResults *aln);

/// Determines the start and end locations for semiglobal alignment, keeping in mind the boundaries of the reference being aligned to. Works with circular alignment as well.
//int GetAlignmentWindowFromRegion(const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const PathGraphEntry *region_results,
//                                 int64_t *win_start, int64_t *win_end, int64_t *win_len);
int GetL1PosInRegion(const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const PathGraphEntry *region_results,
                     int64_t *l1_start, int64_t *l1_end);

// Checks if the region is linear or circular. If it's linear, only a pointer to the beginning of the region (in the index) will be returned. Otherwise, a data array will be created containing the
// concatenated region.
// Returns 0 if the region was linear, otherwise 1. Value of 1 means that manual cleanup of ret_data is required, using free().
int GetAlignmentWindowData(const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const PathGraphEntry *region_results,
                           int8_t** data, int64_t* data_length, int8_t **pos_of_win_start, int8_t **pos_of_win_end, int64_t* offset_from_ref_start, int64_t* pos_of_ref_end, bool *is_cleanup_required);

int FindCircularEnd(const std::vector<uint8_t> &alignment, int64_t pos_of_ref_end,
                    int64_t *ret_end_on_aln, int64_t *ret_end_on_read, int64_t *ret_end_on_ref,
                    int64_t *ret_start_on_aln, int64_t *ret_start_on_read, int64_t *ret_start_on_ref);

int SplitCircularAlignment(const AlignmentResults *aln, int64_t pos_of_ref_end, int64_t ref_start, int64_t ref_len, AlignmentResults *aln_l, AlignmentResults *aln_r);

#endif /* SRC_ALIGNMENT_ALIGNMENT_H_ */
