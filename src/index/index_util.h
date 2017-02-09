/*
 * index_util.h
 *
 *  Created on: Feb 6, 2017
 *      Author: isovic
 */

#ifndef SRC_INDEX_INDEX_UTIL_H_
#define SRC_INDEX_INDEX_UTIL_H_

#include "minimizer_index/minimizer_index.h"
#include "graphmap/transcriptome.h"
#include "../program_parameters.h"

namespace is {

std::string GenerateSAMHeader(std::shared_ptr<is::MinimizerIndex> index, ProgramParameters &parameters);

std::string GenerateSAMHeader(std::shared_ptr<is::Transcriptome> transcriptome);

} /* namespace is */

#endif /* SRC_INDEX_INDEX_UTIL_H_ */
