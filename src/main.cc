//============================================================================
// Name        : graphmap.cpp
// Author      : Ivan Sovic
// Version     :
// Copyright   : Copyright Ivan Sovic, 2014. All rights reserved.
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "sequences/sequence_file.h"
#include "sequences/single_sequence.h"
#include "log_system/log_system.h"
#include "graphmap/graphmap.h"

#include "index/index_sa.h"

#include "program_parameters.h"
#include "utility/utility_general.h"

#include "owler/owler.h"

int main(int argc, char *argv[]) {

  ProgramParameters program_parameters;
  if (ProcessArgs(argc, argv, &program_parameters))
    return 1;

  if (program_parameters.verbose_level == 1) {
    LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_STD;
  } else if (program_parameters.verbose_level > 1) {
    LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_FULL | LOG_VERBOSE_STD;
  }
  fflush(stdout);

  if (program_parameters.alignment_approach == "owler") {
    Owler owler;
    owler.Run(program_parameters);
  } else {
    GraphMap graphmap;
    graphmap.Run(program_parameters);
  }

	return 0;
}
