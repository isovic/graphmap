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
#include "argparser.h"

int main(int argc, char *argv[]) {
  std::string program_name(argv[0]);
  std::string subprogram("");

  ArgumentParser argparser;
  argparser.AddArgument(&subprogram, VALUE_TYPE_STRING, "", "tool", "", "Specifies the tool to run:\n  align - the entire GraphMap pipeline.\n  owler - Overlapping With Long Erroneous Reads.", -1, "");
  argparser.set_program_name(program_name);

  if (argc == 1) {
    fprintf (stderr, "%s", argparser.VerboseUsage().c_str());
    fprintf (stderr, "\n");
    fprintf (stderr, "%s\n", LICENCE_INFORMATION);
    fprintf (stderr, "Version: %s\n", std::string(GRAPHMAP_CURRENT_VERSION).c_str());
    fprintf (stderr, "Build date: %s\n", std::string(GRAPHMAP_CURRENT_VERSION_RELEASE_DATE).c_str());
    fprintf (stderr, "\n");
    exit(1);
  }

  // The ArgumentParser's function for processing arguments is never explicitly called, because it's overly complicated for this purpose.
  // Instead, we just take the value of argv[1] and that's it. ArgumentParser is used only for neat formatting of the usage.
  subprogram = std::string(argv[1]);

  // Remove the 'tools' param to format the command line so it can be seemlesly processed in the next step.
  std::vector<char *> argv2;
  argv2.push_back(argv[0]);
  for (int32_t i=2; i<argc; i++) { argv2.push_back(argv[i]); }
  int32_t argc2 = argv2.size();

  ProgramParameters program_parameters;
  program_parameters.subprogram = subprogram;

  if (subprogram == "align") {
    if (ProcessArgsGraphMap(argc2, &argv2[0], &program_parameters))
      return 1;

    if (program_parameters.verbose_level == 1) {
      LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_STD;
    } else if (program_parameters.verbose_level > 1) {
      LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_FULL | LOG_VERBOSE_STD;
    }
    fflush(stdout);

    GraphMap graphmap;
    graphmap.Run(program_parameters);

  } else if (subprogram == "owler") {
    if (ProcessArgsOwler(argc2, &argv2[0], &program_parameters))
      return 1;

    if (program_parameters.verbose_level == 1) {
      LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_STD;
    } else if (program_parameters.verbose_level > 1) {
      LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_FULL | LOG_VERBOSE_STD;
    }
    fflush(stdout);

    Owler owler;
    owler.Run(program_parameters);

  } else {
    fprintf (stderr, "ERROR: Unknown value of 'tool' parameter. Exiting.\n\n");
    fprintf (stderr, "%s\n", argparser.VerboseUsage().c_str());
    exit(1);

  }


	return 0;
}
