//============================================================================
// Name        : cmdparser.h
// Author      : Ivan Sovic
// Version     : 1.0
// Created on  : Sep 21, 2014
// Copyright   : License: MIT
// Description : Library for parsing command line parameters.
//============================================================================

#ifndef CMDPARSER_H_
#define CMDPARSER_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>

#define VALUE_TYPE_INT    "INT"
#define VALUE_TYPE_FLOAT  "FLT"
#define VALUE_TYPE_STRING "STR"
#define VALUE_TYPE_NONE   " - "

struct Argument {
  std::string arg_short = "";
  std::string arg_long = "";
  std::string value = "";
  std::string default_value = "";
  std::string value_type = VALUE_TYPE_NONE;
  std::string description_short = "";
  std::string arg_group = "";
  int32_t positional = 0;
  int32_t count = 0;
  bool is_set = false;
};

struct positional_less_than_key {
    inline bool operator() (const Argument& op1, const Argument& op2) {
      return (op1.positional < op2.positional);
    }
};

class ArgumentParser {
 public:
  ArgumentParser();
  ArgumentParser(int argc, char* argv[]);
  ~ArgumentParser();

  void AddArgument(std::string arg_short, std::string arg_long,
                   std::string value_type, int64_t default_value,
                   std::string description_short, int32_t positional=0,
                   std::string argument_group="unknown");
  void AddArgument(std::string arg_short, std::string arg_long,
                   std::string value_type, float default_value,
                   std::string description_short, int32_t positional=0,
                   std::string argument_group="unknown");
  void AddArgument(std::string arg_short, std::string arg_long,
                   std::string value_type, std::string default_value,
                   std::string description_short, int32_t positional=0,
                   std::string argument_group="unknown");
  void ProcessArguments(int argc, char* argv[]);
  std::string VerboseArgumentsByGroup();
  void VerboseArguments(FILE *fp);

  Argument* GetArgument(std::string arg_name);
  Argument* GetArgumentByShortName(std::string arg_name);
  Argument* GetArgumentByLongName(std::string arg_name);

  int GetValue(std::string arg_name, int64_t *value);
  int GetValue(std::string arg_name, float *value);
  int GetValue(std::string arg_name, std::string *value);

 private:
  std::string WrapString_(int32_t leading_tab, int32_t wrap_width, std::string text);

  std::map<std::string, int32_t> valid_args_all_;
  std::map<std::string, int32_t> valid_args_short_;
  std::map<std::string, int32_t> valid_args_long_;
  std::map<std::string, std::vector<int32_t>> valid_args_group_;
  std::map<int32_t, int32_t> valid_args_positional_;
  std::vector<Argument> arguments_;
  std::string program_name_;
};

#endif /* CMDPARSER_H_ */
