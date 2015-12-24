//============================================================================
// Name        : argparser.h
// Author      : Ivan Sovic
// Version     : 1.0
// Created on  : Sep 21, 2014
// Copyright   : License: MIT
// Description : Library for parsing command line parameters.
//============================================================================

#ifndef ARGPARSER_H_
#define ARGPARSER_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#define MAX_LINE_LENGTH 120

typedef enum {
  VALUE_TYPE_NONE = 0,
  VALUE_TYPE_BOOL = 1,
  VALUE_TYPE_INT32 = 2,
  VALUE_TYPE_UINT32 = 3,
  VALUE_TYPE_INT64 = 4,
  VALUE_TYPE_UINT64 = 5,
  VALUE_TYPE_FLOAT = 6,
  VALUE_TYPE_DOUBLE = 7,
  VALUE_TYPE_STRING = 8,
  VALUE_TYPE_COMPOSITE = 9
} ValueType;

inline std::string ValueTypeToStr(ValueType value_type) {
  if (value_type == VALUE_TYPE_NONE)
    return std::string(" - ");
  else if (value_type == VALUE_TYPE_BOOL)
    return std::string(" - ");
  else if (value_type == VALUE_TYPE_INT32 || value_type == VALUE_TYPE_UINT32 || value_type == VALUE_TYPE_INT64 || value_type == VALUE_TYPE_UINT64)
    return std::string("INT");
  else if (value_type == VALUE_TYPE_FLOAT || value_type == VALUE_TYPE_DOUBLE)
    return std::string("FLT");
  else if (value_type == VALUE_TYPE_STRING || value_type == VALUE_TYPE_COMPOSITE)
    return std::string("STR");
  return std::string(" - ");
}

struct Argument {
  void *target = NULL;                      /// Pointer to the variable to where the parsed value will be stored.
  std::string arg_short = "";               /// The short name of the argument.
  std::string arg_long = "";                /// The long name of the argument.
  std::string value = "";                   /// The value parsed from the command line.
  std::string default_value = "";           /// Fallback value in case the argument was not specified.
  ValueType value_type = VALUE_TYPE_NONE;   /// The type of the argument (integer, float, string).
  std::string description = "";       /// Description of the argument.
  std::string arg_group = "";               /// The group of the argument.
  int32_t positional = 0;                   /// If 0, the parameter can be placed in an optional position. If > 0 parameter must be located at the front of the command line. If < 0, the parameter must be located at the end of command line.
  int32_t count = 0;                        /// Counts the number of occurances of the arguments.
  bool is_set = false;                      /// If the argument is found, is_set == true.
};

struct positional_less_than_key {
    inline bool operator() (const Argument& op1, const Argument& op2) {
      return (op1.positional < op2.positional);
    }
};

class ArgumentParser {
 public:
  ArgumentParser();
  ~ArgumentParser();

  /// Specify an argument to parse for.
  /// Parameters:
  /// @target Pointer to a variable where the value will be parsed to.
  /// @value_type Specifies the variable type which is expected.
  /// @arg_short Short argument name, single letter. Can be omitted, but either a short or long name needs to be specified.
  /// @arg_long Long argument name, e.g. a word. Can be omitted, but either a short or long name needs to be specified.
  /// @default value The default value of the parameter in case it was never used on the command line.
  /// @description The description that will be used to generate program usage.
  /// @positional If 0, the argument is optional and can be placed anywhere in between required arguments. If > 0 specifies a required argument at this position from the program name. If < 0, the argument is required at the end of the command line.
  /// @argument group Specifies the group of an argument. Used to group arguments when printing usage.
  void AddArgument(void *target,
                   ValueType value_type,
                   std::string arg_short, std::string arg_long,
                   std::string default_value,
                   std::string description, int32_t positional=0,
                   std::string argument_group="unknown");

  /// Defines composite parameters. A composite parameter can change more than one
  /// parameter through a symbolic name. The argument name is literally expanded with the given
  /// arg_expansion string containing other valid command line parameters.
  void AddCompositeArgument(std::string arg_name, std::string arg_expansion);

  /// Processes the command line and parses arguments. Prints usage if argument parameter is invalid.
  /// If offset > 0, then the given number of arguments from the beginning will be skipped (i.e. program name).
  void ProcessArguments(int argc, char* argv[], int offset=1);

  /// Formats and returns the usage of the program according to the specified arguments.
  std::string VerboseUsage();

  /// Formats and returns the values and specifications of all arguments.
  std::string VerboseArguments();

  /// Finds and returns a pointer to the argument by its name. The name can be either a short or a long name.
  Argument* GetArgument(std::string arg_name);

  /// Finds and returns a pointer to the argument by its short name.
  Argument* GetArgumentByShortName(std::string arg_name);

  /// Finds and returns a pointer to the argument by its long name.
  Argument* GetArgumentByLongName(std::string arg_name);

  /// Finds an argument given its name, and parses the integer value.
  int GetValue(std::string arg_name, int64_t *value);

  /// Finds an argument given its name, and parses the float value.
  int GetValue(std::string arg_name, float *value);

  /// Finds an argument given its name, and parses the string value.
  int GetValue(std::string arg_name, std::string *value);

 private:
  /// Wraps a string to a given number of characters in width. Also, assigns a number of tabs at the beginning of each line.
  std::string WrapString_(int32_t leading_tab, int32_t wrap_width, std::string text);

  /// For a given argument value, parses it's value according to the argument_value parameter and stores it into the target variable.
  void SetArgumentTarget_(void *target, ValueType value_type, std::string &argument_value);
  /// Checks if a given arg_name can be found in specified arguments, either long or short ones.
  /// Returns: 0 if it's not specified, 1 if it's a short argument and 2 if it's a long argument.
  /// The value of arg_name must start with a '-', e.g. "-t" or "--threads".
  Argument* CheckArguments_(std::string arg_name);

  std::map<std::string, int32_t> valid_args_all_;
  std::map<std::string, int32_t> valid_args_short_;
  std::map<std::string, int32_t> valid_args_long_;
  std::map<std::string, std::vector<int32_t>> valid_args_group_;
  std::map<int32_t, int32_t> valid_args_positional_;
  std::vector<std::string> arg_groups_in_order_of_appearance_;
  std::vector<Argument> arguments_;
  std::string program_name_;
  std::map<std::string, std::string> composite_args_;
};

#endif /* CMDPARSER_H_ */
