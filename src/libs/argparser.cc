//============================================================================
// Name        : argparser.cc
// Author      : Ivan Sovic
// Version     : 1.0
// Created on  : Sep 21, 2014
// Copyright   : License: MIT
// Description : Library for parsing command line parameters.
//============================================================================

#include "argparser.h"
#include <algorithm>
#include <wordexp.h>

ArgumentParser::ArgumentParser() {
//  AddArgument(NULL, VALUE_TYPE_NONE, "h", "help", "", "Displays this list of commands.", 0, "Help");
}

ArgumentParser::~ArgumentParser() {
  valid_args_all_.clear();
  valid_args_short_.clear();
  valid_args_long_.clear();
  valid_args_group_.clear();
  valid_args_positional_.clear();
  arguments_.clear();
  program_name_.clear();
  composite_args_.clear();
}

void ArgumentParser::AddArgument(void *target,
                                 ValueType value_type,
                                 std::string arg_short, std::string arg_long,
                                 std::string default_value,
                                 std::string description,
                                 int32_t positional,
                                 std::string argument_group) {

  if (arg_short == "" && arg_long == "") {
    fprintf(stderr, "Warning: Argument initialization. Neither a short nor a long name of argument specified. Omitting.\n\n");
    fflush(stderr);
    return;
  }

  /// Check if argument already exists in any of the lists.
  if ((arg_short != "" && valid_args_short_.find(arg_short) != valid_args_short_.end()) ||
      (arg_long != "" && valid_args_long_.find(arg_long) != valid_args_long_.end())) {
    fprintf(stderr, "Warning: Argument '-%s' / '--%s' already defined. Omitting.\n\n", arg_short.c_str(), arg_long.c_str());
    fflush(stderr);
    return;
  }

  Argument arg;
  arg.target = target;
  arg.arg_short = arg_short;
  arg.arg_long = arg_long;
  arg.value = default_value;
  arg.value_type = value_type;
  arg.default_value = default_value;
  arg.description = description;
  arg.positional = positional;
  arg.arg_group = argument_group;
  arg.is_set = false;

  /// In case the argument doesn't require a parameter (but is a switch instead), and the given value is not 0 or 1, set the default value to 0 (False).
  if (arg.value_type == VALUE_TYPE_NONE || arg.value_type == VALUE_TYPE_BOOL) {
    if (arg.default_value != "0" && arg.default_value != "1") {
      arg.value = "0";
      arg.default_value = "0";
    }
  }

  SetArgumentTarget_(target, value_type, default_value);

  arguments_.push_back(arg);
  valid_args_short_[arg.arg_short] = (arguments_.size() - 1);
  valid_args_long_[arg.arg_long] = (arguments_.size() - 1);
  valid_args_all_[arg.arg_short] = (arguments_.size() - 1);
  valid_args_all_[arg.arg_long] = (arguments_.size() - 1);
  if (arg.positional != 0)
    valid_args_positional_[arg.positional] = (arguments_.size() - 1);

  if (valid_args_group_.find(arg.arg_group) == valid_args_group_.end()) {
    std::vector<int32_t> indices;
    indices.push_back((arguments_.size() - 1));
    valid_args_group_[arg.arg_group] = indices;
    arg_groups_in_order_of_appearance_.push_back(arg.arg_group);
  } else {
    valid_args_group_[arg.arg_group].push_back((arguments_.size() - 1));
  }
}

/// Defines composite parameters. A composite parameter can change more than one
/// parameter through a symbolic name. The argument name is literally expanded with the given
/// arg_expansion string containing other valid command line parameters.
void ArgumentParser::AddCompositeArgument(std::string arg_name, std::string arg_expansion) {
  composite_args_[arg_name] = arg_expansion;
}

/// Checks if a given arg_name is in specified arguments.
/// Returns: 0 if it's not specified, 1 if it's a short argument and 2 if it's a long argument.
Argument* ArgumentParser::CheckArguments_(std::string arg_name) {
  if (arg_name.size() < 2) return NULL;
  if (arg_name[0] != '-') return NULL;

  /// The argument has only '-' specified, which means it should be of short type.
  if (arg_name[1] != '-') {
    auto it = valid_args_short_.find(arg_name.substr(1));
    if (it != valid_args_short_.end())
      return &(arguments_[it->second]);
    return NULL;

  } else {
    /// If there is nothing specified after the "--", then something is wrong.
    if (arg_name.size() == 2) return NULL;
    auto it = valid_args_long_.find(arg_name.substr(2));
    if (it != valid_args_long_.end())
      return &(arguments_[it->second]);
    return NULL;
  }

  return NULL;
}

void ArgumentParser::ProcessArguments(int argc, char* argv[], int offset) {
  program_name_ = std::string(argv[0]);

  if ((argc - offset) == 0 && valid_args_positional_.size() > 0) {
    fprintf (stderr, "ERROR: Not all required arguments specified. Exiting.\n\n");
    exit(1);
  }

  for (int32_t i=offset; i<argc; i++) {
    std::string arg = argv[i];
    std::string arg_next = "";

    // Check if for some reason the string is empty.
    if (arg.size() == 0)
      continue;

    // The argument should start with a dash. Otherwise, it can only be positional, or an error.
    if (arg[0] == '-') {
      // If only a dash is given, then something is wrong (or perhaps that was the intention,
      // for instance, to enable pipe processing).
      if (arg.size() == 1) {
        fprintf (stderr, "ERROR: Unspecified parameter '%s'.\n", arg.c_str());
        exit(1);

      } else {
        ///////////////////////////////////////////////////////////////
        /// Find the argument by either the short or the long name.
        //////////////////////////////////////////////////////////////
        int32_t argument_index = -1;
        // If the second character is not '-' then the argument is specified by its short name.
        if (arg[1] != '-') {      // Handle the short argument assignment.
          std::string arg_name = arg.substr(1);

          // If the argument cannot be found, report it.
          if (valid_args_short_.find(arg_name) == valid_args_short_.end()) {
            fprintf (stderr, "ERROR: Unknown parameter '%s'.\n", arg.c_str());
            exit(1);
          } else {
            argument_index = valid_args_short_[arg_name];
          }

        } else if (arg[1] == '-') {                  // Handle the long argument assignment.
          /// Check if there is nothing after '--' in the given argument name.
          if (arg.size() == 2) {
            fprintf (stderr, "ERROR: Unspecified parameter '%s'. Exiting.\n\n", arg.c_str());
            exit(1);
          } else {
            std::string arg_name = arg.substr(2);
            // If the argument cannot be found, report it.
            if (valid_args_long_.find(arg_name) == valid_args_long_.end()) {
              fprintf (stderr, "ERROR: Unknown argument '%s'. Exiting.\n\n", arg.c_str());
              exit(1);
            } else {
              argument_index = valid_args_long_[arg_name];
            }
          }
        }

        ///////////////////////////////////////////////////////////////
        /// Process the argument.
        //////////////////////////////////////////////////////////////
        if (argument_index >= 0) {
          arguments_[argument_index].count += 1;
          arguments_[argument_index].is_set = true;

          if (arguments_[argument_index].value_type != VALUE_TYPE_NONE && arguments_[argument_index].value_type != VALUE_TYPE_BOOL) {      // This argument expects a value.
            if ((i + 1) < argc) {
              arguments_[argument_index].value = argv[i + 1];
              SetArgumentTarget_(arguments_[argument_index].target, arguments_[argument_index].value_type, arguments_[argument_index].value);
              i += 1;
            } else {
              fprintf (stderr, "ERROR: Argument '%s' expects a value. Exiting.\n", arg.c_str());
              exit(1);
            }

            /// Handle the special case when the argument is a composite one. Need to find the composite definition and expand it (recursively).
            if (arguments_[argument_index].value_type == VALUE_TYPE_COMPOSITE) {
              auto it_arg = composite_args_.find(arguments_[argument_index].value);
              if (it_arg == composite_args_.end()) {
                fprintf (stderr, "ERROR: Composite parameter '%s' not defined. Exiting.\n", arguments_[argument_index].value.c_str());
                exit(1);
              }
              wordexp_t p;
              wordexp(it_arg->second.c_str(), &p, 0);
              ProcessArguments((int32_t) p.we_wordc, (char **) p.we_wordv, 0);
              wordfree(&p);
            }

          } else {
            arguments_[argument_index].value = "1";
            if (arguments_[argument_index].value_type == VALUE_TYPE_BOOL) {
              SetArgumentTarget_(arguments_[argument_index].target, arguments_[argument_index].value_type, arguments_[argument_index].value);
            }
          }
        }
       /////////////////////////////////////////////////////////
       /////////////////////////////////////////////////////////
      }
    } else {      // In this case, the argument can either be positional, or misspelled.
      int32_t position_back = i - argc;
      int32_t position_front = i;

      if (valid_args_positional_.find(position_back) != valid_args_positional_.end()) {
        arguments_[valid_args_positional_[position_back]].value = arg;
        arguments_[valid_args_positional_[position_back]].is_set = true;
        SetArgumentTarget_(arguments_[valid_args_positional_[position_back]].target, arguments_[valid_args_positional_[position_back]].value_type, arguments_[valid_args_positional_[position_back]].value);

      } else if (valid_args_positional_.find(position_front) != valid_args_positional_.end()) {
        arguments_[valid_args_positional_[position_front]].value = arg;
        arguments_[valid_args_positional_[position_front]].is_set = true;
        SetArgumentTarget_(arguments_[valid_args_positional_[position_front]].target, arguments_[valid_args_positional_[position_front]].value_type, arguments_[valid_args_positional_[position_front]].value);

      } else {
        fprintf (stderr, "ERROR: Unknown parameter value: '%s'. Exiting.\n\n", arg.c_str());
        exit(1);
      }
    }
  }

  /// Check if there are any positional arguments that were not loaded. These are necessary arguments.
  for (uint32_t i = 0; i < arguments_.size(); i++) {
    if (arguments_[i].positional != 0 && arguments_[i].is_set == false) {
      std::string argument_name = (arguments_[i].arg_long != "") ? arguments_[i].arg_long : arguments_[i].arg_short;
      fprintf (stderr, "ERROR: Not all required arguments specified. Parameter '%s' missing a value. Exiting.\n", argument_name.c_str());
      exit(1);
    }
  }
}

std::string ArgumentParser::VerboseUsage() {
  int32_t max_arg_len = 0;
  for (int32_t i=0; i<arguments_.size(); i++) {
    max_arg_len = std::max(max_arg_len, (int32_t) arguments_[i].arg_long.length());
  }

  std::map<std::string, std::vector<int32_t>>::iterator it;
  const int32_t short_name_starting_char = 4;
  const int32_t value_type_starting_char = short_name_starting_char + 6 + max_arg_len + 1;
  const int32_t description_starting_char = value_type_starting_char + 3 + 3; /// 3 is the length of the type name, and 3 is the spacing between the type name and the description.
  const int32_t wrap_width = MAX_LINE_LENGTH - description_starting_char;

  std::stringstream ret_ss;

  // Printout the usage of the program, if program arguments have actually been processed already.
  // If they haven't been, then program_name is empty (program_name == argv[0]).
  if (program_name_.size() > 0) {
    ret_ss << "Usage:\n";
    ret_ss << "  " << program_name_.c_str();

    // Positional arguments need to be sorted in the ascending order of their position.
    std::vector<Argument> positional_arguments_back;
    std::vector<Argument> positional_arguments_front;
    for (uint32_t i=0; i<arguments_.size(); i++) {
      if (arguments_.at(i).positional < 0) {
        positional_arguments_back.push_back(arguments_.at(i));
      } else if (arguments_.at(i).positional > 0) {
        positional_arguments_front.push_back(arguments_.at(i));
      }
    }
    std::sort(positional_arguments_back.begin(), positional_arguments_back.end(), positional_less_than_key());
    std::sort(positional_arguments_front.begin(), positional_arguments_front.end(), positional_less_than_key());

    for (uint32_t i=0; i<positional_arguments_front.size(); i++) {
      std::string arg_name = positional_arguments_front[i].arg_long;
      if (arg_name == "")
        arg_name = positional_arguments_front[i].arg_short;
      ret_ss << " " << arg_name.c_str();
    }

    ret_ss << " [options]";

    for (uint32_t i=0; i<positional_arguments_back.size(); i++) {
      std::string arg_name = positional_arguments_back[i].arg_long;
      if (arg_name == "")
        arg_name = positional_arguments_back[i].arg_short;
      ret_ss << " " << arg_name.c_str();
    }
    ret_ss << "\n\n";
  }

  ret_ss << "Options\n";

//  for (it = valid_args_group_.begin(); it != valid_args_group_.end(); it++) {
  for (uint32_t group_id = 0; group_id < arg_groups_in_order_of_appearance_.size(); group_id++) {
    it = valid_args_group_.find(arg_groups_in_order_of_appearance_[group_id]);
    if (it == valid_args_group_.end()) { continue; }

    ret_ss << "  " << it->first.c_str() << ":\n";
    for (uint32_t i = 0; i < it->second.size(); i++) {
      std::stringstream ss;
      int32_t num_chars = 0;

//      ss << "    ";
      ss << std::string(short_name_starting_char, ' ');
      num_chars += short_name_starting_char;

      if (arguments_[it->second.at(i)].positional == 0) {
        if (arguments_[it->second.at(i)].arg_short != "") {
          ss << "-" << arguments_[it->second.at(i)].arg_short;
          num_chars += 2;
        } else {
          ss << "  ";
          num_chars += 2;
        }

        if (arguments_[it->second.at(i)].arg_short != "" && arguments_[it->second.at(i)].arg_long != "") {
          ss << ", ";
          num_chars += 2;
        } else {
          ss << "  ";
          num_chars += 2;
        }

        if (arguments_[it->second.at(i)].arg_long != "") {
          ss << "--" << arguments_[it->second.at(i)].arg_long;
          num_chars += 2 + arguments_[it->second.at(i)].arg_long.size();
        }

      } else {
        std::string arg_name = arguments_[it->second.at(i)].arg_long;
        if (arg_name == "")
          arg_name = arguments_[it->second.at(i)].arg_short;

        ss << arg_name;
        num_chars += arg_name.size();
      }

//      if (arguments_[it->second.at(i)].value_type != "") {
      int32_t num_empty_chars = (((value_type_starting_char - num_chars) > 0)?(value_type_starting_char - num_chars):1);
      std::string empty_chars(num_empty_chars, ' ');
      ss << empty_chars;

      ss << ValueTypeToStr(arguments_[it->second.at(i)].value_type);
//      if (arguments_[it->second.at(i)].value_type.size() >= 3)
//        ss << arguments_[it->second.at(i)].value_type.substr(0, 3);
//      else {
//        ss << arguments_[it->second.at(i)].value_type << std::string((3 - arguments_[it->second.at(i)].value_type.size()), ' ');
//      }

      num_chars += num_empty_chars + 3;
//      }

      if (arguments_[it->second.at(i)].description != "") {
        int32_t num_empty_chars = (((description_starting_char - num_chars) > 0)?(description_starting_char - num_chars):1);
        std::string empty_chars(num_empty_chars, ' ');
        ss << empty_chars;
//        printf ("description_starting_char = %d\n", description_starting_char);
//        printf ("wrap_width = %d\n", wrap_width);
//        printf ("description =\n%s\n", arguments_[it->second.at(i)].description.c_str());

        ss << WrapString_(description_starting_char, wrap_width, arguments_[it->second.at(i)].description);
        num_chars += num_empty_chars + arguments_[it->second.at(i)].description.size();
      }
      if (arguments_[it->second.at(i)].value_type != VALUE_TYPE_NONE && arguments_[it->second.at(i)].default_value != "") {
        if (arguments_[it->second.at(i)].value_type != VALUE_TYPE_BOOL) {
          ss << " [" << arguments_[it->second.at(i)].default_value << "]";
        } else {
          ss << " [" << ((arguments_[it->second.at(i)].default_value == "0" || arguments_[it->second.at(i)].default_value == "false") ? "false" : "true") << "]";
        }
      }

      if (arguments_[it->second.at(i)].arg_short != "" || arguments_[it->second.at(i)].arg_long != "") {
        ret_ss << ss.str() << "\n";
      }
    }

    ret_ss << "\n";
  }

  return ret_ss.str();
}

std::string ArgumentParser::VerboseArguments() {
  std::stringstream ss;

  int32_t max_arg_len = 0;
  for (int32_t i=0; i<arguments_.size(); i++) {
    max_arg_len = std::max(max_arg_len, (int32_t) arguments_[i].arg_long.length());
  }

  std::string padding_string(3, ' ');

  for (uint32_t i=0; i<arguments_.size(); i++) {
    ss << "'-" << arguments_[i].arg_short << "'" << std::string(1 + padding_string.size() - arguments_[i].arg_short.size(), ' ');
    ss << "'--" << arguments_[i].arg_long << "'" << std::string(max_arg_len + padding_string.size() - arguments_[i].arg_long.size(), ' ');
    ss << "'" << ValueTypeToStr(arguments_[i].value_type) << "'" << padding_string;
    ss << "value = '" << arguments_[i].value << "'" << padding_string;
    ss << "default = '" << arguments_[i].default_value << "'" << padding_string;
    ss << "count = " << arguments_[i].count << padding_string;
    ss << "positional = " << arguments_[i].positional << "\n";
  }
  return ss.str();
}

Argument* ArgumentParser::GetArgument(std::string arg_name) {
  std::map<std::string, int32_t>::iterator it = valid_args_all_.find(arg_name);
  if (it == valid_args_all_.end())
    return NULL;
  return (&(arguments_.at(it->second)));
}

Argument* ArgumentParser::GetArgumentByShortName(std::string arg_name) {
  std::map<std::string, int32_t>::iterator it = valid_args_short_.find(arg_name);
  if (it == valid_args_short_.end())
    return NULL;
  return (&(arguments_.at(it->second)));
}

Argument* ArgumentParser::GetArgumentByLongName(std::string arg_name) {
  std::map<std::string, int32_t>::iterator it = valid_args_long_.find(arg_name);
  if (it == valid_args_long_.end())
    return NULL;
  return (&(arguments_.at(it->second)));
}

std::string ArgumentParser::WrapString_(int32_t leading_tab, int32_t wrap_width, std::string text) {
  if (text.size() == 0)
    return text;

  std::stringstream ss;
  std::vector<int32_t> spaces;
  std::string prefix(leading_tab, ' ');

  /// Enumerate all whitespaces in the text. This is used to split on word delineations.
  spaces.push_back(-1);
  for (uint32_t i=0; i<text.size(); i++) {
    if (text[i] == ' ' || text[i] == '\t' || text[i] == '\n') {
      spaces.push_back(i);
    }
  }
  spaces.push_back(text.size());

  /// Check for each word whether it contributes with enough characters to exceed the wrap_width.
  int32_t line_length = 0;
  for (uint32_t i=1; i<spaces.size(); i++) {
    int32_t word_length = spaces[i] - (spaces[i - 1] + 1);

    if ((i + 1) == spaces.size() || ((i + 1) < spaces.size() && text[spaces[i-1]] != '\n' && (line_length + word_length) < wrap_width)) {
      if (i > 1) { ss << text[spaces[i-1]]; }
      ss << text.substr((spaces[i - 1] + 1), word_length);
      line_length += word_length + 1;
    } else {
      ss << "\n" << prefix;
      ss << text.substr((spaces[i - 1] + 1), word_length);
      line_length = word_length;
    }
  }

  return ss.str();
}

int ArgumentParser::GetValue(std::string arg_name, int64_t* value) {
  std::map<std::string, int32_t>::iterator it = valid_args_all_.find(arg_name);
  if (it == valid_args_all_.end())
    return 1;
  ValueType value_type = arguments_.at(it->second).value_type;
//  if (value_type != VALUE_TYPE_INT32 || value_type != VALUE_TYPE_INT64 || value_type != VALUE_TYPE_UINT32 || value_type != VALUE_TYPE_UINT64)
  if (value_type != VALUE_TYPE_INT64)
    return 1;
  sscanf(arguments_.at(it->second).value.c_str(), "%ld", value);
  return 0;
}

int ArgumentParser::GetValue(std::string arg_name, float* value) {
  std::map<std::string, int32_t>::iterator it = valid_args_all_.find(arg_name);
  if (it == valid_args_all_.end())
    return 1;
  if (arguments_.at(it->second).value_type != VALUE_TYPE_FLOAT)
    return 1;
  sscanf(arguments_.at(it->second).value.c_str(), "%f", value);
  return 0;
}

int ArgumentParser::GetValue(std::string arg_name, std::string* value) {
  std::map<std::string, int32_t>::iterator it = valid_args_all_.find(arg_name);
  if (it == valid_args_all_.end())
    return 1;
  if (arguments_.at(it->second).value_type != VALUE_TYPE_STRING)
    return 1;
  *value = arguments_.at(it->second).value;
  return 0;
}

void ArgumentParser::SetArgumentTarget_(void *target, ValueType value_type, std::string &argument_value) {
  if (value_type == VALUE_TYPE_NONE)
    return;

  std::stringstream ss(argument_value);
  if (value_type == VALUE_TYPE_BOOL) {
    std::string ss_str = ss.str();
    if (ss_str == "0" || ss_str == "false") {
      *((bool *) target) = false;
    } else {
      *((bool *) target) = true;
    }
  } else if (value_type == VALUE_TYPE_INT32) {
    int32_t value = 0;
    ss >> value;
    *((int32_t *) target) = value;
  } else if (value_type == VALUE_TYPE_UINT32) {
    uint32_t value = 0;
    ss >> value;
    *((uint32_t *) target) = value;
  } else if (value_type == VALUE_TYPE_INT64) {
    int64_t value = 0;
    ss >> value;
    *((int64_t *) target) = value;
  } else if (value_type == VALUE_TYPE_UINT64) {
    uint64_t value = 0;
    ss >> value;
    *((uint64_t *) target) = value;
  } else if (value_type == VALUE_TYPE_FLOAT) {
    float value = 0.0f;
    ss >> value;
    *((float *) target) = value;
  } else if (value_type == VALUE_TYPE_DOUBLE) {
    double value = 0.0f;
    ss >> value;
    *((double *) target) = value;
  } else if (value_type == VALUE_TYPE_STRING) {
    *((std::string *) target) = argument_value;
  }
}
