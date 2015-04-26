//============================================================================
// Name        : cmdparser.cc
// Author      : Ivan Sovic
// Version     : 1.0
// Created on  : Sep 21, 2014
// Copyright   : License: MIT
// Description : Library for parsing command line parameters.
//============================================================================

#include "argumentparser.h"

ArgumentParser::ArgumentParser() {
}

ArgumentParser::ArgumentParser(int argc, char* argv[]) {
  ProcessArguments(argc, argv);
}

ArgumentParser::~ArgumentParser() {
  valid_args_all_.clear();
  valid_args_short_.clear();
  valid_args_long_.clear();
  valid_args_group_.clear();
  valid_args_positional_.clear();
  arguments_.clear();
  program_name_.clear();
}

void ArgumentParser::AddArgument(std::string arg_short, std::string arg_long,
                                 std::string value_type, int64_t default_value,
                                 std::string description_short,
                                 int32_t positional,
                                 std::string argument_group) {
  std::stringstream ss;
  ss << default_value;
  AddArgument(arg_short, arg_long, value_type, ss.str(), description_short, positional, argument_group);
}

void ArgumentParser::AddArgument(std::string arg_short, std::string arg_long,
                                 std::string value_type, float default_value,
                                 std::string description_short,
                                 int32_t positional,
                                 std::string argument_group) {
  std::stringstream ss;
  ss << default_value;
  AddArgument(arg_short, arg_long, value_type, ss.str(), description_short, positional, argument_group);
}

void ArgumentParser::AddArgument(std::string arg_short, std::string arg_long,
                                 std::string value_type, std::string default_value,
                                 std::string description_short,
                                 int32_t positional,
                                 std::string argument_group) {

  if (arg_short == "" && arg_long == "") {
    fprintf(stderr, "ERROR: Argument's name not initialized!\n");
    fflush(stderr);
    return;
  }

  Argument arg;
  arg.arg_short = arg_short;
  arg.arg_long = arg_long;
  arg.value = default_value;
  arg.value_type = value_type;
  arg.default_value = default_value;
  arg.description_short = description_short;
  arg.positional = positional;
  arg.arg_group = argument_group;
  arg.is_set = false;

  if (arg.value_type == VALUE_TYPE_NONE) {
    if (arg.default_value != "0" && arg.default_value != "1") {
      arg.value = "0";
      arg.default_value = "0";
    }
  }

  arguments_.push_back(arg);
  valid_args_short_[arg.arg_short] = (arguments_.size() - 1);
  valid_args_long_[arg.arg_long] = (arguments_.size() - 1);
  valid_args_all_[arg.arg_short] = (arguments_.size() - 1);
  valid_args_all_[arg.arg_long] = (arguments_.size() - 1);
  valid_args_positional_[arg.positional] = (arguments_.size() - 1);

  if (valid_args_group_.find(arg.arg_group) == valid_args_group_.end()) {
    std::vector<int32_t> indices;
    indices.push_back((arguments_.size() - 1));
    valid_args_group_[arg.arg_group] = indices;
  } else {
    valid_args_group_[arg.arg_group].push_back((arguments_.size() - 1));
  }
}

void ArgumentParser::ProcessArguments(int argc, char* argv[]) {
  program_name_ = std::string(argv[0]);

  for (int32_t i=1; i<argc; i++) {
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
        fprintf (stderr, "WARNING: Unspecified parameter '%s'.\n", arg.c_str());

      } else {
        // If the second character is not '-' then the argument is specified by its short name.
        if (arg[1] != '-') {      // Handle the short argument assignment.
          std::string arg_name = arg.substr(1);

          // If the argument cannot be found, report it.
          if (valid_args_short_.find(arg_name) == valid_args_short_.end()) {
            fprintf (stderr, "WARNING: Unknown parameter '%s'.\n", arg.c_str());

          } else {
            int32_t argument_index = valid_args_short_[arg_name];
            arguments_[argument_index].count += 1;
            arguments_[argument_index].is_set = true;

            if (arguments_[argument_index].value_type != VALUE_TYPE_NONE) {      // This argument expects a value.
              if ((i + 1) < argc) {
                arguments_[argument_index].value = argv[i + 1];
                i += 1;
              } else {
                fprintf (stderr, "ERROR: Argument '%s' expects a value. Exiting.\n", arg.c_str());
                exit(1);
              }

            } else {
              arguments_[argument_index].value = "1";

            }
          }

        } else if (arg[1] == '-') {                  // Handle the long argument assignment.
          if (arg.size() == 2) {
            fprintf (stderr, "WARNING: Unspecified parameter '%s'.\n", arg.c_str());
          } else {
            std::string arg_name = arg.substr(2);

            // If the argument cannot be found, report it.
            if (valid_args_long_.find(arg_name) == valid_args_long_.end()) {
              fprintf (stderr, "WARNING: Unknown parameter '%s'.\n", arg.c_str());

            } else {
              int32_t argument_index = valid_args_long_[arg_name];
              arguments_[argument_index].count += 1;
              arguments_[argument_index].is_set = true;

              if (arguments_[argument_index].value_type != VALUE_TYPE_NONE) {      // This argument expects a value.
                if ((i + 1) < argc) {
                  arguments_[argument_index].value = argv[i + 1];
                  i += 1;
                } else {
                  fprintf (stderr, "ERROR: Argument '%s' expects a value. Exiting.\n", arg.c_str());
                  exit(1);
                }
              }
            }
          }
        }
      }
    } else {      // In this case, the argument can either be positional, or misspelled.
      int32_t position = i - argc;

      if (valid_args_positional_.find(position) != valid_args_positional_.end()) {
        arguments_[valid_args_positional_[position]].value = arg;
      } else {
        fprintf (stderr, "ERROR: Unknown parameter value: %s. Exiting.\n", arg.c_str());
        exit(1);
      }
    }
  }
}

std::string ArgumentParser::VerboseArgumentsByGroup() {
  std::map<std::string, std::vector<int32_t>>::iterator it;
  const int32_t value_type_starting_char = 20;
  const int32_t description_starting_char = 25;
  const int32_t wrap_width = 120 - description_starting_char;

  std::stringstream ret_ss;

  // Printout the usage of the program, is program arguments have actually been processed already.
  // If they haven't been, then program_name is empty (program_name == argv[0]).
//  if (program_name_.size() > 0) {
//    ret_ss << "Usage:\n";
//    ret_ss << "  " << program_name_.c_str() << " [options]";
//
//    // Positional arguments need to be sorted in the ascending order of their position.
//    std::vector<Argument> positional_arguments;
//    for (uint32_t i=0; i<arguments_.size(); i++) {
//      if (arguments_.at(i).positional < 0) {
//        positional_arguments.push_back(arguments_.at(i));
//      }
//    }
//    std::sort(positional_arguments.begin(), positional_arguments.end(), positional_less_than_key());
//
//    for (uint32_t i=0; i<positional_arguments.size(); i++) {
//      std::string arg_name = positional_arguments[i].arg_long;
//      if (arg_name == "")
//        arg_name = positional_arguments[i].arg_short;
//      ret_ss << " " << arg_name.c_str();
//    }
//    ret_ss << "\n\n";
//  }

  ret_ss << "Options\n";

  for (it = valid_args_group_.begin(); it != valid_args_group_.end(); it++) {
    ret_ss << "  " << it->first.c_str() << ":\n";
    for (uint32_t i = 0; i < it->second.size(); i++) {
      std::stringstream ss;
      int32_t num_chars = 0;

      ss << "    ";
      num_chars += 4;

      if (arguments_[it->second.at(i)].positional >= 0) {
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

      if (arguments_[it->second.at(i)].value_type != "") {
        int32_t num_empty_chars = (((value_type_starting_char - num_chars) > 0)?(value_type_starting_char - num_chars):1);
        std::string empty_chars(num_empty_chars, ' ');
        ss << empty_chars;

        if (arguments_[it->second.at(i)].value_type.size() >= 3)
          ss << arguments_[it->second.at(i)].value_type.substr(0, 3);
        else {
          ss << arguments_[it->second.at(i)].value_type << std::string((3 - arguments_[it->second.at(i)].value_type.size()), ' ');
        }

        num_chars += num_empty_chars + 3;
      }

      if (arguments_[it->second.at(i)].description_short != "") {
        int32_t num_empty_chars = (((description_starting_char - num_chars) > 0)?(description_starting_char - num_chars):1);
        std::string empty_chars(num_empty_chars, ' ');
        ss << empty_chars;
        ss << WrapString_(description_starting_char, wrap_width, arguments_[it->second.at(i)].description_short);
        num_chars += num_empty_chars + arguments_[it->second.at(i)].description_short.size();
      }
      if (arguments_[it->second.at(i)].value_type != VALUE_TYPE_NONE && arguments_[it->second.at(i)].default_value != "") {
        ss << " [" << arguments_[it->second.at(i)].default_value << "]";
      }

      if (arguments_[it->second.at(i)].arg_short != "" || arguments_[it->second.at(i)].arg_long != "") {
        ret_ss << ss.str() << "\n";
      }
    }

    ret_ss << "\n";
  }

  return ret_ss.str();
}

void ArgumentParser::VerboseArguments(FILE *fp) {
  for (uint32_t i=0; i<arguments_.size(); i++) {
    fprintf (fp, "'-%s'\t'--%s'\t'%s'\tvalue = '%s'\tdefault = '%s'\tpositional = %d\tcount = %d\n", arguments_[i].arg_short.c_str(), arguments_[i].arg_long.c_str(), arguments_[i].value_type.c_str(), arguments_[i].value.c_str(), arguments_[i].default_value.c_str(), arguments_[i].positional, arguments_[i].count);
  }
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

  spaces.push_back(-1);
  for (uint32_t i=0; i<text.size(); i++) {
    if (text[i] == ' ' || text[i] == '\t') {
      spaces.push_back(i);
    }
  }
  spaces.push_back(text.size());

  int32_t line_length = 0;
  for (uint32_t i=1; i<spaces.size(); i++) {
    int32_t word_length = spaces[i] - (spaces[i - 1] + 1);
    if ((i + 1) < spaces.size())
      word_length += 1;

    if ((line_length + word_length + 1) < wrap_width) {

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
  if (arguments_.at(it->second).value_type != VALUE_TYPE_INT)
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
