/*
 * ErrorHandling.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: ivan
 */

#include "log_system/log_system.h"
#include <sstream>
//#include <glog/logging.h>
//#include <glog/log_severity.h>

//extern uint32_t LOG_VERBOSE_TYPE;

LogSystem::LogSystem() {
  LOG_FILE = "graphmap.log";
  LOG_VERBOSE_TYPE = 0;
  PROGRAM_VERBOSE_LEVEL = 0;
}

LogSystem::~LogSystem() {

}

LogSystem& LogSystem::GetInstance() {
    static LogSystem    instance;  // Guaranteed to be destroyed.
    return instance;
}

std::string LogSystem::GenerateErrorMessage(uint32_t error_type,
                                                 const char* additional_message,
                                                 ...) {
  char *formatted_c_string = NULL;
  std::string formatted_string(""), return_message("");

  // Process format arguments given to the function, to format the
  // additional message.
  va_list args;
  va_start(args, additional_message);
  int ret_va = vasprintf(&formatted_c_string, additional_message, args);
  if (ret_va == -1) {
    fprintf (stderr, "ERROR: could not format the string for the error message!\n");
    fprintf (stderr, "Exiting.\n");
    exit(1);
  }
  va_end(args);

  // Check if the memory has been allocated properly. Here we can not
  // call an error handling/reporting function since this is an error
  // handling/reporting function. :)
  if (formatted_c_string == NULL) {
    char error_type_as_string[32];
    sprintf(error_type_as_string, "%d", ERR_MEMORY);

//    LOG(FATAL)<< std::string("#") << std::string(error_type_as_string)
//    << std::string(": ") << "Memory assertion failure. Possible cause - not enough memory or memory not allocated.";
    std::stringstream ss;
    ss << std::string("#") << std::string(error_type_as_string)
    << std::string(": ") << "Memory assertion failure. Possible cause - not enough memory or memory not allocated.";
    Log(SEVERITY_INT_FATAL, __FUNCTION__, ss.str());

    return ((std::string) "");
  }

  formatted_string = std::string(formatted_c_string);

  if (formatted_c_string)
    free(formatted_c_string);
  formatted_c_string = NULL;

  // Format the beginning of the return message, which consists of the
  // error number (error_type variable).
  char error_type_as_string[32];
  sprintf(error_type_as_string, "%d", error_type);
  return_message = std::string("#") + std::string(error_type_as_string)
      + std::string(": ");

  // Attach a custom error message for a specific error type.
  if (error_type == ERR_MEMORY)
    return_message +=
        std::string(
            "Memory assertion failure. Possible cause - not enough memory or memory not allocated.");
  else if (error_type == ERR_OPENING_FILE)
    return_message += std::string("Could not open file.");
  else if (error_type == ERR_CLOSING_FILE)
    return_message += std::string("Could not close file.");
  else if (error_type == ERR_FILE_NOT_FOUND)
    return_message += std::string("File not found.");
  else if (error_type == ERR_FOLDER_NOT_FOUND)
    return_message += std::string("Folder not found.");
  else if (error_type == ERR_OTHER) {
  }
  else if (error_type == ERR_UNEXPECTED_VALUE)
    return_message += std::string("Unexpected value found!");
  else if (error_type == ERR_FILE_WRITE_DATA)
    return_message += std::string("Could not write all data to file!");
  else if (error_type == ERR_FILE_READ_DATA)
    return_message += std::string("Could not get expected data from file, not all data loaded!");
  else if (error_type == ERR_WRONG_FILE_TYPE)
    return_message += std::string("File structure is not in correct format!");
  else if (error_type == ERR_NOT_IMPLEMENTED)
    return_message += std::string("Function not yet implemented!");
  else if (error_type == ERR_SEQUENCE_MISMATCH)
    return_message += std::string("Sequences do not match.");
  else if (error_type == ERR_WRONG_PARAMS)
    return_message += std::string(
        "Wrong function parameters, or parameters out of bounds!");
  else if (error_type == ERR_FILE_DEFORMED_FORMAT)
    return_message += std::string(
        "Deformed file format, format not following specification.");

  // Attach the formatted additional message to the end of the return message.
  return_message += (std::string(" ")) + formatted_string;

  return return_message;
}

int LogSystem::WriteLog(std::string log_entry, bool always_output_to_std) {
  if ((LOG_VERBOSE_TYPE & LOG_VERBOSE_STD) != 0 || always_output_to_std == true) {
    fprintf (stderr, "%s\n", log_entry.c_str());
    fflush(stderr);
  }

  if ((LOG_VERBOSE_TYPE & LOG_VERBOSE_FILE) != 0 && LOG_FILE != "") {
    FILE *fp = fopen(LOG_FILE.c_str(), "a");
    if (fp == NULL) {
      fprintf (stderr, "ERROR: Could not open log file '%s' for writing!\n", LOG_FILE.c_str());
      fprintf (stderr, "Original message that was supposed to be output to the log:\n");
      fprintf (stderr, "%s\n", log_entry.c_str());
      fflush(stderr);
      return 1;
    }
    fprintf (fp, "%s\n", log_entry.c_str());
    fflush(stdout);
    fflush(stderr);
    fclose(fp);
  }

  return 0;
}

int LogSystem::Log(int severity, std::string function, std::string message) {
  std::string timestamp = GetUTCTime();
  bool is_critical = (severity == SEVERITY_INT_ERROR || severity == SEVERITY_INT_FATAL);

  const char *severity_lookup[] = {"INFO", "WARNING", "ERROR", "FATAL"};

  std::stringstream ss;
  ss << "[" << timestamp << " ";
  if (severity == SEVERITY_INT_WARNING)
    ss << severity_lookup[1] << "] ";
  else if (severity == SEVERITY_INT_ERROR)
    ss << severity_lookup[2] << "] ";
  else if (severity == SEVERITY_INT_FATAL)
    ss << severity_lookup[3] << "] ";
  else
    ss << severity_lookup[0] << "] ";

  // << message;

  WriteLog(ss.str() + message, is_critical);

  if (is_critical) {
    WriteLog(ss.str() + std::string("Exiting."), is_critical);
    exit(0);
  }

  return 0;
}

//ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL,
//                                         true,
//                                         FormatString(""));

int LogSystem::VerboseLog(uint32_t verbose_level, bool trigger_condition, std::string message, std::string message_header) {

  if (trigger_condition == false)
    return 1;

  // Although we check for the same condition below, this is a security sanity check, so no debug
  // output escapes in the release version.
  if ((PROGRAM_VERBOSE_LEVEL & VERBOSE_LEVEL_DEBUG) == 0 && (verbose_level & VERBOSE_LEVEL_DEBUG) != 0) {
    return 2;
  }

  // Second security check.
  #ifdef RELEASE_VERSION
    if ((verbose_level & VERBOSE_LEVEL_DEBUG) != 0) {
      return 3;
    }
  #endif

  FILE *fp = stderr;

  std::string timestamp = GetLocalTime();

  if (!(PROGRAM_VERBOSE_LEVEL & VERBOSE_LEVEL_ALL)) {
    if ((verbose_level & VERBOSE_LEVEL_FORCE) != 0) {
      if (message_header == "[]") {
        fprintf (fp, "%s", message.c_str());
      }
      else {
        if (message_header.size() > 0)
          message_header += " ";

        if (message[0] != '\r')
          fprintf (fp, "[%s%s] %s", message_header.c_str(), timestamp.c_str(), message.c_str());
        else
          fprintf (fp, "\r[%s%s] %s", message_header.c_str(), timestamp.c_str(), message.substr(1).c_str());

      }
      fflush(fp);
    }

    return 0;
  }

  if (((verbose_level & VERBOSE_LEVEL_LOW) && (PROGRAM_VERBOSE_LEVEL & VERBOSE_LEVEL_LOW)) ||
      ((verbose_level & VERBOSE_LEVEL_MED) && (PROGRAM_VERBOSE_LEVEL & VERBOSE_LEVEL_MED)) ||
      ((verbose_level & VERBOSE_LEVEL_HIGH) && (PROGRAM_VERBOSE_LEVEL & VERBOSE_LEVEL_HIGH))) {
    if ((PROGRAM_VERBOSE_LEVEL & VERBOSE_LEVEL_DEBUG) ||      // In debug mode output everything.
        ((PROGRAM_VERBOSE_LEVEL & VERBOSE_LEVEL_DEBUG) == (verbose_level & VERBOSE_LEVEL_DEBUG))) { //Otherwise, only output those that aren't debug messages.
      if (message_header == "[]") {
        fprintf (fp, "%s", message.c_str());
      }
      else {
        if (message_header.size() > 0)
          message_header += " ";

        if (message[0] != '\r')
          fprintf (fp, "[%s%s] %s", message_header.c_str(), timestamp.c_str(), message.c_str());
        else
          fprintf (fp, "\r[%s%s] %s", message_header.c_str(), timestamp.c_str(), message.substr(1).c_str());
      }
      fflush(fp);
    }
  }

  return 0;
}

void LogSystem::SetProgramVerboseLevelFromInt(int64_t verbose_level) {
  if (verbose_level == 0)
    PROGRAM_VERBOSE_LEVEL = 0;
  else if (verbose_level == 1)
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_LOW | VERBOSE_FREQ_LOW;
  else if (verbose_level == 2)
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_MED | VERBOSE_FREQ_MED;
  else if (verbose_level == 3)
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_HIGH | VERBOSE_FREQ_HIGH;
  else if (verbose_level == 4)
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_ALL | VERBOSE_FREQ_MED;
  else if (verbose_level == 5)
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_ALL | VERBOSE_FREQ_HIGH;
#ifndef RELEASE_VERSION
  else if (verbose_level == 6)
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_LOW | VERBOSE_LEVEL_DEBUG | VERBOSE_FREQ_ALL;
  else if (verbose_level == 7)
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_MED | VERBOSE_LEVEL_DEBUG | VERBOSE_FREQ_ALL;
  else if (verbose_level == 8)
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_DEBUG | VERBOSE_FREQ_ALL;
  else if (verbose_level == 9)
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_ALL | VERBOSE_LEVEL_DEBUG | VERBOSE_FREQ_ALL;
#endif
  else
    PROGRAM_VERBOSE_LEVEL = VERBOSE_LEVEL_ALL | VERBOSE_FREQ_ALL;
}

//std::string ErrorReporting::GenerateMessage(const char *message, ...) {
//  char *formatted_c_string = NULL;
//
//  // Process format arguments given to the function, to format the
//  // additional message.
//  va_list args;
//  va_start(args, message);
//  int ret_va = vasprintf(&formatted_c_string, message, args);
//  if (ret_va == -1) {
//    fprintf (stderr, "ERROR: could not format the string!\n");
//    va_end(args);
//    return std::string("");
//  }
//  va_end(args);
//
//  // Check if the memory has been allocated properly.
//  if (formatted_c_string == NULL) {
//    fprintf (stderr, "ERROR: Could not allocate memory!\n");
//    return std::string("");
//  }
//
//  std::string ret = std::string(formatted_c_string);
//
//  if (formatted_c_string)
//    free(formatted_c_string);
//  formatted_c_string = NULL;
//
//  return ret;
//}
