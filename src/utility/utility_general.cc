/*
 * utility_general.cc
 *
 *  Created on: Oct 11, 2014
 *      Author: ivan
 */

#include <unistd.h>
#include <ios>
#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include "utility/utility_general.h"
//#include "utility/utility_basic_defines.h"
#include "log_system/log_system.h"


int fileExists(const char *fname) {
     /* from
      * http://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c-cross-platform
      */
     if (access(fname, F_OK) != -1) {
          return 1;
     } else {
          return 0;
     }
}

std::string GetUTCTime(std::string fmt) {
  //  const char* fmt = "%a, %d %b %y %T %z";
  char outstr[200];
  time_t t;
  struct tm *tmp;

  t = time(NULL);
  tmp = gmtime(&t);

  if (tmp == NULL) {
    fprintf (stderr, "ERROR: gmtime returned with error!\n");
    fflush (stderr);
    return std::string("[no_time]");
  }

  if (strftime(outstr, sizeof(outstr), fmt.c_str(), tmp) == 0) {
    fprintf (stderr, "ERROR: Problem formatting time into string!\n");
    fflush (stderr);
    return std::string("[no_time]");
  }

  return std::string(outstr);
}

std::string GetLocalTime() {
  char outstr[200];

  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  sprintf(outstr, "%0.2d:%0.2d:%0.2d", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

  return std::string(outstr);
}



//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
//
// Source: http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c

void ProcessMemUsage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

/*
void PrintRSSUsage() {
//  double vm=0.0f, rss=0.0f;
//  ProcessMemUsage(vm, rss);
//  printf ("[VMsize = %.2f, RSS = %.2f]\n");
  #ifndef RELEASE_VERSION
    printf ("[currentRSS = %ld MB, peakRSS = %ld MB]\n", getCurrentRSS()/(1024*1024), getPeakRSS()/(1024*1024));
    fflush(stdout);
  #endif
}
*/

std::string FormatMemoryConsumptionAsString() {
  std::stringstream ss;
  ss << "[currentRSS = " << getCurrentRSS()/(1024*1024) << " MB, peakRSS = " << getPeakRSS()/(1024*1024) << " MB]";
  return ss.str();
}



/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif





/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}

// The function uses vasprintf for formatting the function arguments.
// This option was chosen because it preallocates the memory needed to store
// the output formatted string, unlike most other *printf alternatives
std::string FormatString(const char* additional_message, ...) {
  char *formatted_string = NULL;
  std::string return_message("");

  va_list args;
  va_start(args, additional_message);
  int ret_va = vasprintf(&formatted_string, additional_message, args);
  if (ret_va == -1) {
    fprintf (stderr, "ERROR: Could not format the string!\n");
    va_end(args);
    return std::string("");
  }
  va_end(args);

  if (formatted_string == NULL) {
//    LOG(FATAL) << ErrorReporting::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Allocation of memory for variable 'formatted_string'.");
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Allocation of memory for variable 'formatted_string'."));
    return ((std::string) "");
  }

  return_message = std::string(formatted_string);
  if (formatted_string)
    free(formatted_string);
  formatted_string = NULL;

  return return_message;
}

std::string FormatStringToLength(std::string original_string, uint32_t length) {
  std::string ret;
  if (original_string.size() > length) {
    ret = original_string.substr(0, length - 3) + "...";
  } else {
    ret = original_string;
  }

  if (ret.size() < length)
    ret += std::string((length - ret.size()), ' ');

  return ret;
}

void PrintSubstring(char *text, int64_t length, FILE *fp) {
  for (int64_t i=0; i<length; i++)
    fprintf (fp, "%c", text[i]);
}

std::string GetSubstring(char *text, int64_t length) {
  std::stringstream ss;
  for (int64_t i=0; i<length; i++)
    ss << text[i];
  return ss.str();
}

/**
 * Returns new sequence that is reverse of given sequence.
 */
unsigned char* CreateReverseCopy(const unsigned char* seq, uint64_t length) {
    unsigned char* reverse_seq = new unsigned char[length];
    for (uint64_t i = 0; i < length; i++) {
      reverse_seq[i] = seq[length - i - 1];
    }
    return reverse_seq;
}

std::string TrimToFirstSpace(std::string original_string) {
  std::string::size_type loc = original_string.find(" ", 0);
  if (loc != std::string::npos) {
    return original_string.substr(0, loc);

  } else {
    // There is no spaces in the string, do nothing and just report it as is.
    return original_string;
  }
  return std::string("");
}

float sigmoid(float x, float mean, float width) {
  float t = ((x - mean) / width) * 6.0f;
  float ret = 1.0f / (1.0f + exp(-t));

  return ret;
}

int GetClippingOpsFromCigar(const std::string &cigar, char *clip_op_front, int64_t *clip_count_front, char *clip_op_back, int64_t *clip_count_back) {
  if (cigar.size() < 2)
    return 1;

  int64_t pos_clip_op_front = -1;
  for (int64_t i=0; i<cigar.size(); i++) {
    if (!isdigit(cigar[i])) {
      if (cigar[i] == 'S' || cigar[i] == 'H') {
        *clip_op_front = cigar[i];
        *clip_count_front = atoi(cigar.substr(0, i).c_str());
        pos_clip_op_front = i;
      }
      break;
    }
  }

  if ((pos_clip_op_front + 1) < cigar.size()) {
    if (cigar.back() == 'S' || cigar.back() == 'H') {
      *clip_op_back = cigar.back();
      for (int64_t i=(cigar.size()-2); i>=0; i--) {
        if (!isdigit(cigar[i])) {
          *clip_count_back = atoi(cigar.substr((i + 1), ((cigar.size()-1) - i - 1)).c_str());
          break;
        }
      }
    }
  }

  return 0;
}
