/** \file  Files.cpp
 *  \brief Provide file routines that terminate on failure
 */

#include "Files.h"
#include "Error.h"

FILE* Files::xfopen(const char *filename, const char *mode) noexcept
{
  FILE* fp =  fopen(filename, mode);
  if (fp == NULL) 
  {
    ERR_CRITICAL_FMT("Error %d opening file %s - %s\n", errno, filename, strerror(errno));
  }
  return fp;
}

FILE* Files::xfopen_fallback(const char* filename, const char* mode, FILE* fallback) noexcept
{
  FILE* fp = fopen(filename, mode);
  if (fp == NULL)
  {
    fp = fallback;
  }
  return fp;
}

FILE* Files::xfopen_if_exists(const char* filename, const char* mode) noexcept
{
  return fopen(filename, mode);
}

void Files::xsscanf(const char* s, int n_expected, const char* format, ...) noexcept
{
  va_list args;
  va_start(args, format);
  int rc = vsscanf(s, format, args);
  va_end(args);
  if (rc != n_expected) {
    ERR_CRITICAL_FMT("Error, xsscanf looking for %s in %s, expected %d matches, got %d\n", format, s, n_expected, rc);
  }
}

void Files::xfscanf(FILE* stream, int n_expected, const char* format, ...) noexcept
{
  va_list args;
  va_start(args, format);
  int rc = vfscanf(stream, format, args);
  va_end(args);
  if (rc != n_expected) {
    ERR_CRITICAL_FMT("Error, fsscanf looking for %s, expected %d matches, got %d\n", format, n_expected, rc);
  }
}
