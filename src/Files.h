/** \file  Files.h
 *  \brief Provide file routines that terminate on failure
 */

#ifndef FILES_H_INCLUDED_
#define FILES_H_INCLUDED_

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <cstdarg>

namespace Files
{
/** \brief           Wrapper around fopen that aborts on error.
 *  \param  filename The filename to open
 *  \param  mode     Mode, as per fopen (eg, "wb")
 *  \return          File handle
 *
 */

  FILE* xfopen(const char *filename, const char *mode) noexcept;

/** \brief            Wrapper around fopen, returning a fallback option if the call fails
 *  \param  filename  The filename to open
 *  \param  mode      Mode, as per fopen (eg, "wb")
 *  \param  fallback  A file stream to return, should the call fail.
 *  \return           File handle
 *
 */

  FILE* xfopen_fallback(const char* filename, const char* mode, FILE* fallback) noexcept;

/** \brief            Wrapper around fopen, safely returning NULL if the file does not exist.
 *  \param  filename  The filename to open
 *  \param  mode      Mode, as per fopen (eg, "wb")
 *  \return           File handle - or NULL if file does not exist
 *
 */

  FILE* xfopen_if_exists(const char* filename, const char* mode) noexcept;

/** \brief             Wrapper around sscanf that aborts on error.
 *  \param  s          The string to be parsed
 *  \param  format     The format string for sscanf
 *  \param  ...        Destinations for each placeholder
 *
 */

  void xsscanf(const char* s, int n_expected, const char* format, ...) noexcept;

/** \brief             Wrapper around fscanf that aborts on error.
 *  \param  stream     The file stream to read from
 *  \param  n_expected Expected number of arguments
 *  \param  format     The format string for fscanf
 *  \param  ...        Destinations for each placeholder
 *
 */

  void xfscanf(FILE* stream, int n_expected, const char* format, ...) noexcept;

/** \brief             Write data with compatibility for large types
 *  \param  buffer     Buffer to write from
 *  \param  size       How many elements
 *  \param  count      Size of each element
 *  \param  stream     File stream to write to
 *  \return            Number of elements written
 *
 */

  size_t fwrite_big(void* buffer, size_t size, size_t count, FILE* stream);
  
  /** \brief             Read data with compatibility for large types
   *  \param  buffer     Buffer to write from
   *  \param  size       How many elements
   *  \param  count      Size of each element
   *  \param  stream     File stream to read from
   *  \return            Number of elements read
   */

  size_t fread_big(void* buffer, size_t size, size_t count, FILE* stream);

  /** \brief             Wrapper around rename that aborts on error.
   *  \param  oldname    Original name of file
   *  \param  newname    How many elements
   *  \param  count      Size of each element
   *  \param  stream     File stream to read from
   *
   */

  void xrename(const char* oldname, const char* newname) noexcept;

/** \brief             Wrapper around fclose that aborts on error.
 *  \param  stream     File stream to close
 */

   void xfclose(FILE* stream) noexcept;

} // namespace Files

#endif // FILES_H_INCLUDED_
