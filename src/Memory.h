/** \file  Memory.h
 *  \brief Provide memory routines that terminate on failure
 */

#ifndef MEMORY_H_INCLUDED_
#define MEMORY_H_INCLUDED_

#include <cstddef>

namespace Memory
{
/** \brief       Wrapper around malloc that aborts on error.
 *  \param  size Amount of memory to allocate
 *  \return      Pointer to allocated memory
 *
 *  If \a size is 0 we return a pointer to a small amount of allocated memory
 *  which should not be used.
 */
void* xmalloc(std::size_t size) noexcept;

/** \brief         Wrapper around calloc that aborts on error.
 *  \param  nelem  Number of elements to allocate
 *  \param  elsize Size of each element
 *  \return        Pointer to allocated memory
 *
 * If \a nelem or \a elsize are 0 we return a pointer to a small amount of
 * allocated memory - which should not be used.
 */
void* xcalloc(std::size_t nelem, std::size_t elsize) noexcept;

/** \brief       Wrapper around realloc that aborts on error.
 *  \param  size Size of new buffer.
 *  \return      Pointer to allocated memory
 *
 *  If \a size is 0 we return a pointer to a small amount of allocated memory
 *  which should not be used.
 */
void* xrealloc(void* ptr, std::size_t size) noexcept;

/** \brief     Wrapper around free that ignores all errors.
 *  \param ptr Pointer to memory to free.
 */
void xfree(void* ptr) noexcept;
} // namespace Memory

#endif // MEMORY_H_INCLUDED_
