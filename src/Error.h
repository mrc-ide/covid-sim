#pragma once

#ifndef SPATIALSIM_ERROR_H_INCLUDED_
#define SPATIALSIM_ERROR_H_INCLUDED_

#define ERR_CRITICAL(msg) ErrorCritical(msg, __FILE__, __LINE__)
#define ERR_CRITICAL_FMT(fmt, ...) ErrorCritical(fmt, __FILE__, __LINE__, __VA_ARGS__)

#ifdef _MSC_VER
__declspec(noreturn)
#endif
void ErrorCritical(const char* msg, const char* file, int line, ...)
#ifdef __GNUC__
                                            __attribute__ ((noreturn)) __attribute__ ((format (printf, 1, 4)))
#endif
;

#endif // SPATIALSIM_ERROR_H_INCLUDED_
