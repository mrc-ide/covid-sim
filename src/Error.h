#ifndef COVIDSIM_ERROR_H_INCLUDED_
#define COVIDSIM_ERROR_H_INCLUDED_

#ifdef _MSC_VER
__declspec(noreturn)
#endif
void ErrorCritical(const char* msg, ...)
#ifdef __GNUC__
                                            __attribute__ ((noreturn)) __attribute__ ((format (printf, 1, 2)))
#endif
;

#endif // COVIDSIM_ERROR_H_INCLUDED_
