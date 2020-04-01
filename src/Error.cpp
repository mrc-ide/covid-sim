#include "Error.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

void ErrorCritical(const char* fmt, const char* file, int line, ...)
{
    va_list args;
    va_start(args, line);
    fprintf(stderr, "[%s line %i] ", file, line);
    vfprintf(stderr, fmt, args);
    va_end(args);
    exit(1);
}
