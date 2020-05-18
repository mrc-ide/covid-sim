#include "Error.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

void ErrorCritical(const char* msg, ...)
{
    const char* file = __FILE__;
    int line = __LINE__;
    va_list args;
    va_start(args, line);
    fprintf(stderr, "[%s line %i] ", file, line);
    vfprintf(stderr, msg, args);
    va_end(args);
    exit(1);
}
