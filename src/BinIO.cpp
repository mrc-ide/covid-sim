#include <cstdio>
#include "BinIO.h"

size_t f_big(void *buffer, size_t size, size_t count, FILE *stream, int ftype)
{
	const unsigned long long mx = 0x80000000;
	unsigned long long j, n, st;
	size_t ret = 0;
	char *buf2;

	st = mx / ((unsigned long long) size);
	n = count / st;
	for (j = 0; j < n; j++)
	{
		buf2 = ((char *) buffer) + j*st*size;
		if (ftype == 0)
			ret += (fwrite(buf2, size, (size_t) st, stream));
		else if (ftype == 1)
			ret += (fread(buf2, size, (size_t) st, stream));
	}
	buf2 = ((char *) buffer) + n*st*size;
	if (ftype == 0)
		ret += (fwrite(buf2, size, (size_t) (count - n*st), stream));
	else if (ftype == 1)
		ret += (fread(buf2, size, (size_t) (count - n*st), stream));
	return ret;
}

size_t fwrite_big(void *buffer, size_t size, size_t count, FILE *stream)
{
	return f_big(buffer, size, count, stream, 0);
}

size_t fread_big(void *buffer, size_t size, size_t count, FILE *stream)
{
	return f_big(buffer, size, count, stream, 1);
}
