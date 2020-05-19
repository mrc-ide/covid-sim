#include <stdio.h>
#include "BinIO.h"

size_t fwrite_big(void *buffer,size_t size,size_t count, FILE *stream)
{
	const unsigned long long mx=0x80000000;
	unsigned long long j,n,st;
	size_t ret=0;
	char *buf2;

	st=mx/((unsigned long long) size);
	n=count/st;
	for(j=0;j<n;j++)
		{
		buf2=((char *)buffer)+j*st*size;
		ret+=(fwrite(buf2,size,(size_t) st,stream));
		}
	buf2=((char *)buffer)+n*st*size;
	ret+=(fwrite(buf2,size,(size_t) (count-n*st),stream));
	return ret;
}

size_t fread_big(void *buffer,size_t size,size_t count, FILE *stream)
{
	const unsigned long long mx=0x80000000;
	unsigned long long j,n,st;
	size_t ret=0;
	char *buf2;

	st=mx/((unsigned long long) size);
	n=count/st;
	for(j=0;j<n;j++)
		{
		buf2=((char *)buffer)+j*st*size;
		ret+=(fread(buf2,size,(size_t) st,stream));
		}
	buf2=((char *)buffer)+n*st*size;
	ret+=(fread(buf2,size,(size_t) (count-n*st),stream));
	return ret;
}
