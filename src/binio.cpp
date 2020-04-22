#include <stdio.h>
#include "binio.h"
#include "Param.h"
/*#include "..\miniz_v111b\miniz.c"
*/
#define COMP_LEVEL 1

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

size_t zfwrite_big(void *buffer,size_t size,size_t count, FILE *stream)
{
	return fwrite_big(buffer,size,count,stream);
	/*
	const unsigned long long mx=0x40000000;
	unsigned long long j,n,st;
	unsigned long cbs,cs,ucs,scs;
	size_t ret=0;
	unsigned char *buf2,*buf3,*tb,*ccbuf,*scbuf,*cbuf1,*cbuf2;

	cbs=mx+0x10000000;
	if(!(cbuf1=(unsigned char *) malloc(((size_t) cbs)*sizeof(unsigned char)))) ERR_CRITICAL("Unable to allocate compression buffer\n");
	if(!(cbuf2=(unsigned char *) malloc(((size_t) cbs)*sizeof(unsigned char)))) ERR_CRITICAL("Unable to allocate compression buffer\n");
	st=mx/((unsigned long long) size);
	n=count/st;
	ucs=(n>0)?((unsigned long) (st*size)):((unsigned long) (count*size));
	ccbuf=cbuf1;
	scbuf=cbuf2;
	buf2=buf3=((unsigned char *)buffer);
	cs=cbs;
	fprintf(stderr,"Compress2 = %i\n",compress2(ccbuf, &cs, buf3, ucs, COMP_LEVEL));
	scs=cs;
	for(j=1;j<=n;j++) 
		{
		tb=scbuf;
		scbuf=ccbuf;
		ccbuf=tb;
		if(j==n) ucs=(unsigned long) ((count-n*st)*size);
		buf3=buf2+ucs;
		cs=cbs;
#pragma omp parallel shared(j,cs,ucs,scs,ret,ccbuf,scbuf,buf3,stream)
			{
#pragma omp sections
				{
#pragma omp section
					{
					fprintf(stderr,"Compressing %i\n",j);
					fprintf(stderr,"Compress2 = %i\n",compress2(ccbuf, &cs, buf3, ucs, COMP_LEVEL));
					}
#pragma omp section
					{
					fprintf(stderr,"Saving %i\n",j);
					fwrite(&scs,sizeof(unsigned long),1,stream);
					ret+=(fwrite(scbuf,sizeof(unsigned char),(size_t) scs,stream));
					}
				}
			}
		scs=cs;
		buf2=buf3;
		}
	scbuf=ccbuf;
	fwrite(&scs,sizeof(unsigned long),1,stream);
	ret+=(fwrite(scbuf,sizeof(unsigned char),(size_t) scs,stream));
	free(cbuf2);
	free(cbuf1);
	return ret;
	*/
}

size_t zfread_big(void *buffer,size_t size,size_t count, FILE *stream)
{
	return fread_big(buffer,size,count,stream);
	/*
	const unsigned long long mx=0x40000000;
	unsigned long long j,n,st;
	unsigned long cbs,cs,ucs,scs,cs2;
	size_t ret=0;
	unsigned char *buf3,*tb,*ccbuf,*scbuf,*cbuf1,*cbuf2;

	cbs=mx+0x10000000;
	if(!(cbuf1=(unsigned char *) malloc(((size_t) cbs)*sizeof(unsigned char)))) ERR_CRITICAL("Unable to allocate compression buffer\n");
	if(!(cbuf2=(unsigned char *) malloc(((size_t) cbs)*sizeof(unsigned char)))) ERR_CRITICAL("Unable to allocate compression buffer\n");
	st=mx/((unsigned long long) size);
	n=count/st;
	ucs=(n>0)?((unsigned long) (st*size)):((unsigned long) (count*size));
	ccbuf=cbuf1;
	scbuf=cbuf2;
	buf3=((unsigned char *)buffer);
	fread(&cs2,sizeof(unsigned long),1,stream);
	fread(scbuf,sizeof(unsigned char),(size_t) cs2,stream);
	scs=cs2;
	for(j=1;j<=n;j++) 
		{
		tb=scbuf;
		scbuf=ccbuf;
		ccbuf=tb;
#pragma omp parallel shared(j,cs,cs2,ucs,scs,ret,ccbuf,scbuf,stream)
			{
#pragma omp sections
				{
#pragma omp section
					{		
					cs=ucs;
					fprintf(stderr,"Uncompress = %i\n",uncompress(buf3,&cs, ccbuf, scs));
					if(cs!=ucs) ERR_CRITICAL("Compressed block size mismatch.\n");
					ret+=ucs;
					}
#pragma omp section
					{
					fread(&cs2,sizeof(unsigned long),1,stream);
					fread(scbuf,sizeof(unsigned char),(size_t) cs2,stream);
					}
				}
			}
		scs=cs2;
		buf3+=ucs;
		}
	ccbuf=scbuf;
	ucs=(unsigned long) ((count-n*st)*size);
	cs=ucs;
	fprintf(stderr,"Uncompress = %i\n",uncompress(buf3,&cs, ccbuf, scs));
	if(cs!=ucs) ERR_CRITICAL("Compressed block size mismatch.\n");
	ret+=ucs;

	free(cbuf2);
	free(cbuf1);
	return ret;
	*/
}
