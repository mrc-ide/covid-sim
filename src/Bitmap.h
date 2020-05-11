#ifndef COVIDSIM_BITMAP_H_INCLUDED_
#define COVIDSIM_BITMAP_H_INCLUDED_

#include <stdint.h>

#ifdef UNIX
#define DIRECTORY_SEPARATOR "/"
#else
#define DIRECTORY_SEPARATOR "\\"
#endif

#define STRICT
#ifdef _WIN32
#define _WIN32_WINNT 0x0400
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <vfw.h>
#include <gdiplus.h>
#endif
#ifdef IMAGE_MAGICK
#include "Magick++.h"
#endif

#define BWCOLS 58

typedef struct BITMAP_HEADER {
	unsigned int filesize;
	unsigned int spare;
	unsigned int boffset;
	unsigned int headersize;
	unsigned int width, height;
	unsigned int PlanesAndBitspp;
	unsigned int compr;
	unsigned int imagesize;
	unsigned int hres, vres;
	unsigned int colours, impcol;
	unsigned char palette[BWCOLS * 4][4];
} bitmap_header;

extern int32_t *bmPopulation, *bmInfected, *bmRecovered, *bmTreated;
extern bitmap_header* bmh;

void CaptureBitmap();
void OutputBitmap(int);
void InitBMHead();

void Bitmap_Finalise();

#endif // COVIDSIM_BITMAP_H_INCLUDED_
