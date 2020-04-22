#pragma once

#ifndef COVIDSIM_BITMAP_H_INCLUDED_
#define COVIDSIM_BITMAP_H_INCLUDED_

#include <stdint.h>

#ifdef UNIX
#define DIRECTORY_SEPARATOR "/"
#else
#define DIRECTORY_SEPARATOR "\\"
#ifndef NO_WIN32_BM
#define WIN32_BM
#endif
#endif
#define STRICT
#ifdef _WIN32
#define _WIN32_WINNT 0x0400
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#ifdef WIN32_BM
#include <vfw.h>
#include <gdiplus.h>
#endif
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

#ifdef WIN32_BM
//DECLARE_HANDLE(HAVI);
//extern HAVI avi;
extern ULONG_PTR m_gdiplusToken;
#endif

void CaptureBitmap();
void OutputBitmap(int);
void InitBMHead();

#endif // COVIDSIM_BITMAP_H_INCLUDED_
