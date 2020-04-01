#pragma once

#ifndef SPATIALSIM_BITMAP_H_INCLUDED_
#define SPATIALSIM_BITMAP_H_INCLUDED_

#ifdef UNIX
#define DIRECTORY_SEPARATOR "/"
#else
#define DIRECTORY_SEPARATOR "\\"
#define WIN32_BM
#endif
#define STRICT
#ifdef WIN32_BM
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

extern float* bmi, *bmi2, *bmi3, *bmi4;
extern bitmap_header* bmh;

#ifdef WIN32_BM
//DECLARE_HANDLE(HAVI);
//extern HAVI avi;
extern ULONG_PTR m_gdiplusToken;
#endif

void CaptureBitmap(int, int);
void OutputBitmap(double, int);
void InitBMHead();

#endif // SPATIALSIM_BITMAP_H_INCLUDED_
