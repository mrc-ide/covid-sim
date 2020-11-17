#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "BinIO.h"
#include "Bitmap.h"
#include "Error.h"
#include "Param.h"
#include "Model.h"
#include "Memory.h"

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

//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
//// **** BITMAP stuff.

#ifdef _WIN32
//HAVI avi;
static ULONG_PTR m_gdiplusToken;
static HBITMAP bmpdib;
static CLSID  encoderClsid;
#endif

#ifndef _WIN32
#include <sys/stat.h> // for mkdir
#endif

namespace CovidSim
{
	namespace BitMap
	{
		int32_t *Builder::population_ = nullptr;
		int32_t *Builder::infected_   = nullptr;
		int32_t *Builder::recovered_  = nullptr;
		int32_t *Builder::treated_    = nullptr;

		unsigned char *Builder::bitmap_ = nullptr;
		unsigned char *Builder::pixels_ = nullptr;
		unsigned char *Builder::info_   = nullptr;
	}
}

// externs from CovidSim.cpp
// TODO: move these to a header files
extern char OutFile[1024], OutFileBase[1024];

void CovidSim::BitMap::Builder::capture()
{
	int x, y, f, mi;
	unsigned j;
	static double logMaxPop;
	static int fst = 1;
	double prev;

	mi = (int)(bounds_.width * bounds_.height);
	if (fst)
	{
		fst = 0;
		int32_t maxPop = 0;
		for (int i = 0; i < mi; i++) population_[i] = 0;
		for (int i = 0; i < P.PopSize; i++)
		{
			x = ((int)(Households[Hosts[i].hh].loc.x * scale_.x)) - min_.x;
			y = ((int)(Households[Hosts[i].hh].loc.y * scale_.y)) - min_.y;
			if ((x >= 0) && (x < bounds_.width) && (y >= 0) && (y < bounds_.height))
			{
				j = y * header_->width + x;
				if ((j < header_->imagesize) && (j >= 0))
				{
					population_[j]++;
					if (population_[j] > maxPop) maxPop = population_[j];
				}
			}
		}
		logMaxPop = log(1.001 * (double)maxPop);
		for (int i = 0; i < P.NMC; i++)
			if (Mcells[i].n > 0)
			{
				f = 0;
				if ((i < P.NMC - 1) && (i / P.total_microcells_high_ == (i + 1) / P.total_microcells_high_) && (Mcells[i + 1].n > 0) && ((mcell_country[i] != mcell_country[i + 1])
					|| ((P.DoAdunitBoundaryOutput) && ((AdUnits[Mcells[i].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor != (AdUnits[Mcells[i + 1].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor)))) f = 1;
				if ((i > 0) && (i / P.total_microcells_high_ == (i - 1) / P.total_microcells_high_) && (Mcells[i - 1].n > 0) && (mcell_country[i] != mcell_country[i - 1])) f = 1;
				if ((i < P.NMC - P.total_microcells_high_) && (Mcells[i + P.total_microcells_high_].n > 0) && ((mcell_country[i] != mcell_country[i + P.total_microcells_high_])
					|| ((P.DoAdunitBoundaryOutput) && ((AdUnits[Mcells[i].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor != (AdUnits[Mcells[i + P.total_microcells_high_].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor)))) f = 1;
				if ((i >= P.total_microcells_high_) && (Mcells[i - P.total_microcells_high_].n > 0) && (mcell_country[i] != mcell_country[i - P.total_microcells_high_])) f = 1;
				if (f)
				{
					x = (int)(P.in_microcells_.width * (((double)(i / P.total_microcells_high_)) + 0.5) * scale_.x) - min_.x;
					y = (int)(P.in_microcells_.height * (((double)(i % P.total_microcells_high_)) + 0.5) * scale_.y) - min_.y;
					if ((x >= 0) && (x < bounds_.width) && (y >= 0) && (y < bounds_.height))
					{
						j = y * header_->width + x;
						if ((j < header_->imagesize) && (j >= 0)) population_[j] = -1;
					}
				}
			}
		for (int i = 0; i < bounds_.width / 2; i++)
		{
			prev = floor(3.99999 * ((double)i) * BWCOLS / ((double)bounds_.width) * 2);
			f = ((int)prev);
			for (j = 0; j < 10; j++)
			{
				pixels_[(j + bounds_.height + 5) * header_->width + bounds_.width / 4 + i] = f;
			}
		}
	}
#pragma omp parallel for schedule(static,5000) default(none) \
		shared(mi, logMaxPop)
	for (int i = 0; i < mi; i++)
	{
		if (population_[i] == -1)
			pixels_[i] = BWCOLS - 1; /* black for country boundary */
		else if (infected_[i] > 0)
			pixels_[i] = (unsigned char)(BWCOLS + BWCOLS * log((double)infected_[i]) / logMaxPop); /* red for infected */
		else if (treated_[i] > 0)
			pixels_[i] = (unsigned char)(2 * BWCOLS + BWCOLS * log((double)treated_[i]) / logMaxPop); /* blue for treated */
		else if (recovered_[i] > 0)
			pixels_[i] = (unsigned char)(3 * BWCOLS + BWCOLS * log((double)recovered_[i]) / logMaxPop);  /* green for recovered */
		else if (population_[i] > 0)
			pixels_[i] = (unsigned char)(BWCOLS * log((double)population_[i]) / logMaxPop); /* grey for just people */
		else
			pixels_[i] = 0;
	}
}

void CovidSim::BitMap::Builder::output(int tp, std::string const& output_file_base)
{
	char buf[3000];
	int j = 0;
	static int cn[4] = {0, 0, 0, 0};
	const char *OutPrefix[4] = {"", "Mean.", "Min.", "Max."};

	// if the output_file_base has a forward or backwards slash, use that as the cutoff
	// point for the leaf_name, otherwise keep the full output_file_base
	auto slash_loc = output_file_base.find_last_of('/');
	auto bslash_loc = output_file_base.find_last_of('\\');
	if (bslash_loc != std::string::npos && (slash_loc == std::string::npos || bslash_loc > slash_loc)) {
		slash_loc = bslash_loc;
	}
	auto leaf_name = output_file_base;
	if (slash_loc != std::string::npos) {
		leaf_name = leaf_name.substr(slash_loc);
	}

	j = cn[tp];
	cn[tp]++;
	auto OutF = output_file_base + ".ge" DIRECTORY_SEPARATOR + OutPrefix[tp] + leaf_name;

	if (format_ == CovidSim::BitMap::Formats::PNG)
	{
#ifdef IMAGE_MAGICK
	  FILE* dat;
	  using namespace Magick;
	  fprintf(stderr, "\noutputing ImageMagick stuff");
	  sprintf(buf, "%s.bmp", OutF);
	  if (!(dat = fopen(buf, "wb"))) ERR_CRITICAL("Unable to open bitmap file\n");
	  fprintf(dat, "BM");
	  //fwrite_big((void *) &bitmap_,sizeof(unsigned char),(sizeof(bitmap_header)/sizeof(unsigned char))+header_->imagesize,dat);
	  fwrite_big((void*)bitmap_, sizeof(bitmap_header), 1, dat);
	  for (int i = 0; i < header_->imagesize; i++) fputc(pixels_[i], dat);
	  fclose(dat);
	  Image bmap(buf);
	  sprintf(buf, "%s.%d.png", OutF, j);
	  ColorRGB white(1.0, 1.0, 1.0);
	  bmap.transparent(white);
	  bmap.write(buf);
#elif defined(_WIN32)
	  //Windows specific bitmap manipulation code - could be recoded using LIBGD or another unix graphics library
	  using namespace Gdiplus;

	  wchar_t wbuf[1024];
	  size_t a;

	  //Add new bitmap to AVI
	  //if ((P.OutputBitmap == 1) && (tp == 0)) AddAviFrame(avi, bmpdib, (unsigned char*)(&bmh->palette[0][0]));

	  //This transfers HBITMAP to GDI+ Bitmap object
	  Bitmap* gdip_bmp = Bitmap::FromHBITMAP(bmpdib, NULL);
	  //Now change White in palette (first entry) to be transparent
	  if ((cn[0] == 1) && (tp == 0))
	  {
		  static UINT palsize;
		  static ColorPalette* palette;
		  palsize = gdip_bmp->GetPaletteSize();
		  palette = (ColorPalette*)Memory::xcalloc(1, palsize);
		  (void)gdip_bmp->GetPalette(palette, palsize);
		  palette->Flags = PaletteFlagsHasAlpha;
		  palette->Entries[0] = 0x00ffffff; // Transparent white
		  gdip_bmp->SetPalette(palette);
	  }
	  //Now save as png
	  sprintf(buf, "%s.%05i.png", OutF.c_str(), j + 1); //sprintf(buf,"%s.ge" DIRECTORY_SEPARATOR "%s.%05i.png",OutFileBase.c_str(),OutF.c_str(),j+1);
	  mbstowcs_s(&a, wbuf, strlen(buf) + 1, buf, _TRUNCATE);
	  gdip_bmp->Save(wbuf, &encoderClsid, NULL);
	  delete gdip_bmp;
#else
	  fprintf(stderr, "Do not know how to output PNG\n");
#endif
	}
	else if (format_ == CovidSim::BitMap::Formats::BMP) {
	  sprintf(buf, "%s.%05i.bmp", OutF.c_str(), j);
	  FILE* dat;
	  if (!(dat = fopen(buf, "wb"))) {
	    char* errMsg = strerror(errno);
	    if (errMsg == nullptr) {
	      ERR_CRITICAL("strerror failed.\n");
	    }
	    ERR_CRITICAL_FMT("Unable to open bitmap file %s (%d): %s\n", buf, errno, errMsg);
	  }
	  fprintf(dat, "BM");
	  fwrite_big((void*)bitmap_, sizeof(unsigned char), sizeof(CovidSim::BitMap::Header) / sizeof(unsigned char) + header_->imagesize, dat);
	  fclose(dat);
	}
	else
	{
	  fprintf(stderr, "Unknown Bitmap format: %d\n", (int)format_);
	}
}
void CovidSim::BitMap::Builder::initialise_header(std::string const& out_file_base)
{
	int i, j, k, k2, value;

	fprintf(stderr, "Initialising bitmap\n");
	k = bounds_.width * height2_;
	k2 = sizeof(CovidSim::BitMap::Header) / sizeof(unsigned char);

	bitmap_ = (unsigned char*)Memory::xcalloc((size_t)k + k2, sizeof(unsigned char));
	pixels_ = &(bitmap_[k2]);
	info_ = &(bitmap_[12]);
	header_ = (CovidSim::BitMap::Header*)bitmap_;
	header_->spare = 0;
	header_->boffset = 2 + sizeof(CovidSim::BitMap::Header);
	header_->headersize = 40; // BITMAPINFOHEADER
	header_->width = bounds_.width;
	header_->height = height2_;
	header_->PlanesAndBitspp = 1 // Number of colour planes; must be 1
	                     + (8 << 16); // Colour depth: 8 bits per pixel
	header_->compr = 0; // No compression (BI_RGB)
	header_->imagesize = header_->width * header_->height;
	header_->filesize = 2 // "BM"
	              + ((unsigned int) sizeof(CovidSim::BitMap::Header)) // BITMAP_HEADER
	              + header_->imagesize; // Image data
	header_->hres = header_->vres = (int)(header_->width * 10); // Resolution, in pixels per metre
	header_->colours = BWCOLS * 4; // Number of colours in the palette
	header_->impcol = 0; // Every colour is important
	for (i = 0; i < BWCOLS * 4; i++)
		header_->palette[i][3] = 0;
	for (j = 0; j < BWCOLS; j++)
	{
		value = 255 - 255 * j / (BWCOLS - 1);
		// Shades of gray:
		header_->palette[j][0] = header_->palette[j][1] = header_->palette[j][2] = (unsigned char)value;
		// Shades of red:
		header_->palette[BWCOLS + j][0] = 0;
		header_->palette[BWCOLS + j][1] = 0;
		header_->palette[BWCOLS + j][2] = (unsigned char)value;
		// Shades of blue:
		header_->palette[2 * BWCOLS + j][0] = (unsigned char)value;
		header_->palette[2 * BWCOLS + j][1] = 0;
		header_->palette[2 * BWCOLS + j][2] = 0;
		// Shades of green:
		header_->palette[3 * BWCOLS + j][0] = 0;
		header_->palette[3 * BWCOLS + j][1] = (unsigned char)value;
		header_->palette[3 * BWCOLS + j][2] = 0;
	}
	population_ = (int32_t*)Memory::xcalloc(header_->imagesize, sizeof(int32_t));
	infected_ = (int32_t*)Memory::xcalloc(header_->imagesize, sizeof(int32_t));
	recovered_ = (int32_t*)Memory::xcalloc(header_->imagesize, sizeof(int32_t));
	treated_ = (int32_t*)Memory::xcalloc(header_->imagesize, sizeof(int32_t));

#ifdef _WIN32
	if (format_ == CovidSim::BitMap::Formats::PNG)
	{
	  bmpdib = CreateDIBSection(GetDC(NULL), (BITMAPINFO*)info_, DIB_RGB_COLORS, (void**)&pixels_, NULL, NULL);
	  Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	  Gdiplus::GdiplusStartup(&m_gdiplusToken, &gdiplusStartupInput, NULL);

	  UINT  num = 0;          // number of image encoders
	  UINT  size = 0;         // size of the image encoder array in bytes

	  Gdiplus::ImageCodecInfo* pImageCodecInfo = NULL;
	  Gdiplus::GetImageEncodersSize(&num, &size);
	  pImageCodecInfo = (Gdiplus::ImageCodecInfo*)Memory::xcalloc(1, size);
	  Gdiplus::GetImageEncoders(num, size, pImageCodecInfo);
	  for (UINT j = 0; j < num; ++j) {
	    // Visual Studio Analyze incorrectly reports this because it doesn't understand Gdiplus::GetImageEncodersSize()
	    // warning C6385: Reading invalid data from 'pImageCodecInfo':  the readable size is 'size' bytes, but '208' bytes may be read.
#pragma warning( suppress: 6385 )
	    const WCHAR* type = pImageCodecInfo[j].MimeType;
	    if (wcscmp(type, L"image/png") == 0) {
	      encoderClsid = pImageCodecInfo[j].Clsid;
	      j = num;
	    }
	  }
	  free(pImageCodecInfo);
	}
#endif

	char buf[1024+3];
	sprintf(buf, "%s.ge", out_file_base.c_str());
#ifdef _WIN32
	if (!(CreateDirectory(buf, NULL))) fprintf(stderr, "Unable to create directory %s\n", buf);
#else
	if (!(mkdir(buf, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))) fprintf(stderr, "Unable to create directory %s\n", buf);
#endif
}

void CovidSim::BitMap::Builder::finalise()
{
#ifdef _WIN32
  if (format_ == CovidSim::BitMap::Formats::PNG)
  {
    Gdiplus::GdiplusShutdown(m_gdiplusToken);
  }
#endif
}
