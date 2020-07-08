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

static unsigned char* bmf, *bmPixels, *bmp;

void CaptureBitmap()
{
	int x, y, f, mi;
	unsigned j;
	static double logMaxPop;
	static int fst = 1;
	double prev;

	mi = (int)(P.b.width * P.b.height);
	if (fst)
	{
		fst = 0;
		int32_t maxPop = 0;
		for (int i = 0; i < mi; i++) bmPopulation[i] = 0;
		for (int i = 0; i < P.PopSize; i++)
		{
			x = ((int)(Households[Hosts[i].hh].loc.x * P.scale.x)) - P.bmin.x;
			y = ((int)(Households[Hosts[i].hh].loc.y * P.scale.y)) - P.bmin.y;
			if ((x >= 0) && (x < P.b.width) && (y >= 0) && (y < P.b.height))
			{
				j = y * bmh->width + x;
				if ((j < bmh->imagesize) && (j >= 0))
				{
					bmPopulation[j]++;
					if (bmPopulation[j] > maxPop) maxPop = bmPopulation[j];
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
					x = (int)(P.in_microcells_.width * (((double)(i / P.total_microcells_high_)) + 0.5) * P.scale.x) - P.bmin.x;
					y = (int)(P.in_microcells_.height * (((double)(i % P.total_microcells_high_)) + 0.5) * P.scale.y) - P.bmin.y;
					if ((x >= 0) && (x < P.b.width) && (y >= 0) && (y < P.b.height))
					{
						j = y * bmh->width + x;
						if ((j < bmh->imagesize) && (j >= 0)) bmPopulation[j] = -1;
					}
				}
			}
		for (int i = 0; i < P.b.width / 2; i++)
		{
			prev = floor(3.99999 * ((double)i) * BWCOLS / ((double)P.b.width) * 2);
			f = ((int)prev);
			for (j = 0; j < 10; j++)
			{
				bmPixels[(j + P.b.height + 5) * bmh->width + P.b.width / 4 + i] = f;
			}
		}
	}
#pragma omp parallel for schedule(static,5000) default(none) \
		shared(mi, bmPixels, bmPopulation, bmInfected, bmTreated, bmRecovered, logMaxPop)
	for (int i = 0; i < mi; i++)
	{
		if (bmPopulation[i] == -1)
			bmPixels[i] = BWCOLS - 1; /* black for country boundary */
		else if (bmInfected[i] > 0)
			bmPixels[i] = (unsigned char)(BWCOLS + BWCOLS * log((double)bmInfected[i]) / logMaxPop); /* red for infected */
		else if (bmTreated[i] > 0)
			bmPixels[i] = (unsigned char)(2 * BWCOLS + BWCOLS * log((double)bmTreated[i]) / logMaxPop); /* blue for treated */
		else if (bmRecovered[i] > 0)
			bmPixels[i] = (unsigned char)(3 * BWCOLS + BWCOLS * log((double)bmRecovered[i]) / logMaxPop);  /* green for recovered */
		else if (bmPopulation[i] > 0)
			bmPixels[i] = (unsigned char)(BWCOLS * log((double)bmPopulation[i]) / logMaxPop); /* grey for just people */
		else
			bmPixels[i] = 0;
	}
}

void OutputBitmap(int tp, std::string const& output_file_base)
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

	if (P.BitmapFormat == BitmapFormats::PNG)
	{
#ifdef IMAGE_MAGICK
	  FILE* dat;
	  using namespace Magick;
	  fprintf(stderr, "\noutputing ImageMagick stuff");
	  sprintf(buf, "%s.bmp", OutF);
	  if (!(dat = fopen(buf, "wb"))) ERR_CRITICAL("Unable to open bitmap file\n");
	  fprintf(dat, "BM");
	  //fwrite_big((void *) &bmf,sizeof(unsigned char),(sizeof(bitmap_header)/sizeof(unsigned char))+bmh->imagesize,dat);
	  fwrite_big((void*)bmf, sizeof(bitmap_header), 1, dat);
	  for (int i = 0; i < bmh->imagesize; i++) fputc(bmPixels[i], dat);
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
	else if (P.BitmapFormat == BitmapFormats::BMP) {
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
	  fwrite_big((void*)bmf, sizeof(unsigned char), sizeof(BitmapHeader) / sizeof(unsigned char) + bmh->imagesize, dat);
	  fclose(dat);
	}
	else
	{
	  fprintf(stderr, "Unknown Bitmap format: %d\n", (int)P.BitmapFormat);
	}
}
void InitBMHead(std::string const& out_file_base)
{
	int i, j, k, k2, value;

	fprintf(stderr, "Initialising bitmap\n");
	k = P.b.width * P.bheight2;
	k2 = sizeof(BitmapHeader) / sizeof(unsigned char);

	bmf = (unsigned char*)Memory::xcalloc((size_t)k + k2, sizeof(unsigned char));
	bmPixels = &(bmf[k2]);
	bmp = &(bmf[12]);
	bmh = (BitmapHeader*)bmf;
	bmh->spare = 0;
	bmh->boffset = 2 + sizeof(BitmapHeader);
	bmh->headersize = 40; // BITMAPINFOHEADER
	bmh->width = P.b.width;
	bmh->height = P.bheight2;
	bmh->PlanesAndBitspp = 1 // Number of colour planes; must be 1
	                     + (8 << 16); // Colour depth: 8 bits per pixel
	bmh->compr = 0; // No compression (BI_RGB)
	bmh->imagesize = bmh->width * bmh->height;
	bmh->filesize = 2 // "BM"
	              + ((unsigned int) sizeof(BitmapHeader)) // BITMAP_HEADER
	              + bmh->imagesize; // Image data
	bmh->hres = bmh->vres = (int)(bmh->width * 10); // Resolution, in pixels per metre
	bmh->colours = BWCOLS * 4; // Number of colours in the palette
	bmh->impcol = 0; // Every colour is important
	for (i = 0; i < BWCOLS * 4; i++)
		bmh->palette[i][3] = 0;
	for (j = 0; j < BWCOLS; j++)
	{
		value = 255 - 255 * j / (BWCOLS - 1);
		// Shades of gray:
		bmh->palette[j][0] = bmh->palette[j][1] = bmh->palette[j][2] = (unsigned char)value;
		// Shades of red:
		bmh->palette[BWCOLS + j][0] = 0;
		bmh->palette[BWCOLS + j][1] = 0;
		bmh->palette[BWCOLS + j][2] = (unsigned char)value;
		// Shades of blue:
		bmh->palette[2 * BWCOLS + j][0] = (unsigned char)value;
		bmh->palette[2 * BWCOLS + j][1] = 0;
		bmh->palette[2 * BWCOLS + j][2] = 0;
		// Shades of green:
		bmh->palette[3 * BWCOLS + j][0] = 0;
		bmh->palette[3 * BWCOLS + j][1] = (unsigned char)value;
		bmh->palette[3 * BWCOLS + j][2] = 0;
	}
	bmPopulation = (int32_t*)Memory::xcalloc(bmh->imagesize, sizeof(int32_t));
	bmInfected = (int32_t*)Memory::xcalloc(bmh->imagesize, sizeof(int32_t));
	bmRecovered = (int32_t*)Memory::xcalloc(bmh->imagesize, sizeof(int32_t));
	bmTreated = (int32_t*)Memory::xcalloc(bmh->imagesize, sizeof(int32_t));

#ifdef _WIN32
	if (P.BitmapFormat == BitmapFormats::PNG)
	{
	  bmpdib = CreateDIBSection(GetDC(NULL), (BITMAPINFO*)bmp, DIB_RGB_COLORS, (void**)&bmPixels, NULL, NULL);
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

void Bitmap_Finalise()
{
#ifdef _WIN32
  if (P.BitmapFormat == BitmapFormats::PNG)
  {
    Gdiplus::GdiplusShutdown(m_gdiplusToken);
  }
#endif
}
