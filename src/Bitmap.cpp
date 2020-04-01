#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "binio.h"
#include "Bitmap.h"
#include "Error.h"
#include "Param.h"
#include "Model.h"

//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
//// **** BITMAP stuff. 

#ifdef WIN32_BM
//HAVI avi;
ULONG_PTR m_gdiplusToken;
static HBITMAP bmpdib;
static CLSID  encoderClsid;
#else
#include <sys/stat.h> // for mkdir
#endif

static unsigned char* bmf, *bmPixels, *bmp;
// externs from SpatialSim.c
// TODO: move these to a header files
extern char OutFile[1024], OutFileBase[1024];

void CaptureBitmap(int ns, int tp)
{
	int i, x, y, f, mi;
	unsigned j;
	static double mx;
	static int fst = 1;
	double prev;

	mi = (int)(P.bwidth * P.bheight);
	if (fst)
	{
		fst = 0;
		mx = 0;
		for (i = 0; i < mi; i++) bmi[i] = 0;
		for (i = 0; i < P.N; i++)
		{
			x = ((int)(Households[Hosts[i].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[i].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				j = y * bmh->width + x;
				if ((j < bmh->imagesize) && (j >= 0))
				{
					bmi[j]++;
					if (bmi[j] > mx) mx = (double)bmi[j];
				}
			}
		}
		mx = log(1.001 * mx);
		for (i = 0; i < P.NMC; i++)
			if (Mcells[i].n > 0)
			{
				f = 0;
				if ((i < P.NMC - 1) && (i / P.nmch == (i + 1) / P.nmch) && (Mcells[i + 1].n > 0) && ((Mcells[i].country != Mcells[i + 1].country)
					|| ((P.DoAdunitBoundaryOutput) && ((AdUnits[Mcells[i].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor != (AdUnits[Mcells[i + 1].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor)))) f = 1;
				if ((i > 0) && (i / P.nmch == (i - 1) / P.nmch) && (Mcells[i - 1].n > 0) && (Mcells[i].country != Mcells[i - 1].country)) f = 1;
				if ((i < P.NMC - P.nmch) && (Mcells[i + P.nmch].n > 0) && ((Mcells[i].country != Mcells[i + P.nmch].country)
					|| ((P.DoAdunitBoundaryOutput) && ((AdUnits[Mcells[i].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor != (AdUnits[Mcells[i + P.nmch].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor)))) f = 1;
				if ((i >= P.nmch) && (Mcells[i - P.nmch].n > 0) && (Mcells[i].country != Mcells[i - P.nmch].country)) f = 1;
				if (f)
				{
					x = (int)(P.mcwidth * (((double)(i / P.nmch)) + 0.5) * P.scalex) - P.bminx;
					y = (int)(P.mcheight * (((double)(i % P.nmch)) + 0.5) * P.scaley) - P.bminy;
					if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
					{
						j = y * bmh->width + x;
						if ((j < bmh->imagesize) && (j >= 0)) bmi[j] = -1;
					}
				}
			}
		for (i = 0; i < P.bwidth / 2; i++)
		{
			prev = floor(3.99999 * ((double)i) * BWCOLS / ((double)P.bwidth) * 2);
			f = ((int)prev);
			for (j = 0; j < 10; j++)
			{
				bmPixels[(j + P.bheight + 5) * bmh->width + P.bwidth / 4 + i] = f;
			}
		}
	}
#pragma omp parallel for private(i) schedule(static,5000)
	for (i = 0; i < mi; i++)
	{
		if (bmi[i] == -1)
			bmPixels[i] = BWCOLS - 1; /* black for country boundary */
		else if ((bmi4[i] == 0) && (bmi2[i] == 0) && (bmi3[i] > 0))
			bmPixels[i] = (unsigned char)(3 * BWCOLS + BWCOLS * log(bmi3[i]) / mx);  /* green for recovered */
		else if (bmi2[i] > 0)
			bmPixels[i] = (unsigned char)(BWCOLS + BWCOLS * log(bmi2[i]) / mx); /* red for infected */
		else if (bmi4[i] > 0)
			bmPixels[i] = (unsigned char)(2 * BWCOLS + BWCOLS * log(bmi[i]) / mx); /* blue for treated */
		else if (bmi[i] > 0)
			bmPixels[i] = (unsigned char)(BWCOLS * log(bmi[i]) / mx); /* grey for just people */
		else
			bmPixels[i] = 0;
	}
}

void OutputBitmap(double t2, int tp)
{
	char buf[3000], OutF[3000];
	int j = 0;
	static int cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0;

	char *OutBaseName = strrchr(OutFile, '/');
	char *OutBaseName2 = strrchr(OutFile, '\\');
	if (OutBaseName2 != nullptr && (OutBaseName == nullptr || OutBaseName2 > OutBaseName)) {
		OutBaseName = OutBaseName2;
	}
	if (OutBaseName == nullptr) {
		OutBaseName = OutFile;
	}

	if (tp == 0)
	{
		j = cn1;
		cn1++;
		sprintf(OutF, "%s.ge" DIRECTORY_SEPARATOR "%s", OutFile, OutBaseName);
	}
	else if (tp == 1)
	{
		j = cn2;
		cn2++;
		sprintf(OutF, "%s.ge" DIRECTORY_SEPARATOR "Mean.%s", OutFile, OutBaseName);
	}
	else if (tp == 2)
	{
		j = cn3;
		cn3++;
		sprintf(OutF, "%s.ge" DIRECTORY_SEPARATOR "Min.%s", OutFile, OutBaseName);
	}
	else if (tp == 3)
	{
		j = cn4;
		cn4++;
		sprintf(OutF, "%s.ge" DIRECTORY_SEPARATOR "Max.%s", OutFile, OutBaseName);
	}

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
#elif defined(WIN32_BM)	
	//Windows specific bitmap manipulation code - could be recoded using LIBGD or another unix graphics library
	using namespace Gdiplus;

	wchar_t wbuf[1024];
	size_t a;

	//Add new bitmap to AVI
	//if ((P.OutputBitmap == 1) && (tp == 0)) AddAviFrame(avi, bmpdib, (unsigned char*)(&bmh->palette[0][0]));

	//This transfers HBITMAP to GDI+ Bitmap object
	Bitmap* gdip_bmp = Bitmap::FromHBITMAP(bmpdib, NULL);
	//Now change White in palette (first entry) to be transparent
	if ((cn1 == 1) && (tp == 0))
	{
		static UINT palsize;
		static ColorPalette* palette;
		palsize = gdip_bmp->GetPaletteSize();
		palette = (ColorPalette*)malloc(palsize);
		if (!palette) ERR_CRITICAL("Unable to allocate palette memory\n");
		(void)gdip_bmp->GetPalette(palette, palsize);
		palette->Flags = PaletteFlagsHasAlpha;
		palette->Entries[0] = 0x00ffffff; // Transparent white 
		gdip_bmp->SetPalette(palette);
	}
	//Now save as png
	sprintf(buf, "%s.%05i.png", OutF, j + 1); //sprintf(buf,"%s.ge" DIRECTORY_SEPARATOR "%s.%05i.png",OutFileBase,OutF,j+1);
	mbstowcs_s(&a, wbuf, strlen(buf) + 1, buf, _TRUNCATE);
	gdip_bmp->Save(wbuf, &encoderClsid, NULL);
	delete gdip_bmp;
#else
	sprintf(buf, "%s.%05i.bmp", OutF, j);
	FILE* dat;
	if (!(dat = fopen(buf, "wb"))) {
		char *errMsg = strerror(errno);
		if (errMsg == nullptr) {
			ERR_CRITICAL("strerror failed.\n");
		}
		ERR_CRITICAL_FMT("Unable to open bitmap file %s (%d): %s\n", buf, errno, errMsg);
	}
	fprintf(dat, "BM");
	fwrite_big((void*)bmf, sizeof(unsigned char), sizeof(bitmap_header) / sizeof(unsigned char) + bmh->imagesize, dat);
	fclose(dat);
#endif
}
void InitBMHead()
{
	int i, j, k, k2, value;

	fprintf(stderr, "Initialising bitmap\n");
	k = P.bwidth * P.bheight2;
	k2 = sizeof(bitmap_header) / sizeof(unsigned char);

	if (!(bmf = (unsigned char*)calloc((size_t)k + k2, sizeof(unsigned char))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	bmPixels = &(bmf[k2]);
	bmp = &(bmf[12]);
	bmh = (bitmap_header*)bmf;
	bmh->spare = 0;
	bmh->boffset = 2 + sizeof(bitmap_header);
	bmh->headersize = 40; // BITMAPINFOHEADER
	bmh->width = P.bwidth;
	bmh->height = P.bheight2;
	bmh->PlanesAndBitspp = 1 // Number of colour planes; must be 1
	                     + (8 << 16); // Colour depth: 8 bits per pixel
	bmh->compr = 0; // No compression (BI_RGB)
	bmh->imagesize = bmh->width * bmh->height;
	bmh->filesize = 2 // "BM"
	              + ((unsigned int) sizeof(bitmap_header)) // BITMAP_HEADER
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
	if (!(bmi = (float*)malloc(bmh->imagesize * sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	if (!(bmi2 = (float*)malloc(bmh->imagesize * sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	if (!(bmi3 = (float*)malloc(bmh->imagesize * sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	if (!(bmi4 = (float*)malloc(bmh->imagesize * sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");

#ifdef WIN32_BM
	bmpdib = CreateDIBSection(GetDC(NULL), (BITMAPINFO*)bmp, DIB_RGB_COLORS, (void**)& bmPixels, NULL, NULL);
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	Gdiplus::GdiplusStartup(&m_gdiplusToken, &gdiplusStartupInput, NULL);

	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes

	Gdiplus::ImageCodecInfo* pImageCodecInfo = NULL;
	Gdiplus::GetImageEncodersSize(&num, &size);
	if (!(pImageCodecInfo = (Gdiplus::ImageCodecInfo*)(malloc(size))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	Gdiplus::GetImageEncoders(num, size, pImageCodecInfo);
	for (UINT j = 0; j < num; ++j)
	{
// Visual Studio Analyze incorrectly reports this because it doesn't understand Gdiplus::GetImageEncodersSize()
// warning C6385: Reading invalid data from 'pImageCodecInfo':  the readable size is 'size' bytes, but '208' bytes may be read.
#pragma warning( suppress: 6385 )
		const WCHAR* type = pImageCodecInfo[j].MimeType;
		if (wcscmp(type, L"image/png") == 0)
		{
			encoderClsid = pImageCodecInfo[j].Clsid;
			j = num;
		}
	}
	free(pImageCodecInfo);
#endif

	char buf[1024+3];
	sprintf(buf, "%s.ge", OutFileBase);
#ifdef WIN32_BM
	if (!(CreateDirectory(buf, NULL))) fprintf(stderr, "Unable to create directory %s\n", buf);
#else
	if (!(mkdir(buf, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))) fprintf(stderr, "Unable to create directory %s\n", buf);
#endif
}

