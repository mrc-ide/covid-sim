#ifndef COVIDSIM_BITMAP_H_INCLUDED_
#define COVIDSIM_BITMAP_H_INCLUDED_

#include <cstdint>
#include <string>

#ifdef UNIX
#define DIRECTORY_SEPARATOR "/"
#else
#define DIRECTORY_SEPARATOR "\\"
#endif

#include "geometry/Size.h"

const int BWCOLS = 58;

namespace CovidSim
{
  namespace BitMap
  {
    struct Header
    {
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
    };

    /// Enumeration of bitmap formats.
    enum struct Formats
    {
      PNG,  ///< default if IMAGE_MAGICK or _WIN32 defined
      BMP   ///< fall-back
    };

    struct Builder
    {
      /// Height in pixels of the entire bitmap output, including both the spectrum at the top and the map area
      int height2_;

      /// Size in pixels of the map area in the bitmap output
      Geometry::Size<int> bounds_;

      Geometry::Vector2<int> min_;

      ///< Format of bitmap (platform dependent and command-line /BM: specified).
      Formats format_;

      /// Whether to output a bitmap
      int output_;

      /// Number of pixels per degree in bitmap output
      Geometry::DiagonalMatrix2d scale_;

      // These allow up to about 2 billion people per pixel, which should be ample.
      /// The population in each bitmap pixel. Special value -1 means "country boundary"
      static int32_t *population_;

      /// The number of infected people in each bitmap pixel.
      static int32_t *infected_;

      /// The number of recovered people in each bitmap pixel.
      static int32_t *recovered_;

      /// The number of treated people in each bitmap pixel.
      static int32_t *treated_;

      Header *header_;

      void capture();
      void output(int, std::string const&);
      void initialise_header(std::string const&);

      void finalise();

      static unsigned char* bitmap_;
      static unsigned char* pixels_;
      static unsigned char* info_;
    };
  }
}

#endif // COVIDSIM_BITMAP_H_INCLUDED_
