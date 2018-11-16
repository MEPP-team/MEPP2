#include <stdio.h>

#ifdef FEVV_USE_JPEG
#define cimg_use_jpeg
#endif

#ifdef FEVV_USE_PNG
#define cimg_use_png
#endif

#ifdef FEVV_USE_TIFF
#define cimg_use_tiff
#endif

#include "CImg.h"
using namespace cimg_library;

int main(int argc, char **argv)
{
  if(argc < 4)
  {
    printf("Usage: ./cimg_tests img.jpg img.png img.tif\n");
    exit(EXIT_FAILURE);
  }

#ifdef FEVV_USE_JPEG
  CImg<unsigned char> image_j(argv[1]);
  //image_j.display("image_jpeg"); // no display under CIs !
#endif

#ifdef FEVV_USE_PNG
  CImg<unsigned char> image_p(argv[2]);
  //image_p.display("image_png"); // no display under CIs !
#endif

#ifdef FEVV_USE_TIFF
  CImg<unsigned char> image_t(argv[3]);
  //image_t.display("image_tiff"); // no display under CIs !
#endif

  return 0;
}
