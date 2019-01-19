// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
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
#endif

#ifdef FEVV_USE_PNG
  CImg<unsigned char> image_p(argv[2]);
#endif

#ifdef FEVV_USE_TIFF
  CImg<unsigned char> image_t(argv[3]);
#endif

  return 0;
}
