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
#include <iostream>

#ifdef FEVV_USE_JPEG
#define cimg_use_jpeg // pb avec vcpkg en Release uniquement, why ?
#endif

#ifdef FEVV_USE_PNG
#define cimg_use_png
#endif

#ifdef FEVV_USE_TIFF
#define cimg_use_tiff
#endif

#undef _PTHREAD_H      // to avoid linking with pthread
#define cimg_display 0 // no display
#include <CImg.h>

int main(int argc, char **argv)
{
  if(argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " <path_to_image_file_1> [<path_to_image_file_2>...]\n";
    return -1; // error
  }

  for(int i = 1; i < argc; i++)
  {
    std::string filename(argv[i]);
    if(filename.substr(filename.size() - 4) == ".jpg")
    {
      // skip jpeg because of a bug in vcpkg jpeg library.
      // TODO: restore when vcpkg jpeg library is fixed.
      continue;
    }

    std::cout << "----------------------------\n";
    std::cout << "loading " << argv[i] << "...\n";
    cimg_library::CImg<unsigned char> image(argv[i]);

    //image.display();
    std::cout << "with     " << image.width()    << '\n';
    std::cout << "height   " << image.height()   << '\n';
    std::cout << "depth    " << image.depth()    << '\n';
    std::cout << "spectrum " << image.spectrum() << '\n';
  }

  return 0; // ok
}
