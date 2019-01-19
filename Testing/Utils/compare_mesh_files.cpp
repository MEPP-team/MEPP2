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

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "Testing/Utils/utils_identical_text_based_files.hpp"

#include <iostream>
#include <string>


int
main(int argc, const char **argv)
{
  if(argc < 3 || argc > 6)
  {
    std::cout << "Compare two mesh files." << std::endl;
    std::cout << "For .OFF files, the comparison can be exact or approximated. "
                 "In the approximated case, the threshold can be a relative "
                 "value (default) or an absolute value."
              << std::endl;
    std::cout << std::endl;
    std::cout << "Usage #1, exact comparison:" << std::endl;
    std::cout << "  " << argv[0] << "  mesh_file_1  mesh_file_2" << std::endl;
    std::cout << "Usage #2, approximate comparison, same threshold for both "
                 "geometry and attributes:"
              << std::endl;
    std::cout << "  " << argv[0]
              << "  mesh_file_1  mesh_file_2  threshold  [ABSOLUTE]"
              << std::endl;
    std::cout << "Usage #3, approximate comparison, separate thresholds for "
                 "geometry and attributes:"
              << std::endl;
    std::cout << "  " << argv[0]
              << "  mesh_file_1  mesh_file_2  geometry_threshold  "
                 "attributes_threshold  [ABSOLUTE]"
              << std::endl;

    std::cout << std::endl;
    std::cout << "Examples: " << std::endl;
    std::cout << "  " << argv[0] << "  airplane1.off  airplane2.off"
              << std::endl;
    std::cout << "  " << argv[0] << "  airplane1.off  airplane2.off  1e-3"
              << std::endl;
    std::cout << "  " << argv[0]
              << "  airplane1.off  airplane2.off  1e-3  ABSOLUTE" << std::endl;
    std::cout << "  " << argv[0] << "  airplane1.off  airplane2.off  1e-3  1e-4"
              << std::endl;
    std::cout << "  " << argv[0]
              << "  airplane1.off  airplane2.off  1e-3  1e-4  ABSOLUTE"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string mesh_file1 = argv[1];
  std::string mesh_file2 = argv[2];

  bool relative_threshold = true;
  if(std::string(argv[argc - 1]) == "ABSOLUTE")
    relative_threshold = false;

  float geom_threshold = 0;
  float attr_threshold = 0;
  if((argc == 4) && relative_threshold)
  {
    // relative mode, same threshold
    geom_threshold = (float)std::atof(argv[3]);
    attr_threshold = (float)std::atof(argv[3]);
  }
  else if((argc == 5) && relative_threshold)
  {
    // relative mode, different thresholds
    geom_threshold = (float)std::atof(argv[3]);
    attr_threshold = (float)std::atof(argv[4]);
  }
  else if((argc == 5) && !relative_threshold)
  {
    // absolute mode, same threshold
    geom_threshold = (float)std::atof(argv[3]);
    attr_threshold = (float)std::atof(argv[3]);
  }
  else if((argc == 6) && !relative_threshold)
  {
    // absolute mode, same threshold
    geom_threshold = (float)std::atof(argv[3]);
    attr_threshold = (float)std::atof(argv[4]);
  }

  //--------------------------------------------------

  std::cout << "Comparing file '" << mesh_file1 << "'" << std::endl;
  std::cout << "     with file '" << mesh_file2 << "'" << std::endl;

  if(FEVV::FileUtils::has_extension(mesh_file1, ".off") ||
     FEVV::FileUtils::has_extension(mesh_file2, ".coff"))
  {
    // use OFF file comparator
    std::cout << "Use OFF file comparison" << std::endl;

    std::cout << "Geometry tolerance: " << geom_threshold << std::endl;
    std::cout << "Attributes tolerance: " << attr_threshold << std::endl;
    if(relative_threshold)
      std::cout << "Mode: relative error" << std::endl;
    else
      std::cout << "Mode: absolute error" << std::endl;

    if(!are_meshes_equal(mesh_file1,
                         mesh_file2,
                         false,
                         geom_threshold,
                         attr_threshold,
                         relative_threshold))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    // use text file comparator
    std::cout << "Use text file comparison" << std::endl;

    if(!identical_text_based_files(mesh_file1, mesh_file2))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Files are identical." << std::endl;

  return 0;
}
