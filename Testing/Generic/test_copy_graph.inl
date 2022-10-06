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
#include <vector>


template< typename MeshT >
int
test_copy_graph(int argc, const char **argv)
{
  // parse arguments

  if(argc < 4)
  {
    std::cout << "Load a mesh, copy it in a new mesh, then save the new mesh."
              << std::endl
              << "Usage: " << argv[0]
              << "  input_mesh_file  output_mesh_file  reference_mesh_file"
              << std::endl
              << "Example: " << argv[0]
              << "  airplane.off  airplane_copy.off  airplane_copy.ref.off"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string input_file_path = argv[1];
  std::string output_file_path = argv[2];
  std::string reference_file_path = argv[3];

  // display parameters summary
  std::cout << "\nParameters summary:" << std::endl
            << " - input file: "       << input_file_path     << std::endl
            << " - output file: "      << output_file_path    << std::endl
            << " - reference file: "   << reference_file_path << std::endl
            << std::endl;

  //--------------------------------------------------

  // read mesh from file
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

  // copy mesh
  MeshT m2;
  FEVV::PMapsContainer pmaps_bag2;
  FEVV::Filters::copy_graph(m, pmaps_bag, m2, pmaps_bag2);

  // save the copied mesh
  FEVV::Filters::write_mesh(output_file_path, m2, pmaps_bag2);

  // check output file
  std::cout << "Comparing output file" << std::endl
            << "  '" << output_file_path << "'" << std::endl
            << "with reference file" << std::endl
            << "  '" << reference_file_path << "'" << std::endl
            << "..." << std::endl;

  if(FEVV::FileUtils::has_extension(output_file_path, ".off") ||
     FEVV::FileUtils::has_extension(output_file_path, ".coff"))
  {
    // use OFF file comparator
    if(! are_meshes_equal(output_file_path,
                          reference_file_path,
                          false))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    // use text file comparator
    std::vector< std::string > skip = { "mtllib" };
    if(! identical_text_based_files(output_file_path,
                                    reference_file_path,
                                    skip))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Files are identical." << std::endl;

  return 0;
}
