// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/Graph_traits_extension_cgal_polyhedron_3.h"

#include "FEVV/Wrappings/properties_polyhedron_3.h"
#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "Testing/Utils/utils_identical_text_based_files.hpp"

#include <iostream>
#include <string>


using CGALKernel = CGAL::Cartesian< double >;
using MeshT =
    CGAL::Polyhedron_3< CGALKernel, CGAL::Polyhedron_items_with_id_3 >;


int
main(int argc, const char **argv)
{
  if(argc != 4)
  {
    std::cout << "Load a mesh from an input file using the generic reader, "
                 "write it to an output file using the generic writer, then "
                 "compare the output file with a reference file."
              << std::endl;
    std::cout << "Usage: " << argv[0]
              << "  input_mesh_file  output_mesh_file  reference_mesh_file"
              << std::endl;
    std::cout << "Example: " << argv[0]
              << "  airplane.off airplane.out.obj airplane.ref.obj"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string input_file_path = argv[1];
  std::string output_file_path = argv[2];
  std::string reference_file_path = argv[3];

  //--------------------------------------------------

  // read mesh from file
  MeshT m;
  FEVV::PMapsContainer pmaps;
  FEVV::Filters::read_mesh(input_file_path, m, pmaps);

  // write mesh to file
  FEVV::Filters::write_mesh(output_file_path, m, pmaps);

  // check output file
  std::cout << "Comparing output file '" << output_file_path
            << "' with reference file '" << reference_file_path << "'..."
            << std::endl;

  if(FEVV::FileUtils::has_extension(output_file_path, ".off") ||
     FEVV::FileUtils::has_extension(output_file_path, ".coff"))
  {
    // use OFF file comparator
    if(!are_meshes_equal(output_file_path, reference_file_path, false))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    // use text file comparator
    if(!identical_text_based_files(output_file_path, reference_file_path))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Files are identical." << std::endl;

  return 0;
}
