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

#include "FEVV/DataStructures/DataStructures_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/properties_polyhedron_3.h"

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"

#include "FEVV/Filters/CGAL/Boolean_Operations/boolean_operations.hpp"

#include "Testing/Utils/utils_are_meshes_identical.hpp"


int usage_and_exit(const char **argv)
{
  std::cout << "Apply boolean union, intersection or subtraction"
               " on two meshes."
            << std::endl;
  std::cout << "Usage:  "
            << argv[0]
            << "  union|inter|minus  mesh_file_1  mesh_file_2"
               "  [reference_mesh_file]"
            << std::endl;
  std::cout << "Example:  "
            << argv[0]
            << "  union"
               "  ../Testing/Data/CubeTriangleFaces.off"
               "  ../Testing/Data/tetra.off"
               "  ../Testing/Data/boolean_operations/mepp1_cube_tetra_union.off"
            << std::endl;

  return EXIT_FAILURE;
}


int main(int argc, const char **argv)
{
  // check arguments
  if(argc < 4  ||  argc > 5)
    return usage_and_exit(argv);

  // check operation
  std::string operation = argv[1];
  if(operation != "union" && operation != "inter" && operation != "minus")
    return usage_and_exit(argv);

  // input and output files
  std::string input_file_1 = argv[2];
  std::string input_file_2 = argv[3];
  std::string reference_file;
  if(argc == 5)
    reference_file = argv[4];
  std::string output_file =
      "test_boolean_operations_polyhedron_" + operation + ".off";

  // display parameters summary
  std::cout << "\nParameters summary:" << std::endl;
  std::cout << " - operation     : " << operation << std::endl;
  std::cout << " - input file 1  : " << input_file_1 << std::endl;
  std::cout << " - input file 2  : " << input_file_2 << std::endl;
  std::cout << " - reference file: " << reference_file << std::endl;
  std::cout << "\nOutput in " << output_file << std::endl;

  // read input meshes from files
  FEVV::MeshPolyhedron m1;
  FEVV::PMapsContainer pmaps_bag_1;
  FEVV::Filters::read_mesh(input_file_1, m1, pmaps_bag_1);
  FEVV::MeshPolyhedron m2;
  FEVV::PMapsContainer pmaps_bag_2;
  FEVV::Filters::read_mesh(input_file_2, m2, pmaps_bag_2);

  // create output mesh
  FEVV::MeshPolyhedron m_out;
  auto pm_out = get(boost::vertex_point, m_out);
  FEVV::PMapsContainer pmaps_bag_out;

  // apply filter, result in m_out
  std::cout << "Running boolean operation" + operation + "..." << std::endl;
  if(operation == "union")
    FEVV::Filters::boolean_union(m1, m2, m_out);
  else if(operation == "inter")
    FEVV::Filters::boolean_inter(m1, m2, m_out);
  else
    FEVV::Filters::boolean_minus(m1, m2, m_out);

  // write result to file
  FEVV::Filters::write_mesh(output_file, m_out, pmaps_bag_out);

  // compare output to reference
  if(reference_file.empty())
  {
    std::cout << "No reference file provided, skipping result check."
              << std::endl;
  }
  else
  {
    std::cout << "Checking result with reference file..." << std::endl;
    bool is_result_ok = are_meshes_equal(reference_file,
                                         output_file,
                                         false, /*verbose*/
                                         1e-5f,
                                         true); /*relative error*/
    if(is_result_ok)
    {
      std::cout << "Boolean operation " + operation + " successfully tested."
                << std::endl;
    }
    else
    {
      std::cout << "Boolean operation " + operation + " FAILED!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  return 0;
}
