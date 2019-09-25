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

#include "FEVV/DataStructures/DataStructures_pcl_point_cloud.h"
#include "FEVV/Wrappings/Geometry_traits_pcl_point_cloud.h"
#include "FEVV/Wrappings/Graph_properties_pcl_point_cloud.h"
#include "FEVV/Wrappings/properties_pcl_point_cloud.h"

#include "FEVV/Filters/PCL/pcl_point_cloud_reader.hpp"
#include "FEVV/Filters/PCL/pcl_point_cloud_writer.hpp"

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "Testing/Utils/utils_identical_text_based_files.hpp"

#include <iostream>
#include <string>


// main
int main(int argc, char *argv[])
{
  // parse arguments
  if(argc != 4)
  {
    std::cout << "Usage:  " << argv[0]
              << "  input_mesh_filename  output_mesh_filename  reference_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0] << "  ../Testing/Data/tetra.xyz  tetra.out.ply  tetra.out.ref.ply"
              << std::endl;
    return EXIT_FAILURE;
  }
  std::string input_file = argv[1];
  std::string output_file = argv[2];
  std::string reference_file = argv[3];

  //----------------------------------

  // create point cloud
  typedef  FEVV::PCLPointCloud  PointCloudT;
  PointCloudT pc;

  // load point cloud
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file, pc, pmaps_bag);

  //----------------------------------

  // save point cloud
  FEVV::Filters::write_mesh(output_file, pc, pmaps_bag);

  //----------------------------------

  // check output file
  std::cout << "Comparing output file '" << output_file
            << "' with reference file '" << reference_file << "'..."
            << std::endl;

  if(FEVV::FileUtils::has_extension(output_file, ".off") ||
     FEVV::FileUtils::has_extension(output_file, ".coff"))
  {
    // use OFF file comparator
    if(!are_meshes_equal(output_file, reference_file, false))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    // use text file comparator
    if(!identical_text_based_files(output_file, reference_file))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Files are identical." << std::endl;
  
  return 0;
}
