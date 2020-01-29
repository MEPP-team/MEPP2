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

#include "FEVV/Filters/Generic/PointCloud/point_cloud_curvature.hpp"

#include "Testing/Utils/utils_identical_text_based_files.hpp"

#include <string>

/**
 * \brief  Generic test of Point Cloud Curvature filter.
 */
template< typename PointCloud >
int
test_point_cloud_curvature(int argc, const char **argv)
{
  // parse arguments

  if(argc < 3)
  {
    std::cout << "Load a point cloud then apply the Point Cloud Curvature"
              << " filter."
              << std::endl;
    std::cout << "Usage: " << argv[0]
              << "  input_point_cloud_file  output_file.ply"
              << "  [nearest_neighbors_number  [reference_file]]"
              << std::endl;
    std::cout << "Examples: \n"
              << argv[0] << "  casting.xyz  casting.curvature.ply\n"
              << argv[0] << "  casting.xyz  casting.curvature.ply  15\n"
              << argv[0] << "  casting.xyz  casting.curvature.ply  15"
                                "  casting.curvature.ref.ply\n"
              << std::endl;

    return EXIT_FAILURE;
  }

  std::string input_file_path = argv[1];
  std::string output_file_path = argv[2];
  std::string reference_file_path;
  unsigned int k = 15;
  if(argc > 3)
    k = static_cast<unsigned int>(std::stoi(argv[3]));
  if(argc > 4)
    reference_file_path = argv[4];

  // display parameters summary
  std::cout << "\nParameters summary:" << std::endl;
  std::cout << " - input file: " << input_file_path << std::endl;
  std::cout << " - output file: " << output_file_path << std::endl;
  std::cout << " - k: " << k << std::endl;

  //--------------------------------------------------

  // read point cloud from file
  PointCloud pc;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file_path, pc, pmaps_bag);

  // retrieve geometry property map
  auto pm = get(boost::vertex_point, pc);

  // create a vertex curvature property map
  auto v_curvm = FEVV::make_vertex_property_map< PointCloud, double >(pc);

  // create a vertex color property map
  auto v_cm = make_property_map(FEVV::vertex_color, pc);
  put_property_map(FEVV::vertex_color, pc, pmaps_bag, v_cm);

  // apply Point Cloud Curvature filter
  FEVV::Filters::point_cloud_curvature(pc, pm, k, v_curvm, v_cm);

  // save the point cloud
  std::cout << "saving point cloud with curvature to " << output_file_path
            << std::endl;
  FEVV::Filters::write_mesh(output_file_path, pc, pmaps_bag);

  // check output file
  if(!reference_file_path.empty())
  {
    std::cout << "Comparing output file" << std::endl;
    std::cout << "  '" << output_file_path << "'" << std::endl;
    std::cout << "with reference file" << std::endl;
    std::cout << "  '" << reference_file_path << "'" << std::endl;
    std::cout << "..." << std::endl;

    // use text file comparator
    if(!identical_text_based_files(output_file_path, reference_file_path))
    {
      std::cout << "Files are different!" << std::endl;
      return EXIT_FAILURE;
    }

    std::cout << "Files are identical." << std::endl;
  }

  return 0;
}

