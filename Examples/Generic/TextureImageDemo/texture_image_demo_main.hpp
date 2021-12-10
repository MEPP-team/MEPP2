// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#include "Examples/Generic/TextureImageDemo/texture_image_demo_filter.hpp"

/**
 * \brief A mesh type templated main(argc, argv) function that
 *         - loads a mesh from a file,
 *         - applies the \ref texture_image_demo_filter generic filter,
 *         - write the resulting mesh to a file
 */
template< typename MeshT >
int
texture_image_demo_main(int argc, const char **argv)
{
  if(argc != 2)
  {
    std::cout << "Apply a dummy texture image filter to the input mesh."
              << std::endl;
    std::cout << "Usage:  " << argv[0] << "  input_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0]
              << "  ../Testing/Data/textures/face-multitexture/face.obj"
              << std::endl;
    return EXIT_FAILURE;
  }

  // input and output files
  std::string input_file_path = argv[1];
  std::string output_file_path = "texture_image_demo.output.obj";

  // read mesh from file
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

  // Note: the property maps must be extracted from the
  //       property maps bag, and explicitely passed as
  //       parameters to the filter, in order to make
  //       clear what property is used by the filter

  // retrieve the materials property map

  if(! has_map(pmaps_bag, FEVV::mesh_materials))
  {
    std::cout << "The mesh has no materials. Aborting." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "retrieve existing mesh_materials property map" << std::endl;
  auto mtl_pm = get_property_map(FEVV::mesh_materials, m, pmaps_bag);

  // apply filter
  texture_image_demo_filter(m, mtl_pm);

  // write mesh to file
  FEVV::Filters::write_mesh(output_file_path, m, pmaps_bag);

  return 0;
}
