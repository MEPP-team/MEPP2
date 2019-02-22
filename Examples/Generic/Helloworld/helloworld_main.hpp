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
#pragma once

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#include "Examples/Generic/Helloworld/helloworld_filter.hpp"

/**
 * \brief A mesh type templated main(argc, argv) function that
 *         - loads a mesh from a file,
 *         - applies the \ref helloworld_filter generic filter,
 *         - write the resulting mesh to a file
 */
template< typename MeshT >
int
helloworld_main(int argc, const char **argv)
{
  if(argc != 2)
  {
    std::cout << "Apply a dummy filter to the input mesh." << std::endl;
    std::cout << "Usage:  " << argv[0] << "  input_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0]
              << "  ../Testing/Data/CubeNonTriangleFaces.off" << std::endl;
    return EXIT_FAILURE;
  }

  // input and output files
  std::string input_file_path = argv[1];
  std::string output_file_path = "helloWorldFilter.output.obj";

  // read mesh from file
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

  // Note: the property maps must be extracted from the
  //       property maps bag, and explicitely passed as
  //       parameters to the filter, in order to make
  //       clear what property is used by the filter

  // retrieve or create vertex-color property map
  using VertexColorMap =
      typename FEVV::PMap_traits< FEVV::vertex_color_t, MeshT >::pmap_type;
  VertexColorMap v_cm;
  if(has_map(pmaps_bag, FEVV::vertex_color))
  {
    std::cout << "use existing vertex-color map" << std::endl;
    v_cm = get_property_map(FEVV::vertex_color, m, pmaps_bag);
  }
  else
  {
    std::cout << "create vertex-color map" << std::endl;
    v_cm = make_property_map(FEVV::vertex_color, m);
    // store property map in property maps bag
    put_property_map(FEVV::vertex_color, m, pmaps_bag, v_cm);
  }

  // retrieve or create vertex-normal property map
  using VertexNormalMap =
      typename FEVV::PMap_traits< FEVV::vertex_normal_t, MeshT >::pmap_type;
  VertexNormalMap v_nm;
  if(has_map(pmaps_bag, FEVV::vertex_normal))
  {
    std::cout << "use existing vertex-normal map" << std::endl;
    v_nm = get_property_map(FEVV::vertex_normal, m, pmaps_bag);
  }
  else
  {
    std::cout << "create vertex-normal map" << std::endl;
    v_nm = make_property_map(FEVV::vertex_normal, m);
    // store property map in property maps bag
    put_property_map(FEVV::vertex_normal, m, pmaps_bag, v_nm);
  }

  // retrieve point property map (aka geometry)
  auto pm = get(boost::vertex_point, m);

  // apply filter
  helloworld_filter(m, pm, v_cm, v_nm);

  // write mesh to file
  FEVV::Filters::write_mesh(output_file_path, m, pmaps_bag);

  return 0;
}
