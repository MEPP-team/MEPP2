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

#include "FEVV/Filters/Generic/generic_writer.hpp"

#include "FEVV/Filters/CGAL/Progressive_Compression/progressive_decompression_filter.hpp"

/**
 * \brief A mesh type templated main(argc, argv) function that
 *         - applies the \ref progressive_decompression_filter generic filter
 *           on a compressed mesh file (binary file) ,
 *         - and write the resulting mesh to a file.
 */
template< typename MeshT >
int
progressive_decompression_main(int argc, const char **argv)
{
  if(argc < 2 || argc > 4)
  {
    std::cout << "Open a binary compressed 3D mesh " << std::endl;
    std::cout << "Usage:  " << argv[0] << "  input_binary_file_name  [output_mesh_filename [nb_max_batches]]" << std::endl;
    std::cout << "Example:  " << argv[0]
              << "  binarymesh.bin" << std::endl;
    return EXIT_FAILURE;
  }
  int nb_max_batches = 10000;
  if (argc == 4){
    nb_max_batches = atoi(argv[3]);
  }
  // Input and output files.
  std::string input_file_path = argv[1];
  std::string output_file_path = (argc >= 3)?argv[2]:
                                  "progressive_decompression_filter_output.obj";

  // Read mesh from file.
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;

  // Note: the property maps must be extracted from the
  //       property maps bag, and explicitely passed as
  //       parameters to the filter, in order to make
  //       clear what property is used by the filter.

  // Retrieve or create vertex-color property map.
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
    // Store property map in property maps bag.
    put_property_map(FEVV::vertex_color, m, pmaps_bag, v_cm);
  }
  // Retrieve or create edge color property map.
  typename FEVV::PMap_traits<FEVV::edge_color_t, MeshT>::pmap_type e_cm;
	if (has_map(pmaps_bag, FEVV::edge_color))
	{
		std::cout << "using edge color property map" << std::endl;
		e_cm = get_property_map(FEVV::edge_color, m, pmaps_bag);
	}
	else
  {
    std::cout << "create vertex-color map" << std::endl;
    e_cm = make_property_map(FEVV::edge_color, m);
    // Store property map in property maps bag.
    put_property_map(FEVV::edge_color, m, pmaps_bag, e_cm);
  }

  // Retrieve or create vertex-normal property map.
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
    // Store property map in property maps bag.
    put_property_map(FEVV::vertex_normal, m, pmaps_bag, v_nm);
  }

  // Retrieve point property map (aka geometry).
  auto pm = get(boost::vertex_point, m);

  ////////////////////////////////////////////////
  // Decompression algorithm call (Algorithm 2) //
  ////////////////////////////////////////////////	
  FEVV::Filters::progressive_decompression_filter(m, 
                                                  pm, 
                                                  v_cm, 
                                                  //v_nm,
                                                  e_cm,
                                                  input_file_path,
                                                  true,
                                                  nb_max_batches);
  // Write mesh to file.
  FEVV::Filters::write_mesh(output_file_path, m, pmaps_bag);

  return 0;
}
