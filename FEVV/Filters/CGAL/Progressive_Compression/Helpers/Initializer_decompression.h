#pragma once

//#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"

#include "FEVV/Filters/CGAL/Surface_mesh/Progressive_Compression_ld/progressive_decompression_filter.hpp"

/**
 * \brief 
 */
template< typename MeshT >
int
progressive_decompression_main(int argc, const char **argv)
{
  if(argc < 2 || argc > 3)
  {
    std::cout << "Open a binary compressed 3D mesh " << std::endl;
    std::cout << "Usage:  " << argv[0] << "  input_mesh_filename  [output_mesh_filename]" << std::endl;
    std::cout << "Example:  " << argv[0]
              << "  binarymesh.bin" << std::endl;
    return EXIT_FAILURE;
  }

  // input and output files
  std::string input_file_path = argv[1];
  std::string output_file_path = (argc == 3)?argv[2]:
                                  "progressive_decompressionFilter.output.obj";

  // read mesh from file
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  //FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

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
  // retrieve or create edge color property map
  typename FEVV::PMap_traits<FEVV::edge_color_t, MeshT>::pmap_type e_cm;
	if (has_map(pmaps_bag, FEVV::edge_color))
	{
		std::cout << "using edge color property map" << std::endl;
		e_cm = get_property_map(FEVV::edge_color, m, pmaps_bag);
		//e_cm_ptr = &e_cm;
	}
	else
  {
    std::cout << "create vertex-color map" << std::endl;
    e_cm = make_property_map(FEVV::edge_color, m);
    // store property map in property maps bag
    put_property_map(FEVV::edge_color, m, pmaps_bag, e_cm);
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
  FEVV::Filters::progressive_decompression_filter(m, 
                                                  pm, 
                                                  v_cm, 
                                                  //v_nm,
                                                  e_cm,
                                                  input_file_path,
                                                  true);
  // write mesh to file
  FEVV::Filters::write_mesh(output_file_path, m, pmaps_bag);

  return 0;
}
