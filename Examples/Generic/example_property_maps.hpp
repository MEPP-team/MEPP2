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


// print property map content
template< typename PointCloud, typename VertexPropertyMap >
void
print_property_map(const PointCloud &pc,
                   const VertexPropertyMap &pmap,
                   const std::string &property_map_name)
{
  std::cout << "property map " << property_map_name << ":" << std::endl;

  int counter = 0;
  auto iterator_pair = vertices(pc);
  auto vi = iterator_pair.first;
  auto vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    counter++;
    auto value = get(pmap, *vi);
    std::cout << "  vertex #" << counter << ": " << value << std::endl;
  }
}


// example of property maps use
template< typename PointCloudT >
int
example_property_maps(int argc, char *argv[])
{
  // parse arguments
  if(argc != 3)
  {
    std::cout << "Usage:  " << argv[0]
              << "  input_mesh_filename  output_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0]
              << "  ../Testing/Data/tetra.xyz  tetra.out.ply" << std::endl;
    return EXIT_FAILURE;
  }
  std::string input_file = argv[1];
  std::string output_file = argv[2];

  //----------------------------------

  // create point cloud
  PointCloudT pc;

  // load point cloud
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file, pc, pmaps_bag);

  // print geometry
  print_property_map(pc, get(boost::vertex_point, pc), "PointMap");
  
  //----------------------------------

  // create vertex-normal property map
  using VertexNormalMap = typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                                      PointCloudT >::pmap_type;
  VertexNormalMap v_nm;
  if(has_map(pmaps_bag, FEVV::vertex_normal))
  {
    std::cout << "use existing vertex-normal map" << std::endl;
    v_nm = get_property_map(FEVV::vertex_normal, pc, pmaps_bag);
    print_property_map(pc, v_nm, "VertexNormalMap");
  }
  else
  {
    std::cout << "create vertex-normal map" << std::endl;
    v_nm = make_property_map(FEVV::vertex_normal, pc);
    // store property map in property maps bag
    put_property_map(FEVV::vertex_normal, pc, pmaps_bag, v_nm);
  }

  // populate vertex normal map
  typedef typename boost::property_traits< VertexNormalMap >::value_type Normal;

  std::cout << "populate vertex-normal map" << std::endl;
  int counter = 0;
  auto iterator_pair = vertices(pc);
  auto vi = iterator_pair.first;
  auto vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    counter++;
    std::cout << "  processing vertex #" << counter << std::endl;

    // write property map data
    Normal new_normal(42.1, 42.2, 42.3);
    put(v_nm, *vi, new_normal);
  }

  //----------------------------------

  // retrieve or create vertex-color property map
  using VertexColorMap = typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                                     PointCloudT >::pmap_type;
  VertexColorMap v_cm;
  if(has_map(pmaps_bag, FEVV::vertex_color))
  {
    std::cout << "use existing vertex-color map" << std::endl;
    v_cm = get_property_map(FEVV::vertex_color, pc, pmaps_bag);
    print_property_map(pc, v_cm, "VertexColorMap");
  }
  else
  {
    std::cout << "create vertex-color map" << std::endl;
    v_cm = make_property_map(FEVV::vertex_color, pc);
    // store property map in property maps bag
    put_property_map(FEVV::vertex_color, pc, pmaps_bag, v_cm);
  }

  // populate vertex color map
  typedef typename boost::property_traits< VertexColorMap >::value_type Color;
  typedef boost::graph_traits< PointCloudT > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;

  std::cout << "populate vertex-color map" << std::endl;
  counter = 0;
  iterator_pair = vertices(pc);
  vi = iterator_pair.first;
  vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    counter++;
    std::cout << "  processing vertex #" << counter << std::endl;

    // write property map data
    Color newcolor(8, 2, 1);
    put(v_cm, *vi, newcolor);
  }

  //----------------------------------

  // create 2 non-standard prop maps
  auto vertex_int_map = FEVV::make_vertex_property_map< PointCloudT, int >(pc);
  auto vertex_color2_map =
      FEVV::make_vertex_property_map< PointCloudT, Color >(pc);

  // populate non-standard prop maps
  std::cout << "populate vertex_int_map and vertex_color2_map (non-standard "
               "prop maps)"
            << std::endl;
  counter = 0;
  iterator_pair = vertices(pc);
  vi = iterator_pair.first;
  vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    counter++;
    std::cout << "  processing vertex #" << counter << std::endl;

    // write property map data
    Color newcolor((counter*3) % 255 , (counter*3 + 1) % 255, (counter*3 + 2) % 255);
    put(vertex_color2_map, *vi, newcolor);

    put(vertex_int_map, *vi, 10 + counter);
  }

  // print non-standard property map content
  print_property_map(pc, vertex_int_map, "non-standard vertex_int_map");
  print_property_map(pc, vertex_color2_map, "non-standard vertex_color2_map");

  //----------------------------------

  // save point cloud with normal and color
  FEVV::Filters::write_mesh(output_file, pc, pmaps_bag);
  
  return 0;
}
