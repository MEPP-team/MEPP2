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

#include "FEVV/DataStructures/DataStructures_cgal_point_set.h"
  // for FEVV::CGALPointSet
#include "FEVV/Wrappings/Graph_traits_cgal_point_set.h"
  // for boost::graph_traits< FEVV::CGALPointSet >
#include "FEVV/Wrappings/Geometry_traits_cgal_point_set.h"
  // for FEVV::RetrieveKernel< FEVV::CGALPointSet >
#include "FEVV/Wrappings/Graph_properties_cgal_point_set.h"
  // for get(FEVV::CGALPointSetPointMap&, ...)
#include "FEVV/Wrappings/properties_cgal_point_set.h"
  // for FEVV::PMap_traits< FEVV::CGALPointSet > and vector-property-maps

#include "FEVV/Filters/CGAL/CGALPointSet/cgal_point_set_reader.hpp"
  // for FEVV::Filters::read_mesh< FEVV::CGALPointSet >
#include "FEVV/Filters/CGAL/CGALPointSet/cgal_point_set_writer.hpp"
  // for FEVV::Filters::write_mesh< FEVV::CGALPointSet >


// print property map content
template< typename PointCloud, typename VertexPropertyMap >
void
print_property_map(const PointCloud &pc,
                   const VertexPropertyMap &pmap,
                   const std::string &property_name)
{
  int counter = 0;
  auto iterator_pair = vertices(pc);
  auto vi = iterator_pair.first;
  auto vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    counter++;
    auto value = get(pmap, *vi);

    std::cout << "vertex #" << counter
              << " has " << property_name
              << " value " << value
              << std::endl;
  }
}


// main
int main(int argc, char *argv[])
{
  // parse arguments
  if(argc != 3)
  {
    std::cout << "Usage:  " << argv[0]
              << "  input_mesh_filename  output_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0] << "  ../Testing/Data/tetra.xyz  tetra.out.ply"
              << std::endl;
    return EXIT_FAILURE;
  }
  std::string input_file = argv[1];
  std::string output_file = argv[2];

  //----------------------------------

  // create point cloud
  typedef  FEVV::CGALPointSet  PointCloudT;
  PointCloudT pc;

  // load point cloud
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file, pc, pmaps_bag);

  // print geometry
  print_property_map(pc, get(boost::vertex_point, pc), "coordinates");
  
  //----------------------------------

  // create vertex-normal property map
  using VertexNormalMap = typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                                      PointCloudT >::pmap_type;
  VertexNormalMap v_nm;
  if(has_map(pmaps_bag, FEVV::vertex_normal))
  {
    std::cout << "use existing vertex-normal map" << std::endl;
    v_nm = get_property_map(FEVV::vertex_normal, pc, pmaps_bag);
    print_property_map(pc, v_nm, "normal");
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
    print_property_map(pc, v_cm, "color");
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
    Color newcolor(0.8, 0.2, 0.1);
    put(v_cm, *vi, newcolor);
  }

  //----------------------------------

  // save point cloud with normal and color
  FEVV::Filters::write_mesh(output_file, pc, pmaps_bag);
  
  return 0;
}
