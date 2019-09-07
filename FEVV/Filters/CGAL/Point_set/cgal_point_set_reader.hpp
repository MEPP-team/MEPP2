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

#include "FEVV/DataStructures/DataStructures_cgal_point_set.h"
#include "FEVV/Wrappings/Wrappings_cgal_point_set.h"

#include "FEVV/Wrappings/properties.h"

#include "FEVV/Filters/Generic/generic_reader.hpp" // specialization of
#include "FEVV/Tools/IO/FileUtilities.hpp" // for FileUtils::has_extension()

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/read_ply_points.h>
#include <stdexcept> // for std::invalid_argument
#include <string> // for std::getline


namespace FEVV {

// parse PLY header to detect normals and colors
inline
void
parse_ply_header(std::ifstream &in, bool &has_normal, bool &has_color)
{
  std::string line;

  in.seekg(0); // rewind
  has_normal = false;
  has_color = false;

  do
  {
    std::getline(in, line);

    // look for "property ... nx" and "property ... red"
    if(line.substr(0, 8) == "property")
    {
      if(line.substr(line.size()-2) == "nx")
      {
        has_normal = true;
        //DBG std::cout << "normal detected in PLY header" << std::endl;
      }
      else if(line.substr(line.size()-3) == "red")
      {
        has_color = true;
        //DBG std::cout << "color detected in PLY header" << std::endl;
      }
    }
  }
  while(line != "end_header");
  
  in.seekg(0); // rewind
}

} // namespace FEVV


namespace FEVV {
namespace Filters {


/**
 * \brief  Load mesh from file.
 *
 * \param  filename    name of the input mesh file
 * \param  g           mesh object to store the mesh being read
 * \param  pmaps       property maps bag to store the properties
 *                     of the mesh being read
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template<>
inline
void
read_mesh< FEVV::CGALPointSet, FEVV::Geometry_traits< FEVV::CGALPointSet > >(
    const std::string &filename,
    FEVV::CGALPointSet &g,
    PMapsContainer &pmaps,
    bool only_pts)
{
#if 1
  bool success = false;

  // open input file
  std::ifstream in(filename);
  if(! in)
  {
    throw std::invalid_argument("read_mesh() error: can not open file " +
                                filename);
  }

  // load point cloud
  if(FEVV::FileUtils::has_extension(filename, ".xyz"))
  {
    // load geometry + normal
    success = CGAL::read_xyz_point_set(in, g);

    if(g.has_normal_map())
      put_property_map(FEVV::vertex_normal, g, pmaps, g.normal_map());
  }
  else if(FEVV::FileUtils::has_extension(filename, ".off"))
  {
    // load geometry + normal
    success = CGAL::read_off_point_set(in, g);

    if(g.has_normal_map())
      put_property_map(FEVV::vertex_normal, g, pmaps, g.normal_map());
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
    // parse ply header to detect normals and colors
    bool has_normal;
    bool has_color;
    parse_ply_header(in, has_normal, has_color);

    // create needed property maps
    using VertexNormalMap =
        typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                    FEVV::CGALPointSet >::pmap_type;
    using VertexColorMap =
        typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                    FEVV::CGALPointSet >::pmap_type;
    VertexNormalMap v_nm;
    VertexColorMap v_cm;

    if(has_normal)
    {
      v_nm = make_property_map(FEVV::vertex_normal, g);
      put_property_map(FEVV::vertex_normal, g, pmaps, v_nm);
    }

    if(has_color)
    {
      v_cm = make_property_map(FEVV::vertex_color, g);
      put_property_map(FEVV::vertex_color, g, pmaps, v_cm);
    }

    if(has_normal && has_color)
    {
      success = CGAL::read_ply_points_with_properties(
          in,
          g.index_back_inserter(),
          CGAL::make_ply_point_reader(g.point_push_map()),
          CGAL::make_ply_normal_reader(g.push_property_map(v_nm)),
          std::make_tuple(g.push_property_map(v_cm),
                          typename CGAL::Kernel_traits<
                              typename VertexColorMap::value_type >::Kernel::
                              Construct_vector_3(),
                          CGAL::PLY_property< double >("red"),
                          CGAL::PLY_property< double >("green"),
                          CGAL::PLY_property< double >("blue")));
    }
    else if(has_normal)
    {
      success = CGAL::read_ply_points_with_properties(
          in,
          g.index_back_inserter(),
          CGAL::make_ply_point_reader(g.point_push_map()),
          CGAL::make_ply_normal_reader(g.push_property_map(v_nm)));
    }
    else if(has_color)
    {
      success = CGAL::read_ply_points_with_properties(
          in,
          g.index_back_inserter(),
          CGAL::make_ply_point_reader(g.point_push_map()),
          std::make_tuple(g.push_property_map(v_cm),
                          typename CGAL::Kernel_traits<
                              typename VertexColorMap::value_type >::Kernel::
                              Construct_vector_3(),
                          CGAL::PLY_property< double >("red"),
                          CGAL::PLY_property< double >("green"),
                          CGAL::PLY_property< double >("blue")));
    }
    else
    {
      success = CGAL::read_ply_points_with_properties(
          in,
          g.index_back_inserter(),
          CGAL::make_ply_point_reader(g.point_push_map()));
    }
  }

  if(! success)
  {
    throw std::invalid_argument("read_mesh() error: can not read file " +
                                filename);
  }
#else
  // DBG
  g.push_back(FEVV::CGALPoint(-1, -1, 2));
  g.push_back(FEVV::CGALPoint(-1,  1, 2));
  g.push_back(FEVV::CGALPoint( 1,  1, 2));
  g.push_back(FEVV::CGALPoint( 1, -1, 2));
#endif  
}


} // namespace Filters
} // namespace FEVV
