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

#if 0
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_ply_points.h>
#else
#include "read_xyz_points_patched.h"
#include "read_off_points_patched.h"
#include "read_ply_points_patched.h"
#endif

#if 0 //ELO-note: doesn't compile, extra dependency needed
#include <CGAL/IO/read_las_points.h>
#endif
#include <stdexcept> // for std::invalid_argument


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
 *
 * \ingroup  GenericManifoldFilters GenericNonManifoldFilters
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

  // create normal and color property maps
  using VertexNormalMap =
      typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                  FEVV::CGALPointSet >::pmap_type;
  using VertexColorMap =
      typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                  FEVV::CGALPointSet >::pmap_type;
  VertexNormalMap v_nm;
  VertexColorMap v_cm;

  // retrieve point map
  auto pm = get(boost::vertex_point, g);

  // load point cloud
  if(FEVV::FileUtils::has_extension(filename, ".xyz"))
  {
    bool normals_found;

    // load geometry + normal
    success = CGAL::read_xyz_points(
        in,
        std::back_inserter(g),
        normals_found,
        CGAL::parameters::point_map(pm).normal_map(v_nm));

    if(normals_found)
      put_property_map(FEVV::vertex_normal, g, pmaps, v_nm);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".off"))
  {
    bool normals_found;

    // load geometry + normal
    success = CGAL::read_off_points(
        in,
        std::back_inserter(g),
        normals_found,
        CGAL::parameters::point_map(pm).normal_map(v_nm));

    if(normals_found)
      put_property_map(FEVV::vertex_normal, g, pmaps, v_nm);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
    bool normals_found, color_found;

    success = CGAL::read_ply_points_with_properties(
        in,
        g,
        normals_found,
        color_found,
        CGAL::make_ply_point_reader(pm),
        CGAL::make_ply_normal_reader(v_nm),
        std::make_tuple(v_cm,
                        typename CGAL::Kernel_traits<
                            typename VertexColorMap::value_type >::Kernel::
                            Construct_vector_3(),
                        CGAL::PLY_property< double >("red"),
                        CGAL::PLY_property< double >("green"),
                        CGAL::PLY_property< double >("blue")));

    if(normals_found)
      put_property_map(FEVV::vertex_normal, g, pmaps, v_nm);

    if(color_found)
      put_property_map(FEVV::vertex_color, g, pmaps, v_cm);
  }
#if 0 //ELO-note: doesn't compile, extra dependency needed
  else if(FEVV::FileUtils::has_extension(filename, ".las") ||
          FEVV::FileUtils::has_extension(filename, ".laz"))
  {
     success = CGAL::read_las_points(
         in,
         std::back_inserter(g));
  }
#endif

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
