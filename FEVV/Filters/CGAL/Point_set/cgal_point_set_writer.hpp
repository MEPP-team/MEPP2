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

#include "FEVV/Filters/Generic/generic_writer.hpp" // specialization of
#include "FEVV/Tools/IO/FileUtilities.hpp" // for FileUtils::has_extension()

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/write_ply_points.h>
#if 0 //ELO-note: doesn't compile, extra dependency needed
#include <CGAL/IO/write_las_points.h>
#endif

#include <stdexcept> // for std::invalid_argument


namespace FEVV {
namespace Filters {


/**
 * \brief  Write mesh to file.
 *
 * \param  filename    name of the output mesh file
 * \param  g           mesh to write to file
 * \param  pmaps       property maps bag of the mesh
 */
template<>
inline
void
write_mesh< FEVV::CGALPointSet, FEVV::Geometry_traits< FEVV::CGALPointSet > >(
    const std::string &filename,
    FEVV::CGALPointSet &g,
    PMapsContainer &pmaps)
{
  bool success = false;

  // init output file
  std::ofstream out(filename);
  if(! out)
  {
    throw std::invalid_argument(
        "write_mesh() error: can not open output file " + filename);
  }
  //out.precision(16);

  // save point cloud
  if(FEVV::FileUtils::has_extension(filename, ".xyz"))
  {
    // colors are not supported with this file format

    // write geometry and normals if present
    success = CGAL::write_xyz_point_set(out, g);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".off"))
  {
    // colors are not supported with this file format

    // write geometry and normals if present
    success = CGAL::write_off_point_set(out, g);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
    // retrieve property maps
    using VertexNormalMap =
        typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                    FEVV::CGALPointSet >::pmap_type;
    using VertexColorMap =
        typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                    FEVV::CGALPointSet >::pmap_type;
    
    VertexNormalMap v_nm;
    bool has_normal = has_map(pmaps, vertex_normal);
    if(has_normal)
      v_nm = get_property_map(FEVV::vertex_normal, g, pmaps);

    VertexColorMap v_cm;
    bool has_color = has_map(pmaps, vertex_color);
    if(has_color)
      v_cm = get_property_map(FEVV::vertex_color, g, pmaps);
    
    auto pm = get(boost::vertex_point, g);

    // set file mode
    CGAL::set_ascii_mode(out);

    // write to file
    if(has_normal && has_color)
    {
      // geometry + normal + color
      success = CGAL::write_ply_points_with_properties(
          out,
          g,
          CGAL::make_ply_point_writer(pm),
          CGAL::make_ply_normal_writer(v_nm),
          std::make_tuple(v_cm,
                          CGAL::PLY_property< unsigned char >("red"),
                          CGAL::PLY_property< unsigned char >("green"),
                          CGAL::PLY_property< unsigned char >("blue")));
    }
    else if(has_normal)
    {
      // geometry + normal
      success = CGAL::write_ply_points_with_properties(
          out,
          g,
          CGAL::make_ply_point_writer(pm),
          CGAL::make_ply_normal_writer(v_nm));
    }
    else if(has_color)
    {
      // geometry + color
      success = CGAL::write_ply_points_with_properties(
          out,
          g,
          CGAL::make_ply_point_writer(pm),
          std::make_tuple(v_cm,
                          CGAL::PLY_property< unsigned char >("red"),
                          CGAL::PLY_property< unsigned char >("green"),
                          CGAL::PLY_property< unsigned char >("blue")));
    }
    else
    {
      // geometry only
      success = CGAL::write_ply_points_with_properties(
          out,
          g,
          CGAL::make_ply_point_writer(pm));
    }
  }
#if 0 //ELO-note: doesn't compile, extra dependency needed
  else if(FEVV::FileUtils::has_extension(filename, ".las") ||
          FEVV::FileUtils::has_extension(filename, ".laz"))
  {
    success = CGAL::write_las_points(out, points);
  }
#endif

  if(! success)
  {
    throw std::invalid_argument(
        "write_mesh() error: failed to write mesh to file " + filename);
  }
}


} // namespace Filters
} // namespace FEVV

