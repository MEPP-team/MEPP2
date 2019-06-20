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

#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_off_points.h>
#if 0 //ELO-note: link error multiple definition of
      // `void CGAL::internal::PLY::property_header_type...`
#include <CGAL/IO/write_ply_points.h>
#endif
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

  std::ofstream out(filename);
  if(! out)
  {
    throw std::invalid_argument(
        "write_mesh() error: can not open output file " + filename);
  }

  out.precision(17);

  if(FEVV::FileUtils::has_extension(filename, ".xyz"))
  {
     success = CGAL::write_xyz_points(out, g);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".off"))
  {
     success = CGAL::write_off_points(out, g);
  }
#if 0 //ELO-note: link error multiple definition of
      // `void CGAL::internal::PLY::property_header_type...`
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
     success = CGAL::write_ply_points(out, g);
  }
#endif
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

