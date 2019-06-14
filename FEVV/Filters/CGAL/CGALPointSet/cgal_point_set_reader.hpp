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

//TODO-elo-rm  #include <boost/graph/graph_traits.hpp>
//TODO-elo-rm  #include <boost/graph/properties.hpp>
//TODO-elo-rm  #include "FEVV/Wrappings/Graph_traits_extension.h"
//TODO-elo-rm  #include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"

#include "FEVV/Filters/Generic/generic_reader.hpp" // specialization of
#include "FEVV/Tools/IO/FileUtilities.hpp" // for FileUtils::has_extension()

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_ply_points.h>
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
  //TODO-elo-clean
  std::cout << "here " << __func__ << std::endl; //TODO-elo-rm
#if 1
  bool success = false;

  std::ifstream in(filename);
  if(! in)
  {
    throw std::invalid_argument("read_mesh() error: can not open file " +
                                filename);
  }

  if(FEVV::FileUtils::has_extension(filename, ".xyz"))
  {
     success = CGAL::read_xyz_points(
         in,
         std::back_inserter(g));
  }
  else if(FEVV::FileUtils::has_extension(filename, ".off"))
  {
     success = CGAL::read_off_points(
         in,
         std::back_inserter(g));
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
     success = CGAL::read_ply_points(
         in,
         std::back_inserter(g));
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
  //TODO-elo-rm DBG
  g.push_back(FEVV::CGALPoint(-1, -1, 2));
  g.push_back(FEVV::CGALPoint(-1,  1, 2));
  g.push_back(FEVV::CGALPoint( 1,  1, 2));
  g.push_back(FEVV::CGALPoint( 1, -1, 2));
#endif  
}


} // namespace Filters
} // namespace FEVV
