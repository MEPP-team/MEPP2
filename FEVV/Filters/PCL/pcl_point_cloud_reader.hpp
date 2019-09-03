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

#include "FEVV/DataStructures/DataStructures_pcl_point_cloud.h"
#include "FEVV/Wrappings/Wrappings_pcl_point_cloud.h"

#include "FEVV/Wrappings/properties.h"

#include "FEVV/Filters/Generic/generic_reader.hpp" // specialization of
#include "FEVV/Tools/IO/FileUtilities.hpp" // for FileUtils::has_extension()

#include <pcl/io/ascii_io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>

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
 */
template<>
inline
void
read_mesh< FEVV::PCLPointCloud, FEVV::Geometry_traits< FEVV::PCLPointCloud > >(
    const std::string &filename,
    FEVV::PCLPointCloud &g,
    PMapsContainer &pmaps,
    bool only_pts)
{
#if 1
  bool success = false;

  if(FEVV::FileUtils::has_extension(filename, ".xyz"))
  {
    pcl::ASCIIReader reader;
    reader.setInputFields< FEVV::PCLPoint >();
    reader.setExtension(".xyz");

    success = (reader.read(filename, g) >= 0);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".pcd"))
  {
    success = (pcl::io::loadPCDFile(filename, g) >= 0);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
    success = (pcl::io::loadPLYFile(filename, g) >= 0);
  }

  if(! success)
  {
    throw std::invalid_argument("read_mesh() error: can not read file " +
                                filename);
  }
#else
  // DBG
  g.points.push_back(FEVV::PCLPoint(-1, -1, 2));
  g.points.push_back(FEVV::PCLPoint(-1,  1, 2));
  g.points.push_back(FEVV::PCLPoint( 1,  1, 2));
  g.points.push_back(FEVV::PCLPoint( 1, -1, 2));
#endif  
}


} // namespace Filters
} // namespace FEVV
