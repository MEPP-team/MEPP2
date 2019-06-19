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

#include "FEVV/Filters/Generic/generic_writer.hpp" // specialization of
#include "FEVV/Tools/IO/FileUtilities.hpp" // for FileUtils::has_extension()

#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>

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
write_mesh< FEVV::PCLPointCloud, FEVV::Geometry_traits< FEVV::PCLPointCloud > >(
    const std::string &filename,
    FEVV::PCLPointCloud &g,
    PMapsContainer &pmaps)
{
  bool success = false;

  if(FEVV::FileUtils::has_extension(filename, ".pcd"))
  {
    success = (pcl::io::savePCDFile(filename, g) >= 0);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
    success = (pcl::io::savePLYFile(filename, g) >= 0);
  }

  if(! success)
  {
    throw std::invalid_argument(
        "write_mesh() error: failed to write mesh to file " + filename);
  }
}


} // namespace Filters
} // namespace FEVV

