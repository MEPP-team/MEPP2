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

#include "FEVV/Wrappings/Graph_traits_pcl_point_cloud.h"
#include <CGAL/boost/graph/properties.h> // for boost::vertex_point_t
#include <CGAL/property_map.h> // for CGAL::Identity_property_map

namespace FEVV {

// see
//   https://github.com/PointCloudLibrary/pcl/blob/39732f5a7c8455ed51fd0f6278d8b25322a68dd9/common/include/pcl/point_cloud.h#L410
// and  
//   https://github.com/PointCloudLibrary/pcl/blob/39732f5a7c8455ed51fd0f6278d8b25322a68dd9/common/include/pcl/point_cloud.h#L426
typedef  FEVV::PCLPointCloud::VectorType  PCLPointCloudPointMap;


/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 */
inline
PCLPointCloudPointMap &
get(const boost::vertex_point_t, FEVV::PCLPointCloud &pc)
{
  // the point map of the PCL Point Cloud is its 'points' attribute,
  // which is a std::vector<pcl::PointXYZ> ;
  // vertex_descriptor is the index of the vertex in the std::vector,
  // so pc.points[vertex_descriptor] is the point coordinates ;
  return pc.points;
}

//! Specialization of get(point_map, key) for PCLPointCloud
inline
const FEVV::PCLPoint &
get(const PCLPointCloudPointMap &pm,
    FEVV::PCLPointCloud::VectorType::size_type key)
{
  return pm[key];
}

//! Specialization of put(point_map, key, value) for PCLPointCloud
inline
void
put(PCLPointCloudPointMap &pm,
    const FEVV::PCLPointCloud::VectorType::size_type &key,
    const FEVV::PCLPoint &value)
{
  pm[key] = value;
}

} // namespace FEVV


namespace boost {

template<>
struct property_traits< FEVV::PCLPointCloudPointMap >
{
  // see https://github.com/PointCloudLibrary/pcl/blob/39732f5a7c8455ed51fd0f6278d8b25322a68dd9/common/include/pcl/point_cloud.h#L433
  typedef FEVV::PCLPointCloud::VectorType::value_type    value_type;
  // see https://github.com/PointCloudLibrary/pcl/blob/39732f5a7c8455ed51fd0f6278d8b25322a68dd9/common/include/pcl/point_cloud.h#L434
  typedef FEVV::PCLPointCloud::VectorType::reference     reference;
};

} // namespace boost
