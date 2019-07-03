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


namespace FEVV {

//note: get(boost::vertex_point_t,...) must return a shallow copy of the point
//      map and NOT a reference, so we must create a real PCLPointCloudPointMap
//      class that:
//       - implement shallow copy
//       - implement operator[]
//      see
//             https://github.com/CGAL/cgal/blob/ea20dfd63fcdec0b98258c9c47b0cbb88cdb356c/BGL/include/CGAL/boost/graph/properties_OpenMesh.h#L185
//      and
//        https://www.boost.org/doc/libs/1_57_0/libs/property_map/doc/property_map.html
class PCLPointCloudPointMap
{
public:
  typedef boost::read_write_property_map_tag   category;
  typedef FEVV::PCLPoint    value_type;
  typedef FEVV::PCLPoint&   reference;
  typedef typename boost::graph_traits< FEVV::PCLPointCloud >::vertex_descriptor
      key_type;

  PCLPointCloudPointMap() : m_pc(NULL) {}

  PCLPointCloudPointMap(/*const*/ FEVV::PCLPointCloud &pc) : m_pc(&pc) {}

  PCLPointCloudPointMap(const FEVV::PCLPointCloudPointMap &pm) : m_pc(pm.m_pc)
  {
  }

  value_type operator[](key_type k) const
  {
    // 'points' attribute of pcl::PointCloud<> is a std::vector<pcl::PointXYZ>
    return m_pc->points[k];
  }

  void set(key_type k, const value_type& v)
  {
    // 'points' attribute of pcl::PointCloud<> is a std::vector<pcl::PointXYZ>
    m_pc->points[k] = v;
  }

private:
  /*const*/ FEVV::PCLPointCloud *m_pc;
};


/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 */
inline
PCLPointCloudPointMap
get(const boost::vertex_point_t, FEVV::PCLPointCloud &pc)
{
  FEVV::PCLPointCloudPointMap pm(pc);
  return pm;
}

//! Specialization of get(point_map, key) for PCLPointCloud
inline
PCLPointCloudPointMap::value_type
get(const PCLPointCloudPointMap &pm,
    PCLPointCloudPointMap::key_type key)
{
  return pm[key];
}

//! Specialization of put(point_map, key, value) for PCLPointCloud
inline
void
put(PCLPointCloudPointMap &pm,
    PCLPointCloudPointMap::key_type key,
    const PCLPointCloudPointMap::value_type &value)
{
  pm.set(key, value);
}

} // namespace FEVV


namespace boost {

template<>
struct property_traits< FEVV::PCLPointCloudPointMap >
{
  typedef FEVV::PCLPointCloudPointMap::value_type    value_type;
  typedef FEVV::PCLPointCloudPointMap::reference     reference;
};

} // namespace boost
