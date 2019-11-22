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

#include <stdexcept> // for std::runtime_error


namespace FEVV {

// --------------- Point map ---------------

//note: get(boost::vertex_point_t,...) must return a shallow copy of the point
//      map and NOT a reference, so we must create a real PCLPointCloudPointMap
//      class that:
//       - implement shallow copy
//       - implement operator[]
//      see
//        https://github.com/CGAL/cgal/blob/ea20dfd63fcdec0b98258c9c47b0cbb88cdb356c/BGL/include/CGAL/boost/graph/properties_OpenMesh.h#L185
//      and
//        https://www.boost.org/doc/libs/1_57_0/libs/property_map/doc/property_map.html
//        https://www.boost.org/doc/libs/1_57_0/libs/property_map/doc/LvaluePropertyMap.html
class PCLPointCloudPointMap
{
public:
  typedef FEVV::PCLPoint                       value_type;
  typedef FEVV::PCLPoint&                      reference;
  typedef typename boost::graph_traits< FEVV::PCLPointCloud >::vertex_descriptor
      key_type;
  typedef boost::read_write_property_map_tag   category;

  PCLPointCloudPointMap(FEVV::PCLPointCloud &pc) : m_pc(&pc) {}

  PCLPointCloudPointMap(const FEVV::PCLPointCloudPointMap &pm) : m_pc(pm.m_pc)
  {
  }

  value_type operator[](key_type k) const
  {
    // 'points' attribute of pcl::PointCloud<>
    // is a std::vector<FEVV::PCLEnrichedPoint>
    const FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
    return FEVV::PCLPoint(enriched_point.x,
                          enriched_point.y,
                          enriched_point.z);
  }

  void set(key_type k, const FEVV::PCLPoint &point)
  {
    // 'points' attribute of pcl::PointCloud<>
    // is a std::vector<FEVV::PCLEnrichedPoint>
    FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
    enriched_point.x = point.x;
    enriched_point.y = point.y;
    enriched_point.z = point.z;
  }

private:
  FEVV::PCLPointCloud *m_pc; //TODO-elo replace by a reference
};

} // namespace FEVV


namespace pcl {

/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 * \note   Must be in 'pcl' namespace because the real type of
 *         FEVV::PCLPointCloud is pcl::PointCloud<...> and the
 *        'get(vertex_point,...)' must be in the same namespace as one
 *         of its parameters.
 */
inline
FEVV::PCLPointCloudPointMap
get(const boost::vertex_point_t, FEVV::PCLPointCloud &pc)
{
  FEVV::PCLPointCloudPointMap pm(pc);
  return pm;
}

} // namespace pcl


namespace FEVV {

//! Specialization of get(point_map, key) for PCLPointCloud
inline
FEVV::PCLPointCloudPointMap::value_type
get(const FEVV::PCLPointCloudPointMap &pm,
    FEVV::PCLPointCloudPointMap::key_type key)
{
  return pm[key];
}

//! Specialization of put(point_map, key, value) for PCLPointCloud
inline
void
put(FEVV::PCLPointCloudPointMap &pm,
    FEVV::PCLPointCloudPointMap::key_type key,
    const FEVV::PCLPointCloudPointMap::value_type &value)
{
  pm.set(key, value);
}

// --------------- Normal map ---------------

// implement a normal map with shallow copy
// normal is stored inside PCLEnrichedPoint structure
class PCLPointCloudNormalMap
{
public:
  typedef FEVV::PCLNormal                      value_type;
  typedef FEVV::PCLNormal&                     reference;
  typedef typename boost::graph_traits< FEVV::PCLPointCloud >::vertex_descriptor
      key_type;
  typedef boost::read_write_property_map_tag   category;

  PCLPointCloudNormalMap() : m_pc(nullptr) {}

  PCLPointCloudNormalMap(FEVV::PCLPointCloud &pc) : m_pc(&pc) {}

  // shallow copy
  PCLPointCloudNormalMap(const FEVV::PCLPointCloudNormalMap &nm) : m_pc(nm.m_pc)
  {
  }

  value_type operator[](key_type k) const
  {
    // 'points' attribute of pcl::PointCloud<> 
    // is a std::vector<FEVV::PCLEnrichedPoint>
    if(m_pc)
    {
      const FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
      return FEVV::PCLNormal(enriched_point.normal_x,
                             enriched_point.normal_y,
                             enriched_point.normal_z);
    }
    else
    {
      throw std::runtime_error(
          "FEVV::PCLPointCloudNormalMap: invalid use of un-initialized map");
      return FEVV::PCLNormal(0, 0, 0); // dummy return to avoid a warning
    }
  }

  void set(key_type k, const FEVV::PCLNormal &normal)
  {
    // 'points' attribute of pcl::PointCloud<>
    // is a std::vector<FEVV::PCLEnrichedPoint>
    if(m_pc)
    {
      FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
      enriched_point.normal_x = normal.x;
      enriched_point.normal_y = normal.y;
      enriched_point.normal_z = normal.z;
    }
    else
    {
      throw std::runtime_error("Invalid un-initialized PCLPointCloudNormalMap");
    }
  }

private:
  FEVV::PCLPointCloud *m_pc;
};

//! Specialization of get(normal_map, key) for PCLPointCloud
inline
FEVV::PCLPointCloudNormalMap::value_type
get(const FEVV::PCLPointCloudNormalMap &nm,
    FEVV::PCLPointCloudNormalMap::key_type key)
{
  return nm[key];
}

//! Specialization of put(normal_map, key, value) for PCLPointCloud
inline
void
put(FEVV::PCLPointCloudNormalMap &nm,
    FEVV::PCLPointCloudNormalMap::key_type key,
    const FEVV::PCLPointCloudNormalMap::value_type &value)
{
  nm.set(key, value);
}

// --------------- Color map ---------------

// implement a color map with shallow copy
// color is stored inside PCLEnrichedPoint structure
class PCLPointCloudColorMap
{
public:
  typedef FEVV::PCLColor                       value_type;
  typedef FEVV::PCLColor&                      reference;
  typedef typename boost::graph_traits< FEVV::PCLPointCloud >::vertex_descriptor
      key_type;
  typedef boost::read_write_property_map_tag   category;

  PCLPointCloudColorMap() : m_pc(nullptr) {}

  PCLPointCloudColorMap(FEVV::PCLPointCloud &pc) : m_pc(&pc) {}

  // shallow copy
  PCLPointCloudColorMap(const FEVV::PCLPointCloudColorMap &cm) : m_pc(cm.m_pc)
  {
  }

  value_type operator[](key_type k) const
  {
    // 'points' attribute of pcl::PointCloud<> 
    // is a std::vector<FEVV::PCLEnrichedPoint>
    if(m_pc)
    {
      const FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
      return FEVV::PCLColor(enriched_point.r,
                            enriched_point.g,
                            enriched_point.b);
    }
    else
    {
      throw std::runtime_error(
          "FEVV::PCLPointCloudColorMap: invalid use of un-initialized map");
      return FEVV::PCLColor(0, 0, 0); // dummy return to avoid a warning
    }
  }

  void set(key_type k, const FEVV::PCLColor &color)
  {
    // 'points' attribute of pcl::PointCloud<>
    // is a std::vector<FEVV::PCLEnrichedPoint>
    if(m_pc)
    {
      FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
      enriched_point.r = color.r;
      enriched_point.g = color.g;
      enriched_point.b = color.b;
    }
    else
    {
      throw std::runtime_error("Invalid un-initialized PCLPointCloudNormalMap");
    }
  }

private:
  FEVV::PCLPointCloud *m_pc;
};

//! Specialization of get(color_map, key) for PCLPointCloud
inline
FEVV::PCLPointCloudColorMap::value_type
get(const FEVV::PCLPointCloudColorMap &cm,
    FEVV::PCLPointCloudColorMap::key_type key)
{
  return cm[key];
}

//! Specialization of put(color_map, key, value) for PCLPointCloud
inline
void
put(FEVV::PCLPointCloudColorMap &cm,
    FEVV::PCLPointCloudColorMap::key_type key,
    const FEVV::PCLPointCloudColorMap::value_type &value)
{
  cm.set(key, value);
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
