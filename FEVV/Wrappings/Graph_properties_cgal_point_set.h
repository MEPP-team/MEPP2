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

#include "FEVV/Wrappings/Graph_traits_cgal_point_set.h"
#include <CGAL/boost/graph/properties.h> // for boost::vertex_point_t


namespace FEVV {

//note: get(boost::vertex_point_t,...) must return a shallow copy of the point
//      map and NOT a reference, so we must create a real CGALPointSetPointMap
//      class that:
//       - implement shallow copy
//       - implement operator[]
//      see
//             https://github.com/CGAL/cgal/blob/ea20dfd63fcdec0b98258c9c47b0cbb88cdb356c/BGL/include/CGAL/boost/graph/properties_OpenMesh.h#L185
//      and
//        https://www.boost.org/doc/libs/1_57_0/libs/property_map/doc/property_map.html
class CGALPointSetPointMap
{
public:
  typedef boost::read_write_property_map_tag   category;
  typedef FEVV::CGALPoint    value_type;
  typedef FEVV::CGALPoint&   reference;
  typedef typename boost::graph_traits< FEVV::CGALPointSet >::vertex_descriptor
      key_type;

  CGALPointSetPointMap() : m_ps(NULL) {}

  CGALPointSetPointMap(/*const*/ FEVV::CGALPointSet &ps) : m_ps(&ps) {}

  CGALPointSetPointMap(const CGALPointSetPointMap &pm) : m_ps(pm.m_ps)
  {
  }

  value_type operator[](key_type k) const
  {
    return (*m_ps)[k];
  }

  void set(key_type k, const value_type& v)
  {
    (*m_ps)[k] = v;
  }

private:
  /*const*/ FEVV::CGALPointSet *m_ps;
};

} // namespace FEVV


namespace std {

/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 */
inline
FEVV::CGALPointSetPointMap
get(const boost::vertex_point_t, FEVV::CGALPointSet &ps)
{
  FEVV::CGALPointSetPointMap pm(ps);
  return pm;
}

} // namespace std


namespace FEVV {

//! Specialization of get(point_map, key) for CGALPointSet
inline
CGALPointSetPointMap::value_type
get(const CGALPointSetPointMap &pm, CGALPointSetPointMap::key_type key)
{
  return pm[key];
}

//! Specialization of put(point_map, key, value) for CGALPointSet
inline
void
put(CGALPointSetPointMap &pm,
    CGALPointSetPointMap::key_type key,
    const CGALPointSetPointMap::value_type &value)
{
  pm.set(key, value);
}

} // namespace FEVV


namespace boost {

template<>
struct property_traits< FEVV::CGALPointSet >
{
  typedef FEVV::CGALPointSetPointMap::value_type    value_type;
  typedef FEVV::CGALPointSetPointMap::reference     reference;
};

} // namespace boost


namespace std {

/**
 * \brief  Returns the vertex index property map of the mesh.
 * \note   This function is used for vector-property-maps,
 *         see properties_cgal_point_set.h ;
 *         it returns a property map that maps an index
 *         to a vertex descriptor ;
 *         see for example
 *          https://github.com/CGAL/cgal/blob/ea20dfd63fcdec0b98258c9c47b0cbb88cdb356c/Surface_mesh/include/CGAL/boost/graph/properties_Surface_mesh.h#L140
 *         for FEVV::CGALPointSet, vertex descriptor are already indices,
 *         so an identity_property_map should do the trick
 */
inline
boost::identity_property_map
get(const boost::vertex_index_t &, const FEVV::CGALPointSet &)
{
  return boost::identity_property_map();
}

} // namespace std

