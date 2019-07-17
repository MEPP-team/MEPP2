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

using CGALPointSetPointMap = CGALPointSet::Point_map;

} // namespace FEVV


namespace CGAL {

/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 */
inline
FEVV::CGALPointSetPointMap
get(const boost::vertex_point_t, FEVV::CGALPointSet &ps)
{
  return ps.point_map();
}

#if 0//TODO-elo-rm if already ok
//! Specialization of get(point_map, key) for CGALPointSet
inline
FEVV::CGALPointSetPointMap::value_type&
get(const FEVV::CGALPointSetPointMap &pm, FEVV::CGALPointSetPointMap::key_type key)
{
  return pm[key];
}

//! Specialization of put(point_map, key, value) for CGALPointSet
inline
void
put(FEVV::CGALPointSetPointMap &pm,
    FEVV::CGALPointSetPointMap::key_type key,
    const FEVV::CGALPointSetPointMap::value_type &value)
{
  pm[key] =value;
}
#endif

} // namespace CGAL


namespace boost {

template<>
struct property_traits< FEVV::CGALPointSet >
{
  typedef FEVV::CGALPointSetPointMap::value_type    value_type;
  typedef FEVV::CGALPointSetPointMap::reference     reference;
};

} // namespace boost

#if 0 //TODO-elo-probably not needed because we will use Point_set_3 internal prop maps
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
#endif

