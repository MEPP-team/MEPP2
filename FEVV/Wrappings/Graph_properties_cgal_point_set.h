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
#include <CGAL/property_map.h> // for CGAL::Identity_property_map

namespace FEVV {

/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 */
inline
FEVV::CGALPointSet &
get(const boost::vertex_point_t, FEVV::CGALPointSet &ps)
{
  // the point map of the Point Set is the Point Set itself
  return ps;
}
#if 0 //TODO-elo-rm
inline
CGAL::Identity_property_map< FEVV::CGALPoint >
get(const boost::vertex_point_t, FEVV::CGALPointSet &ps)
{
  return CGAL::Identity_property_map< FEVV::CGALPoint >();;
}
#endif

//! Specialization of get(point_map, key) for CGALPointSet
inline
const FEVV::CGALPointSet::value_type &
get(const FEVV::CGALPointSet &pm, FEVV::CGALPointSet::size_type key)
{
  return pm[key];
}

//! Specialization of put(point_map, key, value) for CGALPointSet
inline void
put(FEVV::CGALPointSet &pm,
    const FEVV::CGALPointSet::size_type &key,
    const FEVV::CGALPointSet::value_type &value)
{
  pm[key] = value;
}

} // namespace FEVV


namespace boost {

template<>
struct property_traits< FEVV::CGALPointSet >
{
  typedef FEVV::CGALPointSet::value_type    value_type;
  typedef FEVV::CGALPointSet::reference      reference;
};

} // namespace boost
