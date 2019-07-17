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

} // namespace CGAL


namespace boost {

template<>
struct property_traits< FEVV::CGALPointSet >
{
  typedef FEVV::CGALPointSetPointMap::value_type    value_type;
  typedef FEVV::CGALPointSetPointMap::reference     reference;
};

} // namespace boost
