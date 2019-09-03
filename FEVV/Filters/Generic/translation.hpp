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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include "FEVV/Wrappings/Geometry_traits.h"

namespace FEVV {
namespace Filters {

/**
 * \brief  Translate a mesh
 *
 * \param  g                 input mesh
 * \param  pm                output point map
 * \param  offsetX           x translation offset
 * \param  offsetY           y translation offset
 * \param  offsetZ           z translation offset
 * \param  gt                the geometry traits to use
 *
 * \sa      the simplified variant that use the default geometry traits
 *          of the mesh.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
translate(const HalfedgeGraph &g,
          PointMap &pm,
          typename GeometryTraits::Scalar offset_x,
          typename GeometryTraits::Scalar offset_y,
          typename GeometryTraits::Scalar offset_z,
          const GeometryTraits &gt)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename boost::property_traits< PointMap >::value_type Point;
  typedef typename boost::property_traits< PointMap >::reference Reference;

  vertex_iterator vi;
  for(vi = vertices(g).first; vi != vertices(g).second; ++vi)
  {
    Point p = pm[*vi];
    put(pm,
        *vi,
        Point(gt.get_x(p) + offset_x,
              gt.get_y(p) + offset_y,
              gt.get_z(p) + offset_z));
  }
}

/**
 * \brief  Translate a mesh
 *
 * \param  g                 input mesh
 * \param  pm                output point map
 * \param  offsetX           x translation offset
 * \param  offsetY           y translation offset
 * \param  offsetZ           z translation offset
 *
 * \sa      the variant that use the geometry traits provided by the user.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
translate(const HalfedgeGraph &g,
          PointMap &pm,
          typename GeometryTraits::Scalar offset_x,
          typename GeometryTraits::Scalar offset_y,
          typename GeometryTraits::Scalar offset_z)

{
  GeometryTraits gt(g);
  translate< HalfedgeGraph, PointMap, GeometryTraits >(
      g, pm, offset_x, offset_y, offset_z, gt);
}

} // namespace Filters
} // namespace FEVV
