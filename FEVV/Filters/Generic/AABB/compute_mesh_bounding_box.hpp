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
 * \brief Compute the mesh Axis-Aligned Bounding-Box (AABB).
 *
 * \tparam PropertyGraph a Mesh type that provides a Model of the
 *         PropertyGraph Concept through a boost::graph_traits<> specialization.
 * \tparam PointMap A modifiable point map to manage vertex positions.
 * \tparam GeometryTraits The geometric kernel when available. This is defaulted
 *         to FEVV::Geometry_traits<PropertyGraph>.
 * \param  g The PropertyGraph instance from which the vertices will be taken.
 * \param  pm The point map which associate a vertex descriptor with a vertex
 *         position.
 * \param  smoothed_pm  The point map used to replace the vertex positions in
 *         pm.
 * \param  gt The geometry trait object.
 */
template< typename PropertyGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< PropertyGraph > >
void
compute_mesh_bounding_box(
    const PropertyGraph &g,
    const PointMap &pm,
    typename boost::property_traits< PointMap >::value_type &min_aabb,
    typename boost::property_traits< PointMap >::value_type &max_aabb,
    const GeometryTraits &gt)
{
  typedef typename boost::property_traits< PointMap >::value_type Point;
  typedef typename boost::property_traits< PointMap >::reference Reference;

  typedef boost::graph_traits< PropertyGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  auto iterator_pair = vertices(g); // vertices() returns a vertex_iterator pair
  vertex_iterator vi = iterator_pair.first;
  vertex_iterator vi_end = iterator_pair.second;
  if(vi != vi_end)
  {
    Point p = get(pm, *vi);
    min_aabb = Point(gt.get_x(p), gt.get_y(p), gt.get_z(p));
    max_aabb = Point(gt.get_x(p), gt.get_y(p), gt.get_z(p));
    ++vi;
    for(; vi != vi_end; ++vi)
    {
      Point p = get(pm, *vi);
      min_aabb = Point(std::min(gt.get_x(p), gt.get_x(min_aabb)),
                       std::min(gt.get_y(p), gt.get_y(min_aabb)),
                       std::min(gt.get_z(p), gt.get_z(min_aabb)));

      max_aabb = Point(std::max(gt.get_x(p), gt.get_x(max_aabb)),
                       std::max(gt.get_y(p), gt.get_y(max_aabb)),
                       std::max(gt.get_z(p), gt.get_z(max_aabb)));
    }
  }
}

} // namespace Filters
} // namespace FEVV
