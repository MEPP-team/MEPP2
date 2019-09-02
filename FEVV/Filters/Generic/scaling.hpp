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
 * \brief  Scale a mesh
 *
 * \param  g                 input mesh
 * \param  pm                output point map
 * \param  scaleX            x scale factor
 * \param  scaleY            y scale factor
 * \param  scaleZ            z scale factor
 * \param  gt                the geometry traits to use
 *
 * \sa      the simplified variant that use the default geometry traits
 *          of the mesh.
 */
template< typename Graph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< Graph > >
void
calculate_scaling(Graph &g,
                  PointMap &pm,
                  typename GeometryTraits::Scalar scale_x,
                  typename GeometryTraits::Scalar scale_y,
                  typename GeometryTraits::Scalar scale_z,
                  const GeometryTraits &gt)
{
  typedef boost::graph_traits< Graph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename boost::property_traits< PointMap >::value_type Point;
  typedef typename boost::property_traits< PointMap >::reference Reference;

  // std::cout << "Asking to scale mesh to " << scaleX << ";" << scaleY << ";"
  // << scaleZ << std::endl;

  auto iterator_pair = vertices(g); // vertices() returns a vertex_iterator pair
  vertex_iterator vi = iterator_pair.first;
  vertex_iterator vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    Point p = get(pm, *vi);
    put(pm,
        *vi,
        Point(gt.get_x(p) * scale_x,
              gt.get_y(p) * scale_y,
              gt.get_z(p) * scale_z));
  }
}

/**
 * \brief  Scale a mesh
 *
 * \param  g                 input mesh
 * \param  pm                output point map
 * \param  scaleX            x scale factor
 * \param  scaleY            y scale factor
 * \param  scaleZ            z scale factor
 *
 * \sa      the variant that use the geometry traits provided by the user.
 */
template< typename Graph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< Graph > >
void
calculate_scaling(Graph &g,
                  PointMap &pm,
                  typename GeometryTraits::Scalar scale_x,
                  typename GeometryTraits::Scalar scale_y,
                  typename GeometryTraits::Scalar scale_z)

{
  GeometryTraits gt(g);
  calculate_scaling< Graph, PointMap, GeometryTraits >(
      g, pm, scale_x, scale_y, scale_z, gt);
}

} // namespace Filters
} // namespace FEVV
