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

#include <Eigen/Dense>

namespace FEVV {
namespace Filters {


/**
 * \brief  Calculate the homogeneous transformation of the
 *         vertices coordinates using the provided 4x4 homogeneous
 *         transformation matrix.
 *
 * \param[in]  g     The considered mesh
 * \param[in]  pm    The Point Map holding the vertices geometry
 * \param[in]  mat   The 4x4 homogeneous transformation matrix
 * \param[in]  gt    The Geometry Traits used to perfom the geometrical
 *                   calculations
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits >
void
homogeneous_transform(const HalfedgeGraph &g,
                      PointMap &pm,
                      const Eigen::Matrix4d &mat,
                      const GeometryTraits &gt)
{
  typedef boost::graph_traits< HalfedgeGraph >               GraphTraits;
  typedef typename GraphTraits::vertex_iterator          vertex_iterator;
  typedef typename boost::property_traits< PointMap >::value_type  Point;

  auto iterator_pair = vertices(g);
  vertex_iterator vi = iterator_pair.first;
  vertex_iterator vi_end = iterator_pair.second;

  // loop over vertices
  for(; vi != vi_end; ++vi)
  {
    // retrieve vertex coordinates
    Point p = get(pm, *vi);
    double x = gt.get_x(p);
    double y = gt.get_y(p);
    double z = gt.get_z(p);
    Eigen::Vector4d v(x, y, z, 1);

    // compute transformed coordinates
    //   | s.xt |   |            tx | | x |
    //   | s.yt | = |   R3x3     ty |.| y |
    //   | s.zt |   |            tz | | z |
    //   |  s   |   | 0   0   0   1 | | 1 |
    // with s = 1
    Eigen::Vector4d vt = mat * v;
    assert(vt(3) == 1);

    // update point property map
    Point pt(vt(0), vt(1), vt(2));
    put(pm, *vi, pt);
  }
}


/**
 * \brief  Calculate the homogeneous transformation of the
 *         vertices coordinates using the provided 4x4 homogeneous
 *         transformation matrix.
 *
 * \param[in]  g     The considered mesh
 * \param[in]  pm    The Point Map holding the vertices geometry
 * \param[in]  mat   The 4x4 homogeneous transformation matrix
 * \param[in]  gt    The Geometry Traits used to perfom the geometrical
 *                   calculations
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
homogeneous_transform(const HalfedgeGraph &g,
                      PointMap &pm,
                      const Eigen::Matrix4d &mat)
{
  GeometryTraits gt(g);
  return homogeneous_transform< HalfedgeGraph,
                                PointMap,
                                GeometryTraits >(g, pm, mat, gt);
}


} // namespace Filters
} // namespace FEVV

