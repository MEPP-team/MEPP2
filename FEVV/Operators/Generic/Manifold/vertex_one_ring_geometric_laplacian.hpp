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
#include <vector>

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/Generic/Manifold/k_ring.hpp"


namespace FEVV {
namespace Operators {

/**
 * \brief  Geometric Laplacian of the vertex' one ring of vertex positions.
 *
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<> specialization.
 * \tparam  PointMap A modifiable point map to manage vertex positions.
 * \tparam  GeometryTraits The geometric kernel when available. This is
 *          defaulted to FEVV::Geometry_traits<FaceGraph>.
 * \param[in] v The vertex whose one-ring geometric laplacian is computed.
 * \param[in] g The FaceGraph instance from which the v vertex will be taken.
 * \param[in] pm The point map which associate a vertex descriptor with a vertex
 *            position.
 * \param[in] smoothing_factor A (positive/move towards or negative/move away from)
 *            factor used to compute the position of the returned geometric laplacian
 *            V_pos + smoothing_factor x(GL_computed - V_pos).
 * \param[in] gt The geometry trait object.
 * \return The geometric laplacian of the vertex v.
 * \pre     g must be a triangular mesh.
 */
template< typename FaceGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
typename boost::property_traits< PointMap >::value_type
vertex_one_ring_geometric_laplacian(
    typename boost::graph_traits< FaceGraph >::vertex_descriptor v,
    const FaceGraph &g,
    const PointMap &pm,
    const typename GeometryTraits::Scalar smoothing_factor,
    const GeometryTraits &gt)
{
  typedef typename boost::property_traits< PointMap >::value_type
      Point; // this Point type cannot be const because the center Point must be
             // modifiable
  typedef typename boost::property_traits< PointMap >::reference Reference;
  typedef typename GeometryTraits::Scalar Scalar;

  std::vector< typename boost::graph_traits< FaceGraph >::vertex_descriptor >
      qv;
  Operators::extract_1_ring_not_including_v< FaceGraph, GeometryTraits >(v, g, qv);

  Point geom_laplacian(0, 0, 0);
  Scalar sum = 0;
  Point center = get(pm, v); // works for all mesh structures
  unsigned int cpt = 0;
  typename std::vector< typename boost::graph_traits<
      const FaceGraph >::vertex_descriptor >::const_iterator it(qv.begin()),
      ite(qv.end());
  for(; it != ite; ++it, ++cpt)
  {
    Point tmp = get(pm, *it); // works for all mesh structures
    Scalar len = gt.length(center, tmp);
    if(::fabs(len) > std::numeric_limits<typename GeometryTraits::Scalar>::epsilon())
    {
      sum += 1./len ;
      geom_laplacian = Point(gt.get_x(tmp)/len + gt.get_x(geom_laplacian),
                             gt.get_y(tmp)/len + gt.get_y(geom_laplacian),
                             gt.get_z(tmp)/len + gt.get_z(geom_laplacian));
    }
  }
  if(cpt > 0)
  {
    return gt.add_pv(center,
                    gt.scalar_mult(gt.sub(Point(gt.get_x(geom_laplacian) / sum,
                                                gt.get_y(geom_laplacian) / sum,
                                                gt.get_z(geom_laplacian) / sum),
                                          center),
                                   smoothing_factor));
  }
  else
    return get(pm, v);
}

} // namespace Operators
} // namespace FEVV
