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

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/Generic/Manifold/calculate_face_tangent.hpp"
#include "FEVV/Filters/Generic/Manifold/normalize_vertices_tangent.hpp"


namespace FEVV {
namespace Filters {

/**
 * \brief  Calculate the mesh's tangents from its halfedges.
 *
 * \param  g    the halfedge graph from which we take the mesh's informations.
 * \param  pm   vertices' positions.
 * \param  hum  halfedge UV map containing halfedges' texture coordinates.
 * \param  vtm  vertex tangent map to be filled with computed tangents.
 *
 * \sa     calculate_vertices_tangent() to compute tangents from vertices.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename HalfedgeUVMap,
          typename VertexTangentMap >
void
calculate_halfedges_tangent(const HalfedgeGraph &g,
                            const PointMap &pm,
                            const HalfedgeUVMap &hum, // provided halfedge UVs
                            VertexTangentMap &vtm)
{
  using GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >;
  using Vector = typename GeometryTraits::Vector;

  using GraphTraits = typename boost::graph_traits< HalfedgeGraph >;
  using face_iterator = typename GraphTraits::face_iterator;

  using halfedge_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor;

  const GeometryTraits gt(g);

  face_iterator fb, fe;
  for(boost::tie(fb, fe) = faces(g); fb != fe; ++fb)
  {
    halfedge_descriptor half_edge = halfedge(*fb, g);

    const Vector uv1 = get(hum, half_edge);
    const Vector uv2 = get(hum, prev(half_edge, g));
    const Vector uv3 = get(hum, next(half_edge, g));

    Operators::calculate_face_tangent(*fb, g, pm, uv1, uv2, uv3, vtm);
  }

  Filters::normalize_vertices_tangent(g, vtm, gt);
}

} // namespace Filters
} // namespace FEVV
