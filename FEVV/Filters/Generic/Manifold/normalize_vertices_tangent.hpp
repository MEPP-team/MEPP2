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


namespace FEVV {
namespace Filters {

/**
 * \brief  Normalize all computed vertex tangents over the mesh.
 *
 * \param  g    the halfedge graph from which we take the mesh's informations.
 * \param  vtm  vertex tangent map in which we normalize tangents.
 * \param  gt   geometry traits object used to normalize a Vector.
 */
template< typename HalfedgeGraph,
          typename VertexTangentMap,
          typename GeometryTraits >
void
normalize_vertices_tangent(const HalfedgeGraph &g,
                           VertexTangentMap &vtm,
                           const GeometryTraits &gt)
{
  using Vector = typename GeometryTraits::Vector;

  using vertex_iterator =
      typename boost::graph_traits< HalfedgeGraph >::vertex_iterator;

  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb)
  {
    const Vector tan = gt.normalize(get(vtm, *vb));
    put(vtm, *vb, tan);
  }
}

} // namespace Filters
} // namespace FEVV
