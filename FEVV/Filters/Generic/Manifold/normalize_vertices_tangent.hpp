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
 *
 * \ingroup  GenericManifoldFilters
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
