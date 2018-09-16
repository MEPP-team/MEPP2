#pragma once

#include <boost/graph/graph_traits.hpp>

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/Geometry/tangents.hpp"


namespace FEVV {
namespace Operators {

/**
 * \brief  Calculate a face tangent, then accumulate its value on every
 *         vertices it contains to average them; final vertex tangents are
 *         normalized later.
 *
 * \param  f    the face descriptor representing the face on which we compute
 *              the tangent.
 * \param  g    the halfedge graph from which we take the mesh's informations.
 * \param  pm   point map containing vertices' positions.
 * \param  uv1  first vertex's texture coordinates.
 * \param  uv2  second vertex's texture coordinates.
 * \param  uv3  third vertex's texture coordinates.
 * \param  vtm  vertex tangent map to be filled with computed tangents.
 *
 * \ingroup  GenericManifoldFilters
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename Vector,
          typename VertexTangentMap >
void
calculate_face_tangent(
    typename boost::graph_traits< HalfedgeGraph >::face_descriptor f,
    const HalfedgeGraph &g,
    const PointMap &pm,
    const Vector &uv1,
    const Vector &uv2,
    const Vector &uv3,
    VertexTangentMap &vtm)
{
  using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
  using halfedge_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor;
  using vertex_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor;

  // RM: basic test with nth vertices polygons, *this may not work properly!*
  //   It may be needed to have only triangles for this to work (though
  //   unlikely, as long as we recover the plane basis)
  //     -> After having tested multiple meshes, some having 5+ vertices per
  //     face, it is less & less likely that this has an impact.

  /*

  *------* V2
  |     /|
  |    / |
  |   /  |
  |  /   |
  | /    |
  |/     |
  *------* V3
  V1

  Any halfedge is fine, but assuming for example we're getting halfedge (named
  H) V2 -> V1, then points are as such:
    - V1 is the target of H
    - V2 is the source of H
    - V3 is the target of the next halfedge from H

  Actually, we don't care much about which points we get, as long as
  there are 3 of them linked to each other to deduce a plane from it

  */

  halfedge_descriptor half_edge = halfedge(f, g);
  const vertex_descriptor vert1 = target(half_edge, g);          // V1
  const vertex_descriptor vert2 = source(half_edge, g);          // V2
  const vertex_descriptor vert3 = target(next(half_edge, g), g); // V3

  // RM: get positions & texcoords from vertex descriptors
  const Point pos1 = get(pm, vert1);
  const Point pos2 = get(pm, vert2);
  const Point pos3 = get(pm, vert3);

  const Vector tan =
      Operators::calculate_triangle_tangent(pos1, pos2, pos3, uv1, uv2, uv3);

  const halfedge_descriptor first_he = half_edge;
  do
  {
    const vertex_descriptor vert = source(half_edge, g);

    const Vector new_tan = get(vtm, vert) + tan;
    put(vtm, vert, new_tan);

    half_edge = next(half_edge, g);
  } while(half_edge != first_he);
}

} // namespace Operators
} // namespace FEVV
