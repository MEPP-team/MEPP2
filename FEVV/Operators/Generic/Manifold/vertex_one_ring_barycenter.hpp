#pragma once

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <vector>

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/Generic/Manifold/k_ring.hpp"


namespace FEVV {
namespace Operators {

/**
 * \brief  Barycenter of the vertex' one ring of vertex positions.
 *
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<> specialization.
 * \tparam  PointMap A modifiable point map to manage vertex positions.
 * \tparam  GeometryTraits The geometric kernel when available. This is
 *          defaulted to FEVV::Geometry_traits<FaceGraph>.
 * \param[in] v The vertex whose one-ring barycenter is computed.
 * \param[in] g The FaceGraph instance from which the v vertex will be taken.
 * \param[in] pm The point map which associate a vertex descriptor with a vertex
 *            position.
 * \param[in] SMOOTHING_FACTOR A (usually positive) factor used to compute the
 *            position of the returned centroid following
 *            V_pos + SMOOTHING_FACTORx(C_computed - V_pos).
 * \param[in] gt The geometry trait object.
 * \return The barycenter centroid of the vertex v.
 * \pre     g must be a triangular mesh.
 * \ingroup GenericManifoldFilters
 */
template< typename FaceGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
typename boost::property_traits< PointMap >::value_type
vertex_one_ring_barycenter(
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

  std::vector< typename boost::graph_traits< FaceGraph >::vertex_descriptor >
      qv;
  Operators::extract_1_ring_not_including_v< FaceGraph, GeometryTraits >(v, g, qv);

  Point center(0, 0, 0);
  unsigned int cpt = 0;
  typename std::vector< typename boost::graph_traits<
      const FaceGraph >::vertex_descriptor >::const_iterator it(qv.begin()),
      ite(qv.end());
  for(; it != ite; ++it, ++cpt)
  {
    Point tmp = get(pm, *it); // works for all mesh structures
    center = Point(gt.get_x(tmp) + gt.get_x(center),
                   gt.get_y(tmp) + gt.get_y(center),
                   gt.get_z(tmp) + gt.get_z(center));
  }
  if(cpt > 0)
  {
    return gt.add_p(get(pm, v),
                    gt.scalar_mult(smoothing_factor,
                                   gt.sub(Point(gt.get_x(center) / cpt,
                                                gt.get_y(center) / cpt,
                                                gt.get_z(center) / cpt),
                                          get(pm, v))));
  }
  else
    return get(pm, v);
}

} // namespace Operators
} // namespace FEVV
