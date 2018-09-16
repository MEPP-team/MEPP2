#ifndef CALCULATE_FACE_NORMALS_H
#define CALCULATE_FACE_NORMALS_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/Generic/calculate_face_normal.hpp"

namespace FEVV {
namespace Filters {
/**
 * \brief  Calculate "some" normal of all the faces of the considered
 *         mesh and populate the argument provided Face Normal Map
 *         with the resulting normals.
 *
 * \param[in]  g    The considered mesh
 * \param[in]  pm   The Point Map holding the vertices geometry
 * \param[out] fnm  The Property Map that will be populated with the resulting
 *                  normal vectors.
 * \param[in]  gt   The Geometry Traits used to perfom the geometrical
 *                  calculations
 *
 * \ingroup  GenericManifoldFilters GenericNonManifoldFilters
 *
 * \note: this function name is prefixed with 'calculate_' to mimic
 *        CGAL usage of the term 'compute'.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename FaceNormalMap,
          typename GeometryTraits >
void
calculate_face_normals(const HalfedgeGraph &g,
                       const PointMap &pm,
                       FaceNormalMap fnm,
                       const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::face_iterator face_iterator;

  face_iterator fb, fe;
  for(boost::tie(fb, fe) = faces(g); fb != fe; ++fb)
  {
    Vector vec = FEVV::Operators::
        calculate_face_normal< HalfedgeGraph, PointMap, GeometryTraits >(
            *fb, g, pm, gt);
    put(fnm, *fb, vec);
  }
}


/**
 * \brief  Compute all faces normals
 *
 * \param  g     mesh
 * \param  pm    point map
 * \param  fnm   face normal map
 *
 * \ingroup  GenericManifoldFilters GenericNonManifoldFilters
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename FaceNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
calculate_face_normals(const HalfedgeGraph &g,
                       const PointMap &pm,
                       FaceNormalMap fnm)
{
  GeometryTraits gt(g);
  return calculate_face_normals< HalfedgeGraph,
                                 PointMap,
                                 FaceNormalMap,
                                 GeometryTraits >(g, pm, fnm, gt);
}

} // namespace Filters
} // namespace FEVV

#endif // CALCULATE_FACE_NORMAL_H
