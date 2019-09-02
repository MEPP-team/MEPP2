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
#include <CGAL/boost/graph/iterator.h> // for circulators
#include <cmath> // for std::sqrt()

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/Geometry/triangle_rad_angle.hpp"
#include "FEVV/Operators/Generic/is_geometrical_fold.hpp"
#include "FEVV/Operators/Generic/Manifold/calculate_vertex_normal.hpp"


namespace FEVV {
namespace Operators {

/**
 * \brief  Angle-based smoothing of the vertex' one ring of vertex positions.
 *         Increases the min angle and decreases the max angle
 *
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<> specialization.
 * \tparam  PointMap A modifiable point map to manage vertex positions.
 * \tparam  GeometryTraits The geometric kernel when available. This is
 *          defaulted to FEVV::Geometry_traits<FaceGraph>.
 * \param[in] v The vertex whose one-ring angle-based centroid is computed.
 * \param[in] g The FaceGraph instance from which the v vertex will be taken.
 * \param[in] pm The point map which associate a vertex descriptor with a vertex
 *            position.
 * \param[in] smoothing_factor A (usually positive) factor used to compute the
 *            position of the returned centroid following
 *            V_pos + SMOOTHING_FACTORx(C_computed - V_pos).
 * \param[in] gt The geometry trait object.
 * \return The angles-based centroid of the vertex v. 
 * \pre     g must be a triangular mesh.
 */
template< typename FaceGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
typename boost::property_traits< PointMap >::value_type
vertex_one_ring_angles_based_centroid(
    typename boost::graph_traits< FaceGraph >::vertex_descriptor v,
    const FaceGraph &g,
    const PointMap &pm,
    const typename GeometryTraits::Scalar smoothing_factor,
    const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename GeometryTraits::Point Point;
  // if ( boost::optional<typename
  // boost::graph_traits<FaceGraph>::halfedge_descriptor> h = CGAL::is_border(v,
  // g) )
  //{
  //	return get(pm, v); // for the time being do not change border vertex
  //}
  //////////////////////////////////////////////////////////////////////////////////
  Vector vec = gt.NULL_VECTOR;
  double sum_iangles = 0.0f, alpha1, alpha2, iangle;
  // std::cout << "around v " << get(pm, v) << std::endl;
  CGAL::Halfedge_around_target_circulator< FaceGraph > cir_he(v, g),
      cir_he_end(cir_he);
  do
  {
    // The angle must be divided in two to take convex and concave cases into
    // account
    alpha1 = FEVV::Operators::Geometry::triangle_rad_angle< GeometryTraits >(
        get(pm, target(*cir_he, g)),
        get(pm, source(*cir_he, g)),
        get(pm, source(prev(opposite(*cir_he, g), g), g)),
        gt);
    // std::cout << "t1 = " << get(pm, target(*cir_he, g)) << "; " << get(pm,
    // source(*cir_he, g)) << "; " << get(pm, source(prev(opposite(*cir_he, g),
    // g), g)) << std::endl;
    alpha2 = FEVV::Operators::Geometry::triangle_rad_angle< GeometryTraits >(
        get(pm, target(next(*cir_he, g), g)),
        get(pm, source(*cir_he, g)),
        get(pm, target(*cir_he, g)),
        gt);
    // std::cout << "t2 = " << get(pm, target(next(*cir_he, g), g)) << "; " <<
    // get(pm, source(*cir_he, g)) << "; " << get(pm, target(*cir_he, g)) <<
    // std::endl; std::cout << "alpha1 = " << alpha1 << " alpha2 = " << alpha2 <<
    // std::endl;
    if(alpha1 < 0.0 || alpha2 < 0.0)
      assert(false);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Vector d1, d2, c, dc;
    d1 = gt.sub(get(pm, target(next(opposite(*cir_he, g), g), g)),
                get(pm, source(*cir_he, g)));
    // d1 = pm[target(next(opposite(*cir_he, g), g), g)] - pm[source(*cir_he,
    // g)]; // more warnings
    d1 = gt.normalize(d1);
    d2 = gt.sub(get(pm, target(next(*cir_he, g), g)),
                get(pm, source(*cir_he, g)));
    d2 = gt.normalize(d2);
    c = 0.5 * (d1 + d2);

    if(FEVV::Operators::
           is_geometrical_fold< FaceGraph, PointMap, GeometryTraits >(
               prev(prev(opposite(*cir_he, g), g), g),
                 // prev prev instead of next to handle polygonal cases
               next(*cir_he, g),
               g,
               pm,
               gt))
      continue;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(gt.dot_product(c, c) > MACH_EPS_DOUBLE)
    {
      c = c * (1.0 / std::sqrt(gt.dot_product(c, c)));
      Point new_ideal_point = gt.add_p(
          get(pm, source(*cir_he, g)),
          c * std::sqrt(gt.dot_product(gt.sub(get(pm, target(*cir_he, g)),
                                              get(pm, source(*cir_he, g))),
                                       gt.sub(get(pm, target(*cir_he, g)),
                                              get(pm, source(*cir_he, g))))));

      // Calculate the difference between two adjacent angles...
      // Small angles (alpha1+alpha2) must be favorized to avoid foldover
      iangle = alpha1 + alpha2;
      iangle = 1.0 / (iangle * iangle);
      dc = gt.sub(new_ideal_point, get(pm, v));
      // normalize_Vector(dc); // unit direction towards the bissector of that
      // angle
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // the current best results:
      vec = vec + dc * iangle;
      sum_iangles += iangle;
    }
  } while(++cir_he != cir_he_end);
  if(::fabs(sum_iangles) > MACH_EPS_DOUBLE)
    vec = vec * (1.0 / sum_iangles);

  Vector n;
#if 0
		// MORE FIDELITY!!!!!
		// we do that to not change the NRJ being minimized
		N = (pVertex->closestObs())->normal();
      // to be more sensitive to the initial curvature
		//N = CGAL::cross_product(pVertex->closestObs()->VKmin,
    //    pVertex->closestObs()->VKmax); // too much geometric error
#endif
  // more geometric error, since we use the initial surface as reference
  n = FEVV::Operators::
      calculate_vertex_normal< FaceGraph, PointMap, GeometryTraits >(
          v, g, pm, gt); // to avoid introducing noise in flat region

  vec = vec - gt.dot_product(vec, n) * n;
  // std::cout << "vec = (" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")"
  // << std::endl;
  return gt.add_p(get(pm, v), smoothing_factor * vec);
}

} // namespace Operators
} // namespace FEVV
