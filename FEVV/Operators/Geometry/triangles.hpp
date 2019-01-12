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

#include "FEVV/Operators/Geometry/AngleOperations.hpp"
// For are_collinear() and MACH_EPS_DOUBLE :
#include "FEVV/Tools/Math/MatrixOperations.hpp"

namespace FEVV {
namespace Operators {
namespace Geometry {
/**
 * \brief   Compute the barycenter/mean position of a triangle (given by 3
 * points).
 *
 * \tparam GeometryTraits The geometric kernel.
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] gt The geometry trait object.
 * \return The barycenter/mean position of the 3 input points (Point).
 */
template< typename GeometryTraits >
inline typename GeometryTraits::Point
triangle_barycenter(const typename GeometryTraits::Point &a,
                    const typename GeometryTraits::Point &b,
                    const typename GeometryTraits::Point &c,
                    const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Point Point;

  return Point((a[0] + b[0] + c[0]) / 3,
               (a[1] + b[1] + c[1]) / 3,
               (a[2] + b[2] + c[2]) / 3);
}
/**
 * \brief   Compute the perimeter of a triangle (given by 3 points).
 *
 * \tparam GeometryTraits The geometric kernel.
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] gt The geometry trait object.
 * \return The triangle perimeter (double Scalar).
 */
template< typename GeometryTraits >
inline double
triangle_perimeter(const typename GeometryTraits::Point &a,
                   const typename GeometryTraits::Point &b,
                   const typename GeometryTraits::Point &c,
                   const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  Vector ab(b - a), bc(c - b), ca(a - c);

  return (gt.length(ab) + gt.length(bc) + gt.length(ca));
}
/**
 * \brief   Compute the unnormalized normal of a triangle (given by 3 points).
 *
 * \tparam GeometryTraits The geometric kernel.
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] gt The geometry trait object.
 * \return The triangle unnormalized normal (Vector).
 */
template< typename GeometryTraits >
inline typename GeometryTraits::Vector
triangle_normal_unnormalized(const typename GeometryTraits::Point &a,
                             const typename GeometryTraits::Point &b,
                             const typename GeometryTraits::Point &c,
                             const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  Vector cb(b[0] - c[0], b[1] - c[1], b[2] - c[2]),
      ca(a[0] - c[0], a[1] - c[1], a[2] - c[2]);
  return gt.cross_product(ca, cb);
}
/**
 * \brief   Compute the area of a triangle (given by 3 points).
 *
 * \tparam GeometryTraits The geometric kernel.
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param gt The geometry trait object.
 * \return The triangle area (double Scalar).
 */
template< typename GeometryTraits >
inline double
triangle_area(const typename GeometryTraits::Point &a,
              const typename GeometryTraits::Point &b,
              const typename GeometryTraits::Point &c,
              const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  Vector n = triangle_normal_unnormalized< GeometryTraits >(a, b, c, gt);
  return gt.length(n) * 0.5;
}
/**
 * \brief   Compute the shape potential of a triangle (given by 3 points).
 *          Triangle "equilateralness" measure.
 *
 * \tparam  FaceGraph  a Mesh type that provides a Model of the FaceGraph
 *                     Concept through a boost::graph_traits<> specialization.
 * \tparam  GeometryTraits  The geometric kernel when available. This is
                            defaulted to FEVV::Geometry_traits<FaceGraph>.
 * \param[in]   a          The first point.
 * \param[in]   b          The second point.
 * \param[in]   c          The third point.
 * \param[in]   PROHIB_ENERGY  The maximum authorized value, e.g. returns that
                           value if smaller than the computed shape potential.
                           This parameter is needed because when we compute the
                           mesh quality using the sum of shape potential of
                           each triangle, the sum may overflow.
 * \param[in]   gt         The geometry trait object.
 * \return  The triangle area (double Scalar).
 */
template< typename FaceGraph,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
inline double
triangle_shape_potential(
    const typename GeometryTraits::Point &a,
    const typename GeometryTraits::Point
        &b, // the central point for which we give the gradient
    const typename GeometryTraits::Point &c,
    const double prohib_energy,
    const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  double d[3], circumradius, potential = prohib_energy, den;

  Vector ab(b[0] - a[0], b[1] - a[1], b[2] - a[2]),
      cb(b[0] - c[0], b[1] - c[1], b[2] - c[2]),
      ca(a[0] - c[0], a[1] - c[1], a[2] - c[2]);

  // Length of the triangle sides
  d[0] = gt.length(ab);
  d[1] = gt.length(ca);
  d[2] = gt.length(cb);

  den = std::min< double >(d[2], std::min< double >(d[0], d[1]));
  if((den > 1e-8) &&
     !FEVV::Math::Vector::are_collinear< GeometryTraits >(ab, ca))
  {
    circumradius = d[0] * d[1] * d[2] /
                   (4.0 * triangle_area< GeometryTraits >(a, b, c, gt));
    potential = circumradius / den;
  }

  if(potential > prohib_energy)
    potential =
        prohib_energy; // to truncate the shape error to a managable max value
  // std::cout << potential << std::endl;
  return potential;
}

} // namespace Geometry
} // namespace Operators
} // namespace FEVV
