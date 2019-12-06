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
 * \brief   Compute the barycenter/mean position of a quad (given by 4
 * points).
 *
 * \tparam GeometryTraits The geometric kernel.
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] d The fourth point. 
 * \param[in] gt The geometry trait object.
 * \return The barycenter/mean position of the 4 input points (Point).
 */
template< typename GeometryTraits >
inline typename GeometryTraits::Point
quad_barycenter(const typename GeometryTraits::Point &a,
                const typename GeometryTraits::Point &b,
                const typename GeometryTraits::Point &c,
                const typename GeometryTraits::Point &d,				
                const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Point Point;

  return Point((a[0] + b[0] + c[0] + d[0]) / 4,
               (a[1] + b[1] + c[1] + d[1]) / 4,
               (a[2] + b[2] + c[2] + d[2]) / 4);
}
/**
 * \brief   Compute the perimeter of a quad (given by 4 points).
 *
 * \tparam GeometryTraits The geometric kernel.
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] d The fourth point.  
 * \param[in] gt The geometry trait object.
 * \return The quad perimeter (double Scalar).
 */
template< typename GeometryTraits >
inline double
quad_perimeter(const typename GeometryTraits::Point &a,
               const typename GeometryTraits::Point &b,
               const typename GeometryTraits::Point &c,
               const typename GeometryTraits::Point &d,			   
               const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  Vector ab(b - a), bc(c - b), cd(d - c), da(a - d);

  return (gt.length(ab) + gt.length(bc) + gt.length(cd)) + gt.length(da));
}

/**
 * \brief   Compute the area of a quad (given by 4 points).
 *
 * \tparam GeometryTraits The geometric kernel.
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] d The fourth point. 
 * \param gt The geometry trait object.
 * \return The quad area (double Scalar).
 * \remark This works for any side length.
 */
template< typename GeometryTraits >
inline double
quad_area(const typename GeometryTraits::Point &a,
          const typename GeometryTraits::Point &b,
          const typename GeometryTraits::Point &c,
          const typename GeometryTraits::Point &d,		  
          const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  Vector n1 = triangle_normal_unnormalized< GeometryTraits >(a, b, c, gt);
  Vector n2 = triangle_normal_unnormalized< GeometryTraits >(c, d, a, gt);
  return (gt.length(n1) + gt.length(n2)) * 0.5;
}
/**
 * \brief   Compute the shape potential of a quad (given by 4 points).
 *          Quad "quality" measure.
 *
 * \tparam  FaceGraph  a Mesh type that provides a Model of the FaceGraph
 *                     Concept through a boost::graph_traits<> specialization.
 * \tparam  GeometryTraits  The geometric kernel when available. This is
                            defaulted to FEVV::Geometry_traits<FaceGraph>.
 * \param[in]   a          The first point.
 * \param[in]   b          The second point.
 * \param[in]   c          The third point.
 * \param[in]   d          The fourth point.  
 * \param[in]   PROHIB_ENERGY  The maximum authorized value, e.g. returns that
                           value if smaller than the computed shape potential.
                           This parameter is needed because when we compute the
                           mesh quality using the sum of shape potential of
                           each face, the sum may overflow.
 * \param[in]   gt         The geometry trait object.
 * \return  The shape potential (double Scalar).
 */
template< typename FaceGraph,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
inline double
quad_shape_potential(
    const typename GeometryTraits::Point &a,
    const typename GeometryTraits::Point
        &b, // the central point for which we give the gradient
    const typename GeometryTraits::Point &c,
    const typename GeometryTraits::Point &d,	
    const double prohib_energy,
    const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  double d[5], circumradius, potential = prohib_energy, den;

  Vector ab(b[0] - a[0], b[1] - a[1], b[2] - a[2]),
         cb(b[0] - c[0], b[1] - c[1], b[2] - c[2]),
         ca(a[0] - c[0], a[1] - c[1], a[2] - c[2]),
		 cd(d[0] - c[0], d[1] - c[1], d[2] - c[2]),
		 da(a[0] - d[0], a[1] - d[1], a[2] - d[2]);

  // Length of the sides
  d[0] = gt.length(ab);
  d[1] = gt.length(ca);
  d[2] = gt.length(cb);
  d[3] = gt.length(cd);  
  d[4] = gt.length(da);   

  assert(false);// not finished yet

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
