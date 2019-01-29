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

#include "FEVV/Tools/Math/MatrixOperations.hpp"

namespace FEVV {
namespace Operators {

namespace Geometry {

 /**
 * \brief  Compute the intersection of a sphere (center + radius)
 *         with a ray/line (starting point + direction vector).
 *         Sphere equation is given by (x-Cx)^2+(y-Cy)^2+(z-Cz)^2=r^2
 *         Line equation is given by (x, y, z) = p + txV
 * \tparam  GeometryTraits The geometric kernel. 
 * \param[in] center The sphere center.
 * \param[in] r The sphere radius.
 * \param[in] p The starting point/origin of the ray/line.
 * \param[in,out] v The ray/line direction to clip by the sphere when 
 *                an intersection occurs in the direction of v.
 * \param[in] gt The geometry trait object.			  
 * \return True if there is an intersection between the sphere and
 *         the ray starting at p position and going along the v direction,
 *         not farther than the v' length.
 * \pre    p must be inside the sphere. 
 */	
template< typename GeometryTraits >
static bool
sphere_clip_vector(
    const typename GeometryTraits::Point &center, 
    double                               r,                                     
    const typename GeometryTraits::Point &p,
    typename GeometryTraits::Vector      &v,
    const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  if(r < 0)
    r = -r; // the radius cannot be negative

  Vector w = gt.sub(p, center); // direction towards P from the sphere center
  double a = gt.dot_product(v, v);

  if(fabs(a) < std::numeric_limits< double >::epsilon())
    return false;

  double b = 2.0 * gt.dot_product(v, w);
  double c = gt.dot_product(w, w) - r * r;
  double delta = b * b - 4. * a * c;
  if(delta < 0.)
  {
    // Should not happen when p is inside the sphere, but happens sometimes 
	// (numerical precision)
    return true;
  }

  double t1 = (-b + ::sqrt(delta)) /
              (2.0 * a) /*, t2 = (-b - ::sqrt(delta)) / (2.0 * a)*/;

  if(t1 >= 1.)
  {
    // Inside the sphere
    return false;
  }

  // if (t1 < 0.) {
  // Should not happen, but happens sometimes (numerical precision or P not
  // inside the sphere)
  //	return true;
  //}

  v = gt.scalar_mult(v, t1);

  return true;
}
} // namespace Geometry

} // namespace Operators
} // namespace FEVV
