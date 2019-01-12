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

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/Geometry/AngleOperations.hpp"
#include "FEVV/Tools/Math/MatrixOperations.hpp" // MACH_EPS_DOUBLE

namespace FEVV {
namespace Operators {
namespace Geometry {
/**
 * \brief   Compute the angle of a triangle (given by 3 points).
 *
 * \tparam GeometryTraits The geometric kernel.
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] gt The geometry trait object.
 * \return The triangle angle in rad (Scalar).
 */
template< typename GeometryTraits >
inline typename GeometryTraits::Scalar
triangle_rad_angle(const typename GeometryTraits::Point &a,
                   const typename GeometryTraits::Point
                       &b, // the central point on which we compute the angle
                   const typename GeometryTraits::Point &c,
                   const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename GeometryTraits::Scalar Scalar;
  Vector ba(a[0] - b[0], a[1] - b[1], a[2] - b[2]),
      bc(c[0] - b[0], c[1] - b[1], c[2] - b[2]);
  Scalar lba = gt.length(ba), lbc = gt.length(bc);

  if(lba > MACH_EPS_DOUBLE)
    ba = ba / lba;
  else
    ba = gt.NULL_VECTOR;
  if(lbc > MACH_EPS_DOUBLE)
    bc = bc / lbc;
  else
    bc = gt.NULL_VECTOR;

  Scalar cosabc = gt.dot_product(ba, bc);

  return FEVV::Operators::Geometry::acos(cosabc);
}

} // namespace Geometry
} // namespace Operators
} // namespace FEVV
