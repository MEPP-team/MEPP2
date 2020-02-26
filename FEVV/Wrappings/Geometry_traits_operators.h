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

#include <iostream> // due to some std::cout

namespace FEVV {
namespace GeometryTraits {

/**
 * \brief  Returns the sum of vectors V1 and V2.
 */
template< typename GeometryTraits >
typename GeometryTraits::Vector
add_v(const typename GeometryTraits::Vector &v1,
      const typename GeometryTraits::Vector &v2)
{
  typedef typename GeometryTraits::Vector Vector;
  Vector result(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
  return result;
};

/**
 * \brief  Returns the sum of point P and vector V.
 */
template< typename GeometryTraits >
typename GeometryTraits::Point
add_pv(const typename GeometryTraits::Point &p,
       const typename GeometryTraits::Vector &v)
{
  typedef typename GeometryTraits::Point Point;
  Point result(p[0] + v[0], p[1] + v[1], p[2] + v[2]);
  return result;
};

/**
 * \brief  Returns point P minus vector V.
 */
template< typename GeometryTraits >
typename GeometryTraits::Point
sub_pv(const typename GeometryTraits::Point &p,
       const typename GeometryTraits::Vector &v)
{
  typedef typename GeometryTraits::Point Point;
  Point result(p[0] - v[0], p[1] - v[1], p[2] - v[2]);
  return result;
};

/**
 * \brief  Returns point P1 minus point P2.
 */
template< typename GeometryTraits >
typename GeometryTraits::Vector
sub_p(const typename GeometryTraits::Point &p1,
      const typename GeometryTraits::Point &p2)
{
  typedef typename GeometryTraits::Vector Vector;
  Vector result(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
  return result;
};

/**
 * \brief  Returns the dot product of vectors V1 and V2.
 */
template< typename GeometryTraits >
typename GeometryTraits::Scalar
dot_product(const typename GeometryTraits::Vector &v1,
            const typename GeometryTraits::Vector &v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
};

} // namespace GeometryTraits
} // namespace FEVV

