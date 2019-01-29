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
add_p(const typename GeometryTraits::Point &p,
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
sub_p(const typename GeometryTraits::Point &p,
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
sub(const typename GeometryTraits::Point &p1,
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

/**
 * \brief  Checks that geometry traits vectors have a size of 3.
 */
template< typename GeometryTraits >
int
verify_assumption(GeometryTraits gt)
{
  typedef typename GeometryTraits::Vector Vector;
  Vector v;

  if(v.size() != 3)
  {
    std::cout << "Dimension of ambiant space vector is " << v.size()
              << " instead of the expected 3" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template< typename GeometryTraits >
int
verify_operator_results(GeometryTraits gt)
{
  typedef typename GeometryTraits::Vector Vector;

  if(gt.dot_product(Vector(5, 8, 4), Vector(10, -2, -3)) != 22)
  {
    std::cout << "Erroneous result for DotProduct. Expected " << 22
              << " and instead got "
              << gt.dot_product(Vector(5, 8, 4), Vector(10, -2, -3))
              << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


} // namespace GeometryTraits
} // namespace FEVV

