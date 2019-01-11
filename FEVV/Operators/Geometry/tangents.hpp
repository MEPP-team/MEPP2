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

#include <cmath>

namespace FEVV {
namespace Operators {

/**
 * \brief  Calculate the actual face tangent from three connected vertices'
 * positions & texture coordinates.
 *
 * \param[in]  pt1                 first vertex's position.
 * \param[in]  pt2                 second vertex's position.
 * \param[in]  pt3                 third vertex's position.
 * \param[in]  uv1                 first vertex's texture coordinates.
 * \param[in]  uv2                 second vertex's texture coordinates.
 * \param[in]  uv3                 third vertex's texture coordinates.
 *
 * \return  calculated face tangent.
 */
template< typename Point, typename Vector >
Vector
calculate_triangle_tangent(const Point &pt1,
                           const Point &pt2,
                           const Point &pt3,
                           const Vector &uv1,
                           const Vector &uv2,
                           const Vector &uv3)
{
  const Vector pos1(pt1[0], pt1[1], pt1[2]);
  const Vector pos2(pt2[0], pt2[1], pt2[2]);
  const Vector pos3(pt3[0], pt3[1], pt3[2]);

  const Vector first_edge = pos2 - pos1;
  const Vector second_edge = pos3 - pos1;

  Vector first_uv_diff = uv2 - uv1;
  Vector second_uv_diff = uv3 - uv1;

  // RM: fix for black areas, which happens when UVs are equal;
  // manually forcing UVs to be different from each other
  if(std::abs(first_uv_diff[0]) <= 0.01f && std::abs(first_uv_diff[1]) <= 0.01f)
  {
    const Vector uv2_fix(uv2[0] + 0.0001, uv2[1], uv2[2]);
    first_uv_diff = uv2_fix - uv1;
  }

  if(std::abs(second_uv_diff[0]) <= 0.01f &&
     std::abs(second_uv_diff[1]) <= 0.01f)
  {
    const Vector uv3_fix(uv3[0], uv3[1] + 0.0001, uv3[2]);
    second_uv_diff = uv3_fix - uv1;
  }

  // RM: calculate tangent from positions & UV space
  const float inversion_factor = 1.f / (first_uv_diff[0] * second_uv_diff[1] -
                                        second_uv_diff[0] * first_uv_diff[1]);

  const Vector tangent =
      (first_edge * second_uv_diff[1] - second_edge * first_uv_diff[1]) *
      inversion_factor;

  return tangent;
}

} // namespace Operators
} // namespace FEVV
