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

#include <cmath> // M_PI


namespace FEVV {
namespace Math {

template< typename T >
inline T deg2rad(T deg)
{
  return static_cast< T >(deg * (M_PI / 180.f));
}

template< typename T >
inline T rad2deg(T rad)
{
  return static_cast< T >(rad * (180.f / M_PI));
}

} // namespace Math
} // namespace FEVV
