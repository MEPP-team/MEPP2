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

#include <cmath> // M_PI, asin

#define DEG2RAD(deg) (deg * (M_PI / 180.f))
#define RAD2DEG(rad) (rad * (180.f / M_PI))

namespace FEVV {
namespace Operators {

namespace Geometry {

/**
 * \brief  Safe call to the std::asin function.
 *         
 * \param[in] sine The sinus value.  
 * \return The arcsin (in rad) of sine value, sine value being 
 *         truncated to [-1,1] range to avoid input domain error.
 */
template< class T >
static T
asin(T sine)
{
  if(sine >= 1.f)
    return M_PI * 0.5f;
  else if(sine <= -1.f)
    return -M_PI * 0.5f;
  else
    return std::asin(sine);
}

/**
 * \brief  Safe call to the std::asin function.
 *         
 * \param[in] sine The sinus value.  
 * \return The arcsin (in degree) of sine value, sine value being 
 *         truncated to [-1,1] range to avoid input domain error.
 */
template< class T >
static T
asin_degree(T sine)
{
  return RAD2DEG(Geometry::asin< T >(sine));
}

/**
 * \brief  Safe call to the std::acos function.
 *         
 * \param[in] cosine The cosine value.  
 * \return The arccos (in rad) of cosine value, cosine value being 
 *         truncated to [-1,1] range to avoid input domain error.
 */
template< class T >
inline T
acos(T cosine)
{
  if(cosine < -1.f)
    cosine = -1.f;
  else if(cosine > 1.f)
    cosine = 1.f;

  return std::acos(cosine);
}

/**
 * \brief  Safe call to the std::acos function.
 *         
 * \param[in] cosine The cosine value.  
 * \return The arccos (in degree) of cosine value, cosine value being 
 *         truncated to [-1,1] range to avoid input domain error.
 */
template< class T >
inline T
acos_degree(T cosine)
{
  return RAD2DEG(Geometry::acos< T >(cosine));
}
} // namespace Geometry

} // namespace Operators
} // namespace FEVV
