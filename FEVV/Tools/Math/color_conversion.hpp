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

#include <cstdint> // for uint8_t


namespace FEVV {
namespace Math {


// Note: it is assumed that color value is either a 8 bits unsigned int value in
// range [0; 255], or a floating point value in range [0.0; 1.0].


/**
 * Convert color value given in [0; 255] to another type.
 * General case, no range change, casting is enough.
 */
template< typename OutT >
OutT convert_color_value(uint8_t color_0_255)
{
  return static_cast< OutT >(color_0_255);
}

/**
 * Convert color value in [0; 255] to color value in [0.0; 1.0] (float).
 */
template<>
inline
float convert_color_value< float >(uint8_t color_0_255)
{
  return color_0_255 / 255.0f;
}

/**
 * Convert color value in [0; 255] to color value in [0.0; 1.0] (double).
 */
template<>
inline
double convert_color_value< double >(uint8_t color_0_255)
{
  return color_0_255 / 255.0;
}

/**
 * Convert color value given in [0.0; 1.0] to another type.
 * General case, no range change, casting is enough.
 */
template< typename OutT >
OutT convert_color_value(double color_0_1)
{
  return static_cast< OutT >(color_0_1);
}

/**
 * Convert color value in [0.0; 1.0] to color value in [0; 255].
 */
template<>
inline
uint8_t convert_color_value< uint8_t >(double color_0_1)
{
  // note: applies also to color_0_1 of type float
  return static_cast< uint8_t >(color_0_1 * 255.0);
}


} // namespace Math
} // namespace FEVV
