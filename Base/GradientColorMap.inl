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
#include <cstdlib>
#include <vector>

template< typename PValue >
inline FEVV::GradientColorMap< PValue >::GradientColorMap(
    const Value &_minV,
    const Value &_maxV,
    const Color &_firstColor,
    const Color &_lastColor)
    : myMin(_minV), myMax(_maxV)
{
  // ASSERT_MSG(myMin < myMax, "Max should be strictly greather than Min in a
  // colormap.");
  if(_firstColor != Color::None() && _lastColor != Color::None())
  {
    myColors.push_back(_firstColor);
    myColors.push_back(_lastColor);
  }
}

template< typename PValue >
inline FEVV::GradientColorMap< PValue >::GradientColorMap(
    const GradientColorMap< Value > &_other)
    : myMin(_other.myMin), myMax(_other.myMax), myColors(_other.myColors)
{
  // ASSERT_MSG(myMin < myMax, "Max should be strictly greather than Min in a
  // colormap.");
}

template< typename PValue >
FEVV::GradientColorMap< PValue > &
FEVV::GradientColorMap< PValue >::
operator=(const FEVV::GradientColorMap< Value > &_other)
{
  if(&_other != this)
  {
    myMin = _other.myMin;
    myMax = _other.myMax;
    myColors = _other.myColors;
    // ASSERT_MSG(myMin < myMax, "Max should be strictly greather than Min in a
    // colormap.");
  }
  return *this;
}

template< typename PValue >
inline FEVV::Color
FEVV::GradientColorMap< PValue >::operator()(const Value &_value) const
{
  return GradientColorMap< Value >::getColor(myColors, myMin, myMax, _value);
}

template< typename PValue >
inline void
FEVV::GradientColorMap< PValue >::addColor(const Color &_color)
{
  myColors.push_back(_color);
}

template< typename PValue >
inline void
FEVV::GradientColorMap< PValue >::clear()
{
  myColors.clear();
}

template< typename PValue >
inline const PValue &
FEVV::GradientColorMap< PValue >::minValue() const
{
  return myMin;
}

template< typename PValue >
inline const PValue &
FEVV::GradientColorMap< PValue >::maxValue() const
{
  return myMax;
}

template< typename PValue >
inline FEVV::Color
FEVV::GradientColorMap< PValue >::getColor(const std::vector< Color > &_colors,
                                           const Value &_min,
                                           const Value &_max,
                                           const Value &_value)
{
  // ASSERT_MSG(min < max, "Max should be strictly greather than Min in a
  // colormap.");
  if(_colors.size() < 2)
  {
    return Color::None();
  }

  double scale = static_cast< double >(_value - _min) / (_max - _min);
  const int intervals = (const int)_colors.size() - 1;
  int upper_index = static_cast< int >(ceil(intervals * scale));
  if(!upper_index) // Special case when value == min.
  {
    upper_index = 1;
  }
  const Color &firstColor = _colors[upper_index - 1];
  const Color &lastColor = _colors[upper_index];
  scale = (scale * intervals) - (upper_index - 1);

  const unsigned char red = static_cast< unsigned char >(
      firstColor.red() + scale * (lastColor.red() - firstColor.red()));
  const unsigned char green = static_cast< unsigned char >(
      firstColor.green() + scale * (lastColor.green() - firstColor.green()));
  const unsigned char blue = static_cast< unsigned char >(
      firstColor.blue() + scale * (lastColor.blue() - firstColor.blue()));
  return Color(red, green, blue);
}
