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
#include <boost/algorithm/clamp.hpp>
#include <algorithm>


inline
FEVV::Color::Color()
    : redComponent((unsigned char)0), greenComponent((unsigned char)0),
      blueComponent((unsigned char)0), alphaComponent((unsigned char)255)
{
}

inline
FEVV::Color::Color(const unsigned char _redValue,
                   const unsigned char _greenValue,
                   const unsigned char _blueValue,
                   const unsigned char _alphaValue)
    : redComponent(_redValue), greenComponent(_greenValue),
      blueComponent(_blueValue), alphaComponent(_alphaValue)
{
}

inline
FEVV::Color::Color(const unsigned char _grayValue,
                   const unsigned char _alphaValue)
    : redComponent(_grayValue), greenComponent(_grayValue),
      blueComponent(_grayValue), alphaComponent(_alphaValue)
{
}

inline
FEVV::Color::Color(const Color &_color)
    : redComponent(_color.redComponent), greenComponent(_color.greenComponent),
      blueComponent(_color.blueComponent), alphaComponent(_color.alphaComponent)
{
}

///////////////////////////////

inline
FEVV::Color &
FEVV::Color::setRGB(const unsigned char _redValue,
                    const unsigned char _greenValue,
                    const unsigned char _blueValue)
{
  redComponent = _redValue;
  greenComponent = _greenValue;
  blueComponent = _blueValue;
  alphaComponent = 255;

  return *this;
}

inline
FEVV::Color &
FEVV::Color::setRGBA(const unsigned char _redValue,
                     const unsigned char _greenValue,
                     const unsigned char _blueValue,
                     const unsigned char _alphaValue)
{
  redComponent = _redValue;
  greenComponent = _greenValue;
  blueComponent = _blueValue;
  alphaComponent = _alphaValue;

  return *this;
}

inline
void
FEVV::Color::red(const unsigned char _redValue)
{
  redComponent = _redValue;
}

inline
void
FEVV::Color::green(const unsigned char _greenValue)
{
  greenComponent = _greenValue;
}

inline
void
FEVV::Color::blue(const unsigned char _blueValue)
{
  blueComponent = _blueValue;
}

inline
void
FEVV::Color::alpha(const unsigned char _alphaValue)
{
  alphaComponent = _alphaValue;
}

inline
unsigned char
FEVV::Color::red() const
{
  return redComponent;
}

inline
unsigned char
FEVV::Color::green() const
{
  return greenComponent;
}

inline
unsigned char
FEVV::Color::blue() const
{
  return blueComponent;
}

inline
unsigned char
FEVV::Color::alpha() const
{
  return alphaComponent;
}

///////////////////////////////

inline
bool
FEVV::Color::operator==(const Color &_color) const
{
  return redComponent == _color.redComponent &&
         greenComponent == _color.greenComponent &&
         blueComponent == _color.blueComponent &&
         alphaComponent == _color.alphaComponent;
}

inline
bool
FEVV::Color::operator!=(const Color &_color) const
{
  return redComponent != _color.redComponent ||
         greenComponent != _color.greenComponent ||
         blueComponent != _color.blueComponent ||
         alphaComponent != _color.alphaComponent;
}

inline
bool
FEVV::Color::operator<(const Color &_color) const
{
  if(redComponent < _color.redComponent)
  {
    return true;
  }
  else if(redComponent == _color.redComponent)
  {
    if(greenComponent < _color.greenComponent)
    {
      return true;
    }
    else if(greenComponent == _color.greenComponent)
    {
      if(blueComponent < _color.blueComponent)
      {
        return true;
      }
      else if(blueComponent == _color.blueComponent)
      {
        return alphaComponent < _color.alphaComponent;
      }
    }
  }
  return false;
}

inline
FEVV::Color &
FEVV::Color::operator=(const Color &_color)
{
  redComponent = _color.redComponent;
  greenComponent = _color.greenComponent;
  blueComponent = _color.blueComponent;
  alphaComponent = _color.alphaComponent;

  return *this;
}

inline
FEVV::Color
FEVV::Color::operator+(const Color &_color)
{
  Color c;

  c.redComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)redComponent + (unsigned int)_color.redComponent, 0, 255);
  c.greenComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)greenComponent + (unsigned int)_color.greenComponent,
      0,
      255);
  c.blueComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)blueComponent + (unsigned int)_color.blueComponent, 0, 255);
  c.alphaComponent = (unsigned char)std::max(
      (unsigned int)alphaComponent, (unsigned int)_color.alphaComponent);

  return c;
}

inline
FEVV::Color &
FEVV::Color::operator+=(const Color &_color)
{
  redComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)redComponent + (unsigned int)_color.redComponent, 0, 255);
  greenComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)greenComponent + (unsigned int)_color.greenComponent,
      0,
      255);
  blueComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)blueComponent + (unsigned int)_color.blueComponent, 0, 255);
  alphaComponent = (unsigned char)std::max((unsigned int)alphaComponent,
                                           (unsigned int)_color.alphaComponent);

  return *this;
}

inline
FEVV::Color
FEVV::Color::operator-(const Color &_color)
{
  Color c;

  c.redComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)redComponent - (unsigned int)_color.redComponent, 0, 255);
  c.greenComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)greenComponent - (unsigned int)_color.greenComponent,
      0,
      255);
  c.blueComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)blueComponent - (unsigned int)_color.blueComponent, 0, 255);
  c.alphaComponent = (unsigned char)std::min(
      (unsigned int)alphaComponent, (unsigned int)_color.alphaComponent);

  return c;
}

inline
FEVV::Color &
FEVV::Color::operator-=(const Color &_color)
{
  redComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)redComponent - (unsigned int)_color.redComponent, 0, 255);
  greenComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)greenComponent - (unsigned int)_color.greenComponent,
      0,
      255);
  blueComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)blueComponent - (unsigned int)_color.blueComponent, 0, 255);
  alphaComponent = (unsigned char)std::min((unsigned int)alphaComponent,
                                           (unsigned int)_color.alphaComponent);

  return *this;
}

inline
FEVV::Color FEVV::Color::operator*(const double _scale)
{
  Color c;

  c.redComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)redComponent * _scale, 0, 255);
  c.greenComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)greenComponent * _scale, 0, 255);
  c.blueComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)blueComponent * _scale, 0, 255);
  c.alphaComponent = alphaComponent;

  return c;
}

inline
FEVV::Color &
FEVV::Color::operator*=(const double _scale)
{
  redComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)redComponent * _scale, 0, 255);
  greenComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)greenComponent * _scale, 0, 255);
  blueComponent = (unsigned char)boost::algorithm::clamp(
      (unsigned int)blueComponent * _scale, 0, 255);

  return *this;
}

///////////////////////////////

inline
void
FEVV::Color::selfDisplay(std::ostream &_out) const
{
  _out << "[Color] RGBA(" << (unsigned int)redComponent << ","
       << (unsigned int)greenComponent << "," << (unsigned int)blueComponent
       << "," << (unsigned int)alphaComponent << ")";
}

///////////////////////////////

inline
FEVV::Color FEVV::Color::None(void)
{
  return Color(0, 0, 0, 0);
}

inline
FEVV::Color FEVV::Color::Black(void)
{
  return Color(0, 0, 0);
}

inline
FEVV::Color FEVV::Color::White(void)
{
  return Color(255, 255, 255);
}

inline
FEVV::Color FEVV::Color::Red(void)
{
  return Color(244, 67, 54);
}

inline
FEVV::Color FEVV::Color::Pink(void)
{
  return Color(233, 30, 99);
}

inline
FEVV::Color FEVV::Color::Purple(void)
{
  return Color(156, 39, 176);
}

inline
FEVV::Color FEVV::Color::DeepPurple(void)
{
  return Color(103, 58, 183);
}

inline
FEVV::Color FEVV::Color::Indigo(void)
{
  return Color(63, 81, 181);
}

inline
FEVV::Color FEVV::Color::Blue(void)
{
  return Color(33, 150, 243);
}

inline
FEVV::Color FEVV::Color::LightBlue(void)
{
  return Color(3, 169, 244);
}

inline
FEVV::Color FEVV::Color::Cyan(void)
{
  return Color(0, 188, 212);
}

inline
FEVV::Color FEVV::Color::Teal(void)
{
  return Color(0, 150, 136);
}

inline
FEVV::Color FEVV::Color::Green(void)
{
  return Color(76, 175, 80);
}

inline
FEVV::Color FEVV::Color::LightGreen(void)
{
  return Color(139, 195, 74);
}

inline
FEVV::Color FEVV::Color::Lime(void)
{
  return Color(205, 220, 57);
}

inline
FEVV::Color FEVV::Color::Yellow(void)
{
  return Color(255, 235, 59);
}

inline
FEVV::Color FEVV::Color::Amber(void)
{
  return Color(255, 193, 7);
}

inline
FEVV::Color FEVV::Color::Orange(void)
{
  return Color(243, 156, 18);
}

inline
FEVV::Color FEVV::Color::DeepOrange(void)
{
  return Color(255, 87, 34);
}

inline
FEVV::Color FEVV::Color::Brown(void)
{
  return Color(121, 85, 72);
}

inline
FEVV::Color FEVV::Color::Grey(void)
{
  return Color(158, 158, 158);
}

inline
FEVV::Color FEVV::Color::BlueGrey(void)
{
  return Color(96, 125, 139);
}

///////////////////////////////


inline
FEVV::Color FEVV::Color::Turquoise(void)
{
  return Color(26, 188, 156);
}

inline
FEVV::Color FEVV::Color::GreenSea(void)
{
  return Color(22, 160, 133);
}

inline
FEVV::Color FEVV::Color::Emerald(void)
{
  return Color(46, 204, 113);
}

inline
FEVV::Color FEVV::Color::Nephritis(void)
{
  return Color(39, 174, 96);
}

inline
FEVV::Color FEVV::Color::PeterRiver(void)
{
  return Color(41, 128, 185);
}

inline
FEVV::Color FEVV::Color::BelizeHole(void)
{
  return Color(41, 128, 185);
}

inline
FEVV::Color FEVV::Color::Amethyst(void)
{
  return Color(155, 89, 182);
}

inline
FEVV::Color FEVV::Color::Wisteria(void)
{
  return Color(142, 68, 173);
}

inline
FEVV::Color FEVV::Color::WetAsphalt(void)
{
  return Color(52, 73, 94);
}

inline
FEVV::Color FEVV::Color::MidnightBlue(void)
{
  return Color(44, 62, 80);
}

inline
FEVV::Color FEVV::Color::SunFlower(void)
{
  return Color(241, 196, 15);
}

inline
FEVV::Color FEVV::Color::Carrot(void)
{
  return Color(230, 126, 34);
}

inline
FEVV::Color FEVV::Color::Pumpkin(void)
{
  return Color(211, 84, 0);
}

inline
FEVV::Color FEVV::Color::Alizarin(void)
{
  return Color(231, 76, 60);
}

inline
FEVV::Color FEVV::Color::Pomegranate(void)
{
  return Color(192, 57, 43);
}

inline
FEVV::Color FEVV::Color::Clouds(void)
{
  return Color(236, 240, 241);
}

inline
FEVV::Color FEVV::Color::Silver(void)
{
  return Color(189, 195, 199);
}

inline
FEVV::Color FEVV::Color::Concrete(void)
{
  return Color(149, 165, 166);
}

inline
FEVV::Color FEVV::Color::Asbestos(void)
{
  return Color(127, 140, 141);
}

