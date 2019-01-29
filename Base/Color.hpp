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

#include <iostream>

namespace FEVV {

class Color
{
public:
  /**
   * Destructor.
   */
  ~Color() {}

  /**
   * Default constructor.
   * @note Set the color to BLACK.
   */
  Color();

  /**
   * Constructor.
   *
   * @param[in]    _redValue     the red component. Value must be into [0,255]
   * interval.
   * @param[in]    _greenValue   the green component. Value must be into [0,255]
   * interval.
   * @param[in]    _blueValue    the blue component. Value must be into [0,255]
   * interval.
   * @param[in]    _alphaValue   the alpha transparency. Value must be into
   * [0,255] interval. (Default value = 255)
   */
  Color(const unsigned char _redValue,
        const unsigned char _greenValue,
        const unsigned char _blueValue,
        const unsigned char _alphaValue = 255);

  /**
   * Constructor.
   *
   * @param[in]    _grayValue    the grey component. Value must be into [0,255]
   * interval.
   * @param[in]    _alphaValue   the alpha transparency. Value must be into
   * [0,255] interval. (Default value = 255)
   */
  Color(const unsigned char _grayValue, const unsigned char _alphaValue = 255);

  /**
   * Copy Constructor.
   *
   * @param[in] _color the color to copy.
   */
  Color(const Color &_color);

  ///////////////////////////////

  /**
   * Set the RGB value to the current Color.
   *
   * @note The alpha component will be set to 255.
   *
   * @param[in]    _redValue     the red component. Value must be into [0,255]
   * interval.
   * @param[in]    _greenValue   the green component. Value must be into [0,255]
   * interval.
   * @param[in]    _blueValue    the blue component. Value must be into [0,255]
   * interval.
   */
  Color &setRGB(const unsigned char _redValue,
                const unsigned char _greenValue,
                const unsigned char _blueValue);

  /**
   * Set the RGBA value to the current Color.
   *
   * @param[in]    _redValue     the red component. Value must be into [0,255]
   * interval.
   * @param[in]    _greenValue   the green component. Value must be into [0,255]
   * interval.
   * @param[in]    _blueValue    the blue component. Value must be into [0,255]
   * interval.
   * @param[in]    _alphaValue   the alpha transparency. Value must be into
   * [0,255] interval. (Default value = 255)
   */
  Color &setRGBA(const unsigned char _redValue,
                 const unsigned char _greenValue,
                 const unsigned char _blueValue,
                 const unsigned char _alphaValue = 255);

  /**
   * Set the red component.
   * @param[in]  _redValue     the red component.
   */
  void red(const unsigned char _redValue);

  /**
   * Set the green component.
   * @param[in]  _greenValue   the green component.
   */
  void green(const unsigned char _greenValue);

  /**
   * Set the blue component.
   * @param[in]  _blueValue    the blue component.
   */
  void blue(const unsigned char _blueValue);

  /**
   * Set the alpha component.
   * @param[in]  _alphaValue   the alpha component.
   */
  void alpha(const unsigned char _alphaValue);

  /**
   * Get the red component.
   * @return the red component.
   */
  unsigned char red() const;

  /**
   * Get the green component.
   * @return the green component.
   */
  unsigned char green() const;

  /**
   * Get the blue component.
   * @return the blue component.
   */
  unsigned char blue() const;

  /**
   * Get the alpha component.
   * @return the alpha component.
   */
  unsigned char alpha() const;

  ///////////////////////////////

  bool operator==(const Color &_color) const;

  bool operator!=(const Color &_color) const;

  bool operator<(const Color &_color) const;

  Color &operator=(const Color &_color);

  Color operator+(const Color &_color);

  Color &operator+=(const Color &_color);

  Color operator-(const Color &_color);

  Color &operator-=(const Color &_color);

  Color operator*(const double _scale);

  Color &operator*=(const double _scale);

  ///////////////////////////////

  void selfDisplay(std::ostream &out) const;

  ///////////////////////////////

  // Please refer to
  // https://www.google.com/design/spec/style/color.html#color-color-palette
  static Color None(void);
  static Color Black(void);
  static Color White(void);
  static Color Red(void);
  static Color Pink(void);
  static Color Purple(void);
  static Color DeepPurple(void);
  static Color Indigo(void);
  static Color Blue(void);
  static Color LightBlue(void);
  static Color Cyan(void);
  static Color Teal(void);
  static Color Green(void);
  static Color LightGreen(void);
  static Color Lime(void);
  static Color Yellow(void);
  static Color Amber(void);
  static Color Orange(void);
  static Color DeepOrange(void);
  static Color Brown(void);
  static Color Grey(void);
  static Color BlueGrey(void);

  ///////////////////////////////

  // Please refer to https://flatuicolors.com/
  static Color Turquoise(void);
  static Color GreenSea(void);
  static Color Emerald(void);
  static Color Nephritis(void);
  static Color PeterRiver(void);
  static Color BelizeHole(void);
  static Color Amethyst(void);
  static Color Wisteria(void);
  static Color WetAsphalt(void);
  static Color MidnightBlue(void);
  static Color SunFlower(void);
  static Color Carrot(void);
  static Color Pumpkin(void);
  static Color Alizarin(void);
  static Color Pomegranate(void);
  static Color Clouds(void);
  static Color Silver(void);
  static Color Concrete(void);
  static Color Asbestos(void);

private:
  unsigned char redComponent;
  unsigned char greenComponent;
  unsigned char blueComponent;
  unsigned char alphaComponent;
};

} // namespace FEVV


#ifndef Q_MOC_RUN
// implementation
#include "Color.inl"
#endif // Q_MOC_RUN
