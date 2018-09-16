#pragma once

#if defined(Color_RECURSES)
#error Recursive header files inclusion detected in Color.h
#else // defined(Color_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Color_RECURSES

#if !defined Color_h
/** Prevents repeated inclusion of headers. */
#define Color_h

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
  static const Color None;
  static const Color Black;
  static const Color White;
  static const Color Red;
  static const Color Pink;
  static const Color Purple;
  static const Color DeepPurple;
  static const Color Indigo;
  static const Color Blue;
  static const Color LightBlue;
  static const Color Cyan;
  static const Color Teal;
  static const Color Green;
  static const Color LightGreen;
  static const Color Lime;
  static const Color Yellow;
  static const Color Amber;
  static const Color Orange;
  static const Color DeepOrange;
  static const Color Brown;
  static const Color Grey;
  static const Color BlueGrey;

  ///////////////////////////////

  // Please refer to https://flatuicolors.com/
  static const Color Turquoise;
  static const Color GreenSea;
  static const Color Emerald;
  static const Color Nephritis;
  static const Color PeterRiver;
  static const Color BelizeHole;
  static const Color Amethyst;
  static const Color Wisteria;
  static const Color WetAsphalt;
  static const Color MidnightBlue;
  static const Color SunFlower;
  static const Color Carrot;
  static const Color Pumpkin;
  static const Color Alizarin;
  static const Color Pomegranate;
  static const Color Clouds;
  static const Color Silver;
  static const Color Concrete;
  static const Color Asbestos;

private:
  unsigned char redComponent;
  unsigned char greenComponent;
  unsigned char blueComponent;
  unsigned char alphaComponent;
};

} // namespace FEVV

#endif // !defined Color_h

#undef Color_RECURSES
#endif // else defined(Color_RECURSES)
