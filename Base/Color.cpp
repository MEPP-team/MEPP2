#include "Base/Color.h"

#include <boost/algorithm/clamp.hpp>
#include <algorithm>

const FEVV::Color FEVV::Color::None(0, 0, 0, 0);
const FEVV::Color
    FEVV::Color::Black((unsigned char)0, (unsigned char)0, (unsigned char)0);
const FEVV::Color FEVV::Color::White((unsigned char)255,
                                     (unsigned char)255,
                                     (unsigned char)255);
const FEVV::Color
    FEVV::Color::Red((unsigned char)244, (unsigned char)67, (unsigned char)54);
const FEVV::Color
    FEVV::Color::Pink((unsigned char)233, (unsigned char)30, (unsigned char)99);
const FEVV::Color FEVV::Color::Purple((unsigned char)156,
                                      (unsigned char)39,
                                      (unsigned char)176);
const FEVV::Color FEVV::Color::DeepPurple((unsigned char)103,
                                          (unsigned char)58,
                                          (unsigned char)183);
const FEVV::Color FEVV::Color::Indigo((unsigned char)63,
                                      (unsigned char)81,
                                      (unsigned char)181);
const FEVV::Color FEVV::Color::Blue((unsigned char)33,
                                    (unsigned char)150,
                                    (unsigned char)243);
const FEVV::Color FEVV::Color::LightBlue((unsigned char)3,
                                         (unsigned char)169,
                                         (unsigned char)244);
const FEVV::Color
    FEVV::Color::Cyan((unsigned char)0, (unsigned char)188, (unsigned char)212);
const FEVV::Color
    FEVV::Color::Teal((unsigned char)0, (unsigned char)150, (unsigned char)136);
const FEVV::Color FEVV::Color::Green((unsigned char)76,
                                     (unsigned char)175,
                                     (unsigned char)80);
const FEVV::Color FEVV::Color::LightGreen((unsigned char)139,
                                          (unsigned char)195,
                                          (unsigned char)74);
const FEVV::Color FEVV::Color::Lime((unsigned char)205,
                                    (unsigned char)220,
                                    (unsigned char)57);
const FEVV::Color FEVV::Color::Yellow((unsigned char)255,
                                      (unsigned char)235,
                                      (unsigned char)59);
const FEVV::Color FEVV::Color::Amber((unsigned char)255,
                                     (unsigned char)193,
                                     (unsigned char)7);
const FEVV::Color FEVV::Color::Orange((unsigned char)243,
                                      (unsigned char)156,
                                      (unsigned char)18);
const FEVV::Color FEVV::Color::DeepOrange((unsigned char)255,
                                          (unsigned char)87,
                                          (unsigned char)34);
const FEVV::Color FEVV::Color::Brown((unsigned char)121,
                                     (unsigned char)85,
                                     (unsigned char)72);
const FEVV::Color FEVV::Color::Grey((unsigned char)158,
                                    (unsigned char)158,
                                    (unsigned char)158);
const FEVV::Color FEVV::Color::BlueGrey((unsigned char)96,
                                        (unsigned char)125,
                                        (unsigned char)139);

///////////////////////////////

const FEVV::Color FEVV::Color::Turquoise((unsigned char)26,
                                         (unsigned char)188,
                                         (unsigned char)156);
const FEVV::Color FEVV::Color::GreenSea((unsigned char)22,
                                        (unsigned char)160,
                                        (unsigned char)133);
const FEVV::Color FEVV::Color::Emerald((unsigned char)46,
                                       (unsigned char)204,
                                       (unsigned char)113);
const FEVV::Color FEVV::Color::Nephritis((unsigned char)39,
                                         (unsigned char)174,
                                         (unsigned char)96);
const FEVV::Color FEVV::Color::PeterRiver((unsigned char)41,
                                          (unsigned char)128,
                                          (unsigned char)185);
const FEVV::Color FEVV::Color::BelizeHole((unsigned char)41,
                                          (unsigned char)128,
                                          (unsigned char)185);
const FEVV::Color FEVV::Color::Amethyst((unsigned char)155,
                                        (unsigned char)89,
                                        (unsigned char)182);
const FEVV::Color FEVV::Color::Wisteria((unsigned char)142,
                                        (unsigned char)68,
                                        (unsigned char)173);
const FEVV::Color FEVV::Color::WetAsphalt((unsigned char)52,
                                          (unsigned char)73,
                                          (unsigned char)94);
const FEVV::Color FEVV::Color::MidnightBlue((unsigned char)44,
                                            (unsigned char)62,
                                            (unsigned char)80);
const FEVV::Color FEVV::Color::SunFlower((unsigned char)241,
                                         (unsigned char)196,
                                         (unsigned char)15);
const FEVV::Color FEVV::Color::Carrot((unsigned char)230,
                                      (unsigned char)126,
                                      (unsigned char)34);
const FEVV::Color FEVV::Color::Pumpkin((unsigned char)211,
                                       (unsigned char)84,
                                       (unsigned char)0);
const FEVV::Color FEVV::Color::Alizarin((unsigned char)231,
                                        (unsigned char)76,
                                        (unsigned char)60);
const FEVV::Color FEVV::Color::Pomegranate((unsigned char)192,
                                           (unsigned char)57,
                                           (unsigned char)43);
const FEVV::Color FEVV::Color::Clouds((unsigned char)236,
                                      (unsigned char)240,
                                      (unsigned char)241);
const FEVV::Color FEVV::Color::Silver((unsigned char)189,
                                      (unsigned char)195,
                                      (unsigned char)199);
const FEVV::Color FEVV::Color::Concrete((unsigned char)149,
                                        (unsigned char)165,
                                        (unsigned char)166);
const FEVV::Color FEVV::Color::Asbestos((unsigned char)127,
                                        (unsigned char)140,
                                        (unsigned char)141);

///////////////////////////////

FEVV::Color::Color()
    : redComponent((unsigned char)0), greenComponent((unsigned char)0),
      blueComponent((unsigned char)0), alphaComponent((unsigned char)255)
{
}

FEVV::Color::Color(const unsigned char _redValue,
                   const unsigned char _greenValue,
                   const unsigned char _blueValue,
                   const unsigned char _alphaValue)
    : redComponent(_redValue), greenComponent(_greenValue),
      blueComponent(_blueValue), alphaComponent(_alphaValue)
{
}

FEVV::Color::Color(const unsigned char _grayValue,
                   const unsigned char _alphaValue)
    : redComponent(_grayValue), greenComponent(_grayValue),
      blueComponent(_grayValue), alphaComponent(_alphaValue)
{
}

FEVV::Color::Color(const Color &_color)
    : redComponent(_color.redComponent), greenComponent(_color.greenComponent),
      blueComponent(_color.blueComponent), alphaComponent(_color.alphaComponent)
{
}

///////////////////////////////

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

void
FEVV::Color::red(const unsigned char _redValue)
{
  redComponent = _redValue;
}

void
FEVV::Color::green(const unsigned char _greenValue)
{
  greenComponent = _greenValue;
}

void
FEVV::Color::blue(const unsigned char _blueValue)
{
  blueComponent = _blueValue;
}

void
FEVV::Color::alpha(const unsigned char _alphaValue)
{
  alphaComponent = _alphaValue;
}

unsigned char
FEVV::Color::red() const
{
  return redComponent;
}

unsigned char
FEVV::Color::green() const
{
  return greenComponent;
}

unsigned char
FEVV::Color::blue() const
{
  return blueComponent;
}

unsigned char
FEVV::Color::alpha() const
{
  return alphaComponent;
}

///////////////////////////////

bool
FEVV::Color::operator==(const Color &_color) const
{
  return redComponent == _color.redComponent &&
         greenComponent == _color.greenComponent &&
         blueComponent == _color.blueComponent &&
         alphaComponent == _color.alphaComponent;
}

bool
FEVV::Color::operator!=(const Color &_color) const
{
  return redComponent != _color.redComponent ||
         greenComponent != _color.greenComponent ||
         blueComponent != _color.blueComponent ||
         alphaComponent != _color.alphaComponent;
}

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

FEVV::Color &
FEVV::Color::operator=(const Color &_color)
{
  redComponent = _color.redComponent;
  greenComponent = _color.greenComponent;
  blueComponent = _color.blueComponent;
  alphaComponent = _color.alphaComponent;

  return *this;
}

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

void
FEVV::Color::selfDisplay(std::ostream &_out) const
{
  _out << "[Color] RGBA(" << (unsigned int)redComponent << ","
       << (unsigned int)greenComponent << "," << (unsigned int)blueComponent
       << "," << (unsigned int)alphaComponent << ")";
}
