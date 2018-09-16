
#ifndef __AngleOperations_hxx
#define __AngleOperations_hxx

#include <cmath> // M_PI, asin

#define DEG2RAD(deg) (deg * (M_PI / 180.f))
#define RAD2DEG(rad) (rad * (180.f / M_PI))

namespace FEVV {
namespace Math {

namespace Angle {

/// sine can be either a float, a double or a long double
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

/// sine can be either a float, a double or a long double
template< class T >
static T
asin_degree(T sine)
{
  return RAD2DEG(Angle::asin< T >(sine));
}

///////////////////////////////////////////////////////////////////////////
/// cosine can be either a float, a double or a long double
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
/// cosine can be either a float, a double or a long double
template< class T >
inline T
acos_degree(T cosine)
{
  return RAD2DEG(Angle::acos< T >(cosine));
}
} // namespace Angle

} // namespace Math
} // namespace FEVV

#endif /* __AngleOperations_hxx */
