#ifndef CLIPPING_AND_INTERSECTIONS_HPP
#define CLIPPING_AND_INTERSECTIONS_HPP

#include "FEVV/Tools/Math/MatrixOperations.hpp"

namespace FEVV {
namespace Operators {

namespace Geometry {

/**
 * Sphere equation is given by (x-Cx)^2+(y-Cy)^2+(z-Cz)^2=r^2
 * Line equation is given by (x, y, z) = P + txV
 */
template< typename GeometryTraits >
static bool
sphere_clip_vector(
    const typename GeometryTraits::Point &center, /// Sphere center
    double r,                                     /// Sphere radius
    const typename GeometryTraits::Point
        &p, /// Starting point/origin of the line (we used it currently for P
            /// inside the sphere)
    typename GeometryTraits::Vector
        &v, /// vector/ray/line direction to clip by the sphere when needed
            /// [starting from P]
    const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  if(r < 0)
    r = -r; // the radius cannot be negative

  Vector w = gt.sub(p, center); // direction towards P from the sphere center
  double a = gt.dot_product(v, v);

  if(fabs(a) < std::numeric_limits< double >::epsilon())
    return false;

  double b = 2.0 * gt.dot_product(v, w);
  double c = gt.dot_product(w, w) - r * r;
  double delta = b * b - 4. * a * c;
  if(delta < 0.)
  {
    // Should not happen, but happens sometimes (numerical precision)
    return true;
  }

  double t1 = (-b + ::sqrt(delta)) /
              (2.0 * a) /*, t2 = (-b - ::sqrt(delta)) / (2.0 * a)*/;

  if(t1 >= 1.)
  {
    // Inside the sphere
    return false;
  }

  // if (t1 < 0.) {
  // Should not happen, but happens sometimes (numerical precision or P not
  // inside the sphere)
  //	return true;
  //}

  v = gt.scalar_mult(v, t1);

  return true;
}
} // namespace Geometry

} // namespace Operators
} // namespace FEVV

#endif // CLIPPING_AND_INTERSECTIONS_HPP
