/**
 * \file		Interpolation.hxx
 * \author		Vincent Vidal
 * \date     	2015-04-14
 */

#ifndef __Interpolation_hxx
#define __Interpolation_hxx
#include <vector>

namespace FEVV {
namespace Math {

namespace Barycentric {
namespace Tetrahedron {
template< class P, class V >
static inline bool
isInTetrahedron(const P &a, const P &b, const P &c, const P &d, const P &Pos)
{
  V vap = Pos - a, vbp = Pos - b,
    // vcp = Pos - c,
    // vdp = Pos - d,
      vab = b - a, vac = c - a, vad = d - a, vbc = c - b, vbd = d - b;

  V temp_ac_ad(vac[1] * vad[2] - vac[2] * vad[1],
               vac[2] * vad[0] - vac[0] * vad[2],
               vac[0] * vad[1] - vac[1] * vad[0]);

  V temp(vbd[1] * vbc[2] - vbd[2] * vbc[1],
         vbd[2] * vbc[0] - vbd[0] * vbc[2],
         vbd[0] * vbc[1] - vbd[1] * vbc[0]);
#ifdef _KERNEL_EXACT_
  double va = CGAL::to_double(vbp * temp) * 1. / 6;
#else
  double va = (vbp * temp) * 1. / 6;
#endif

#ifdef _KERNEL_EXACT_
  double vb = CGAL::to_double(vap * temp_ac_ad) * 1. / 6;
#else
  double vb = (vap * temp_ac_ad) * 1. / 6;
#endif

  temp = V(vad[1] * vab[2] - vad[2] * vab[1],
           vad[2] * vab[0] - vad[0] * vab[2],
           vad[0] * vab[1] - vad[1] * vab[0]);
#ifdef _KERNEL_EXACT_
  double vc = CGAL::to_double(vap * temp) * 1. / 6;
#else
  double vc = (vap * temp) * 1. / 6;
#endif

  temp = V(vab[1] * vac[2] - vab[2] * vac[1],
           vab[2] * vac[0] - vab[0] * vac[2],
           vab[0] * vac[1] - vab[1] * vac[0]);
#ifdef _KERNEL_EXACT_
  double vd = CGAL::to_double(vap * temp) * 1. / 6;
#else
  double vd = (vap * temp) * 1. / 6;
#endif

#ifdef _KERNEL_EXACT_
  double v = CGAL::to_double(vab * temp_ac_ad) * 1. / 6;
#else
  double v = (vab * temp_ac_ad) * 1. / 6;
#endif
  double alpha = va / v, beta = vb / v, gamma = vc / v, delta = vd / v;
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // std::cout << "barycentric coordinates: (" << alpha << ", " << beta << ", "
  // << gamma << ", " << delta << ")" << std::endl ;

  return !((alpha < 0.0) || (beta < 0.0) || (gamma < 0.0) || (delta < 0.0));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class P, class V >
static inline double
trilinear_interpolation(const P &a,
                        const P &b,
                        const P &c,
                        const P &d,
                        const P &Pos,
                        double vala,
                        double valb,
                        double valc,
                        double vald)
{
  V vap = Pos - a, vbp = Pos - b,
    // vcp = Pos - c,
    // vdp = Pos - d,
      vab = b - a, vac = c - a, vad = d - a, vbc = c - b, vbd = d - b;

  V temp_ac_ad(vac[1] * vad[2] - vac[2] * vad[1],
               vac[2] * vad[0] - vac[0] * vad[2],
               vac[0] * vad[1] - vac[1] * vad[0]);

  V temp(vbd[1] * vbc[2] - vbd[2] * vbc[1],
         vbd[2] * vbc[0] - vbd[0] * vbc[2],
         vbd[0] * vbc[1] - vbd[1] * vbc[0]);
#ifdef _KERNEL_EXACT_
  double va = CGAL::to_double(vbp * temp) * 1. / 6;
#else
  double va = (vbp * temp) * 1. / 6;
#endif

#ifdef _KERNEL_EXACT_
  double vb = CGAL::to_double(vap * temp_ac_ad) * 1. / 6;
#else
  double vb = (vap * temp_ac_ad) * 1. / 6;
#endif

  temp = V(vad[1] * vab[2] - vad[2] * vab[1],
           vad[2] * vab[0] - vad[0] * vab[2],
           vad[0] * vab[1] - vad[1] * vab[0]);
#ifdef _KERNEL_EXACT_
  double vc = CGAL::to_double(vap * temp) * 1. / 6;
#else
  double vc = (vap * temp) * 1. / 6;
#endif

  temp = V(vab[1] * vac[2] - vab[2] * vac[1],
           vab[2] * vac[0] - vab[0] * vac[2],
           vab[0] * vac[1] - vab[1] * vac[0]);
#ifdef _KERNEL_EXACT_
  double vd = CGAL::to_double(vap * temp) * 1. / 6;
#else
  double vd = (vap * temp) * 1. / 6;
#endif

#ifdef _KERNEL_EXACT_
  double v = CGAL::to_double(vab * temp_ac_ad) * 1. / 6;
#else
  double v = (vab * temp_ac_ad) * 1. / 6;
#endif
  double alpha = va / v, beta = vb / v, gamma = vc / v, delta = vd / v;
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if(alpha < 0)
    alpha = 0.0;
  if(alpha > 1)
    alpha = 1.0;

  if(beta < 0)
    beta = 0.0;
  if(beta > 1)
    beta = 1.0;

  if(gamma < 0)
    gamma = 0.0;
  if(gamma > 1)
    gamma = 1.0;

  if(delta < 0)
    delta = 0.0;
  if(delta > 1)
    delta = 1.0;
  /////////////////////////////////////////////////////////////////////////////////////////////////
  return alpha * vala + beta * valb + gamma * valc + delta * vald;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class P, class V >
static inline std::vector< double >
trilinear_interpolation(const P &a,
                        const P &b,
                        const P &c,
                        const P &d,
                        const P &Pos,
                        const std::vector< double > &vala,
                        const std::vector< double > &valb,
                        const std::vector< double > &valc,
                        const std::vector< double > &vald)
{
  size_t nbe = vala.size();
  if(valb.size() < nbe)
    nbe = valb.size();
  if(valc.size() < nbe)
    nbe = valc.size();
  if(vald.size() < nbe)
    nbe = vald.size();
  ///////////////////////////////////////////////////////////////////
  std::vector< double > res(nbe, 0.);

  for(size_t i = 0; i < nbe; ++i)
  {
    res[i] = trilinear_interpolation< P, V >(
        a, b, c, d, Pos, vala[i], valb[i], valc[i], vald[i]);
  }

  return res;
};
} // namespace Tetrahedron
namespace Triangle {
//------------------------------------------------------------------------------------------------
template< typename GeometryTraits >
static inline bool
isOverSegment(const typename GeometryTraits::Point &a,
              const typename GeometryTraits::Point &b,
              const typename GeometryTraits::Point &Pos,
              typename GeometryTraits::Vector &diff,
              const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename GeometryTraits::Point Point;
  Vector vab(gt.sub(b, a)), vaPos(gt.sub(Pos, a)), inter, N, Ntmp;
  Point interp;
  inter = gt.cross_product(vab, vaPos);
  // std::cout << " inter: (" << inter[0] << " " << inter[1] << " " << inter[2]
  // << ")" << std::endl;
  Ntmp = gt.cross_product(inter, vab);
  // std::cout << " Ntmp: (" << Ntmp[0] << " " << Ntmp[1] << " " << Ntmp[2] <<
  // ")" << std::endl;
  N = gt.normalize(Ntmp);
  // std::cout << " N: (" << N[0] << " " << N[1] << " " << N[2] << ")" <<
  // std::endl;
  diff = gt.scalar_mult(-gt.dot_product(vaPos, N),
                        N); // "-" to go towards the surface!!

  interp = Point(Pos[0] + diff[0], Pos[1] + diff[1], Pos[2] + diff[2]);
  // std::cout << " interp: (" << interp[0] << " " << interp[1] << " " <<
  // interp[2] << ")" << std::endl;
  return (gt.dot_product(
              Vector(interp[0] - a[0], interp[1] - a[1], interp[2] - a[2]),
              vab) >= 0) &&
         (gt.dot_product(
              Vector(b[0] - interp[0], b[1] - interp[1], b[2] - interp[2]),
              vab) >= 0);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
// is the projection of point Pos on the triangle plane in the triangle (a,b,c)?
template< typename GeometryTraits >
static inline bool
isInTriangle(const typename GeometryTraits::Point &a,
             const typename GeometryTraits::Point &b,
             const typename GeometryTraits::Point &c,
             const typename GeometryTraits::Point &Pos,
             const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename GeometryTraits::Scalar Scalar;
  // Pos =
  // we are going to compute the barycentric coordinates of Pos   (from
  // numerical recipes)
  Vector p1p(a[0] - c[0], a[1] - c[1], a[2] - c[2]),
      p2p(b[0] - c[0], b[1] - c[1], b[2] - c[2]),
      Posp(Pos[0] - c[0], Pos[1] - c[1], Pos[2] - c[2]);
#ifdef _KERNEL_EXACT_
  double p1p2 = CGAL::to_double(p1p * p1p), p2p2 = CGAL::to_double(p2p * p2p),
         p1pp2p = CGAL::to_double(p1p * p2p),
         p1pPosp = CGAL::to_double(p1p * Posp),
         p2pPosp = CGAL::to_double(p2p * Posp);
#else
  Scalar p1p2 = gt.dot_product(p1p, p1p), p2p2 = gt.dot_product(p2p, p2p),
         p1pp2p = gt.dot_product(p1p, p2p), p1pPosp = gt.dot_product(p1p, Posp),
         p2pPosp = gt.dot_product(p2p, Posp);
#endif
  /////////////////////////////////////////////////////////////////////////////////////////////////
  Scalar den = p1p2 * p2p2 - p1pp2p * p1pp2p;
  Scalar alpha = (p2p2 * p1pPosp - p1pp2p * p2pPosp) / den,
         beta = (p1p2 * p2pPosp - p1pp2p * p1pPosp) / den, gamma;
  gamma = 1 - (alpha + beta);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  return !((alpha < 0.0) || (beta < 0.0) || (gamma < 0.0) || (alpha > 1.0) ||
           (beta > 1.0) || (gamma > 1.0));
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename GeometryTraits >
static inline typename GeometryTraits::Scalar
bilinear_interpolation(const typename GeometryTraits::Point &a,
                       const typename GeometryTraits::Point &b,
                       const typename GeometryTraits::Point &c,
                       const typename GeometryTraits::Point &Pos,
                       typename GeometryTraits::Scalar vala,
                       typename GeometryTraits::Scalar valb,
                       typename GeometryTraits::Scalar valc,
                       const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename GeometryTraits::Scalar Scalar;

  Vector p1p = a - c, p2p = b - c, Posp = Pos - c;
#ifdef _KERNEL_EXACT_
  double p1p2 = CGAL::to_double(gt.dot_product(p1p, p1p)),
         p2p2 = CGAL::to_double(gt.dot_product(p2p, p2p)),
         p1pp2p = CGAL::to_double(gt.dot_product(p1p, p2p)),
         p1pPosp = CGAL::to_double(gt.dot_product(p1p, Posp)),
         p2pPosp = CGAL::to_double(gt.dot_product(p2p, Posp));
#else
  Scalar p1p2 = gt.dot_product(p1p, p1p), p2p2 = gt.dot_product(p2p, p2p),
         p1pp2p = gt.dot_product(p1p, p2p), p1pPosp = gt.dot_product(p1p, Posp),
         p2pPosp = gt.dot_product(p2p, Posp);
#endif
  /////////////////////////////////////////////////////////////////////////////////////////////////
  Scalar den = p1p2 * p2p2 - p1pp2p * p1pp2p;
  Scalar alpha = (p2p2 * p1pPosp - p1pp2p * p2pPosp) / den,
         beta = (p1p2 * p2pPosp - p1pp2p * p1pPosp) / den, gamma;
  gamma = 1 - (alpha + beta);
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if(alpha < 0)
    alpha = 0;
  if(alpha > 1)
    alpha = 1;

  if(beta < 0)
    beta = 0;
  if(beta > 1)
    beta = 1;

  if(gamma < 0)
    gamma = 0;
  if(gamma > 1)
    gamma = 1;
  /////////////////////////////////////////////////////////////////////////////////////////////////
  return alpha * vala + beta * valb + gamma * valc;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename GeometryTraits >
static inline std::vector< double >
bilinear_interpolation(
    const typename GeometryTraits::Point &a,
    const typename GeometryTraits::Point &b,
    const typename GeometryTraits::Point &c,
    const typename GeometryTraits::Point &Pos,
    const std::vector< typename GeometryTraits::Scalar > &vala,
    const std::vector< typename GeometryTraits::Scalar > &valb,
    const std::vector< typename GeometryTraits::Scalar > &valc,
    const GeometryTraits &gt)
{
  size_t nbe = vala.size();
  if(valb.size() < nbe)
    nbe = valb.size();
  if(valc.size() < nbe)
    nbe = valc.size();
  ///////////////////////////////////////////////////////////////////
  std::vector< typename GeometryTraits::Scalar > res(nbe, 0);

  for(size_t i = 0; i < nbe; ++i)
  {
    res[i] = bilinear_interpolation< GeometryTraits >(
        a, b, c, Pos, vala[i], valb[i], valc[i], gt);
  }

  return res;
};
} // namespace Triangle
} // namespace Barycentric

} // namespace Math
} // namespace FEVV
#endif
