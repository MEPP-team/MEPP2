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

#include <vector>

namespace FEVV {
namespace Operators {
namespace Geometry {

/**
 * \brief  Predicate fonction to know whether the 
 *         Pos position is inside the abcd tetrahedron.
 *
 * \tparam  GeometryTraits The geometric kernel. 
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] d The fourth point.
 * \param[in] Pos The point that is either inside 
              or outside the tetrahedron (abcd). 
 * \param[in] gt The geometry trait object.			  
 * \return True if Pos is inside the tetrahedron (abcd),
                else false.
 */	
template< typename GeometryTraits >
static inline bool
is_in_tetrahedron(const typename GeometryTraits::Point &a, 
                  const typename GeometryTraits::Point &b, 
                  const typename GeometryTraits::Point &c, 
                  const typename GeometryTraits::Point &d, 
                  const typename GeometryTraits::Point &Pos,
				  const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename GeometryTraits::Scalar Scalar;
  
  Vector vap = Pos - a, vbp = Pos - b,
    // vcp = Pos - c,
    // vdp = Pos - d,
      vab = b - a, vac = c - a, vad = d - a, vbc = c - b, vbd = d - b;

  Vector temp_ac_ad(vac[1] * vad[2] - vac[2] * vad[1],
                    vac[2] * vad[0] - vac[0] * vad[2],
                    vac[0] * vad[1] - vac[1] * vad[0]);

  Vector temp(vbd[1] * vbc[2] - vbd[2] * vbc[1],
              vbd[2] * vbc[0] - vbd[0] * vbc[2],
              vbd[0] * vbc[1] - vbd[1] * vbc[0]);

  Scalar va = (vbp * temp) * 1.f / 6;

  Scalar vb = (vap * temp_ac_ad) * 1.f / 6;

  temp = Vector(vad[1] * vab[2] - vad[2] * vab[1],
           vad[2] * vab[0] - vad[0] * vab[2],
           vad[0] * vab[1] - vad[1] * vab[0]);

  Scalar vc = (vap * temp) * 1.f / 6;

  temp = Vector(vab[1] * vac[2] - vab[2] * vac[1],
                vab[2] * vac[0] - vab[0] * vac[2],
                vab[0] * vac[1] - vab[1] * vac[0]);

  Scalar vd = (vap * temp) * 1.f / 6;

  Scalar v = (vab * temp_ac_ad) * 1.f / 6;

  Scalar alpha = va / v, beta = vb / v, gamma = vc / v, delta = vd / v;
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // std::cout << "barycentric coordinates: (" << alpha << ", " << beta << ", "
  // << gamma << ", " << delta << ")" << std::endl ;

  return !((alpha < 0.f) || (beta < 0.f) || (gamma < 0.f) || (delta < 0.f));
}

/**
 * \brief  Trilinear interpolation of values at tetrahedron
 *         points.
 *
 * \tparam  GeometryTraits The geometric kernel. 
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] d The fourth point.
 * \param[in] Pos The point for which the trilinear 
              interpolation of values is computed.	
 * \param[in] vala The value of the first point.
 * \param[in] valb The value of the second point.
 * \param[in] valc The value of the third point.
 * \param[in] vald The value of the fourth point. 			  
 * \param[in] gt The geometry trait object.			  
 * \return The trilinear interpolation of values vala, valb, valc,
           and vald, at Pos position.
 * \pre    Pos must be inside the tetrahedron (abcd), else 
	       the trilinear interpolation coefficients are
	       truncated to [0,1].			  
 */
template< typename GeometryTraits >
static inline 
typename GeometryTraits::Scalar
trilinear_interpolation(const typename GeometryTraits::Point &a,
                        const typename GeometryTraits::Point &b,
                        const typename GeometryTraits::Point &c,
                        const typename GeometryTraits::Point &d,
                        const typename GeometryTraits::Point &Pos,
                        typename GeometryTraits::Scalar vala,
                        typename GeometryTraits::Scalar valb,
                        typename GeometryTraits::Scalar valc,
                        typename GeometryTraits::Scalar vald,
						const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename GeometryTraits::Scalar Scalar;
  
  Vector vap = Pos - a, vbp = Pos - b,
    // vcp = Pos - c,
    // vdp = Pos - d,
      vab = b - a, vac = c - a, vad = d - a, vbc = c - b, vbd = d - b;

  Vector temp_ac_ad(vac[1] * vad[2] - vac[2] * vad[1],
               vac[2] * vad[0] - vac[0] * vad[2],
               vac[0] * vad[1] - vac[1] * vad[0]);

  Vector temp(vbd[1] * vbc[2] - vbd[2] * vbc[1],
         vbd[2] * vbc[0] - vbd[0] * vbc[2],
         vbd[0] * vbc[1] - vbd[1] * vbc[0]);

  Scalar va = (vbp * temp) * 1.f / 6;

  Scalar vb = (vap * temp_ac_ad) * 1.f / 6;

  temp = Vector(vad[1] * vab[2] - vad[2] * vab[1],
                vad[2] * vab[0] - vad[0] * vab[2],
                vad[0] * vab[1] - vad[1] * vab[0]);

  Scalar vc = (vap * temp) * 1.f / 6;

  temp = Vector(vab[1] * vac[2] - vab[2] * vac[1],
                vab[2] * vac[0] - vab[0] * vac[2],
                vab[0] * vac[1] - vab[1] * vac[0]);

  Scalar vd = (vap * temp) * 1.f / 6;

  Scalar v = (vab * temp_ac_ad) * 1.f / 6;

  Scalar alpha = va / v, beta = vb / v, gamma = vc / v, delta = vd / v;
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if(alpha < 0.f)
    alpha = 0.f;
  if(alpha > 1.f)
    alpha = 1.f;

  if(beta < 0.f)
    beta = 0.f;
  if(beta > 1.f)
    beta = 1.f;

  if(gamma < 0.f)
    gamma = 0.f;
  if(gamma > 1.f)
    gamma = 1.f;

  if(delta < 0.f)
    delta = 0.f;
  if(delta > 1.f)
    delta = 1.f;
  /////////////////////////////////////////////////////////////////////////////////////////////////
  return alpha * vala + beta * valb + gamma * valc + delta * vald;
};

/**
 * \brief  Trilinear interpolation of values at tetrahedron
 *         points.
 *
 * \tparam  GeometryTraits The geometric kernel. 
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] d The fourth point.
 * \param[in] Pos The point for which the trilinear 
              interpolations of values are computed.	
 * \param[in] vala The vector of values at the first point.
 * \param[in] valb The vector of values at the second point.
 * \param[in] valc The vector of values at the third point.
 * \param[in] vald The vector of values at the fourth point. 			  
 * \param[in] gt The geometry trait object.			  
 * \return The vector of the trilinear interpolation of values 
           at a, b, c, d positions for Pos position.
 * \pre    Pos must be inside the tetrahedron (abcd), else 
	       the trilinear interpolation coefficients are
	       truncated to [0,1].			  
 */
template< typename GeometryTraits >
static inline std::vector< typename GeometryTraits::Scalar >
trilinear_interpolation(const typename GeometryTraits::Point &a,
                        const typename GeometryTraits::Point &b,
                        const typename GeometryTraits::Point &c,
                        const typename GeometryTraits::Point &d,
                        const typename GeometryTraits::Point &Pos,
                        const std::vector< typename GeometryTraits::Scalar > &vala,
                        const std::vector< typename GeometryTraits::Scalar > &valb,
                        const std::vector< typename GeometryTraits::Scalar > &valc,
                        const std::vector< typename GeometryTraits::Scalar > &vald,
						const GeometryTraits &gt)
{
  size_t nbe = vala.size();
  if(valb.size() < nbe)
    nbe = valb.size();
  if(valc.size() < nbe)
    nbe = valc.size();
  if(vald.size() < nbe)
    nbe = vald.size();
  ///////////////////////////////////////////////////////////////////
  std::vector< typename GeometryTraits::Scalar > res(nbe, 0.f);

  for(size_t i = 0; i < nbe; ++i)
  {
    res[i] = trilinear_interpolation< GeometryTraits >(
        a, b, c, d, Pos, vala[i], valb[i], valc[i], vald[i], gt);
  }

  return res;
};

/**
 * \brief  Predicate fonction to know whether the 
 *         Pos position is over the ab segment.
 *
 * \tparam  GeometryTraits The geometric kernel.
 * \param[in]  a The first point.
 * \param[in]  b The second point.
 * \param[in]  Pos The point that is either over 
               or beside the segment (ab). 
 * \param[out] diff	The orthogonal vector to translate
 *             Pos onto the segment (ab) (or line when
 *             Pos is not over the segment). 
 * \param[in]  gt The geometry trait object.	 
 * \return True if Pos is over the segment (ab),
                else false.
 */	
template< typename GeometryTraits >
static inline bool
is_over_segment(const typename GeometryTraits::Point &a,
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

/**
 * \brief  Predicate fonction to know whether the 
 *         projection of point Pos on the triangle 
 *         plane is in the triangle (abc).
 *
 * \tparam  GeometryTraits The geometric kernel. 
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] Pos The point whose projection is either  
              inside or outside the triangle (abc).
 * \param[in] gt The geometry trait object.			  
 * \return True if the projection of point Pos is inside
 *         the triangle (abc), else false.
 */
template< typename GeometryTraits >
static inline bool
is_in_triangle(const typename GeometryTraits::Point &a,
               const typename GeometryTraits::Point &b,
               const typename GeometryTraits::Point &c,
               const typename GeometryTraits::Point &Pos,
               const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename GeometryTraits::Scalar Scalar;
  // Pos =
  // we are going to compute the barycentric coordinates of Pos
  Vector p1p(a[0] - c[0], a[1] - c[1], a[2] - c[2]),
      p2p(b[0] - c[0], b[1] - c[1], b[2] - c[2]),
      Posp(Pos[0] - c[0], Pos[1] - c[1], Pos[2] - c[2]);

  Scalar p1p2 = gt.dot_product(p1p, p1p), p2p2 = gt.dot_product(p2p, p2p),
         p1pp2p = gt.dot_product(p1p, p2p), p1pPosp = gt.dot_product(p1p, Posp),
         p2pPosp = gt.dot_product(p2p, Posp);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  Scalar den = p1p2 * p2p2 - p1pp2p * p1pp2p;
  Scalar alpha = (p2p2 * p1pPosp - p1pp2p * p2pPosp) / den,
         beta = (p1p2 * p2pPosp - p1pp2p * p1pPosp) / den, gamma;
  gamma = 1 - (alpha + beta);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  return !((alpha < 0.0) || (beta < 0.0) || (gamma < 0.0) || (alpha > 1.0) ||
           (beta > 1.0) || (gamma > 1.0));
};

/**
 * \brief  Bilinear interpolation of values at triangle
 *         points.
 *
 * \tparam  GeometryTraits The geometric kernel. 
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] Pos The point for which the bilinear 
              interpolation of values is computed.	
 * \param[in] vala The value of the first point.
 * \param[in] valb The value of the second point.
 * \param[in] valc The value of the third point.		  
 * \param[in] gt The geometry trait object.			  
 * \return The bilinear interpolation of values vala, valb, valc,
              and vald, at Pos position.
 * \pre    Pos must be inside the triangle (abc), else 
	       the bilinear interpolation coefficients are
	       truncated to [0,1].			  
 */
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
  Scalar p1p2 = gt.dot_product(p1p, p1p), p2p2 = gt.dot_product(p2p, p2p),
         p1pp2p = gt.dot_product(p1p, p2p), p1pPosp = gt.dot_product(p1p, Posp),
         p2pPosp = gt.dot_product(p2p, Posp);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  Scalar den = p1p2 * p2p2 - p1pp2p * p1pp2p;
  Scalar alpha = (p2p2 * p1pPosp - p1pp2p * p2pPosp) / den,
         beta = (p1p2 * p2pPosp - p1pp2p * p1pPosp) / den, gamma;
  gamma = 1 - (alpha + beta);
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if(alpha < 0.f)
    alpha = 0.f;
  if(alpha > 1.f)
    alpha = 1.f;

  if(beta < 0.f)
    beta = 0.f;
  if(beta > 1.f)
    beta = 1.f;

  if(gamma < 0.f)
    gamma = 0.f;
  if(gamma > 1.f)
    gamma = 1.f;
  /////////////////////////////////////////////////////////////////////////////////////////////////
  return alpha * vala + beta * valb + gamma * valc;
};
/**
 * \brief  Bilinear interpolation of values at triangle
 *         points.
 *
 * \tparam  GeometryTraits The geometric kernel. 
 * \param[in] a The first point.
 * \param[in] b The second point.
 * \param[in] c The third point.
 * \param[in] Pos The point for which the bilinear 
              interpolations of values are computed.	
 * \param[in] vala The vector of values at the first point.
 * \param[in] valb The vector of values at the second point.
 * \param[in] valc The vector of values at the third point.			  
 * \param[in] gt The geometry trait object.			  
 * \return The vector of the bilinear interpolation of values 
           at a, b, c positions for Pos position.
 * \pre    Pos must be inside the triangle (abc), else 
	       the bilinear interpolation coefficients are
	       truncated to [0,1].			  
 */
template< typename GeometryTraits >
static inline std::vector< typename GeometryTraits::Scalar >
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
} // namespace Geometry
} // namespace Operators
} // namespace FEVV
