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

#include "degree_rad_conversion.h" // rad2deg()

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm> // std::sort
#include <cassert>
#include <cmath> // floor

#ifndef ABS_GUY
#define ABS_GUY fabs
#endif
#ifndef SQRT_GUY
#define SQRT_GUY sqrt
#endif
#ifndef MACH_EPS_DOUBLE
#define MACH_EPS_DOUBLE 2e-16
#endif
#ifndef MACH_EPS_FLOAT
#define MACH_EPS_FLOAT 2e-7
#endif

namespace FEVV {
namespace Math {

namespace Vector {
namespace Stats {
///////////////////////////////////////////////////////////////////////////
template< class ElementType >
static std::vector< ElementType >
unique(const std::vector< ElementType > &v)
{
  std::vector< ElementType > ret;
  size_t nb_e = v.size();
  if(nb_e == 0)
    return ret;

  ret.push_back(v[0]);
  for(size_t i = 1; i < nb_e; ++i)
  {
    bool is_unique = true;
    for(size_t j = 0; j < ret.size(); j++)
    {
      if(fabs(v[i] - ret[j]) < MACH_EPS_FLOAT)
      {
        is_unique = false;
        break;
      }
    }
    if(is_unique)
    {
      ret.push_back(v[i]);
    }
  }

  return ret;
}
///////////////////////////////////////////////////////////////////////////
template< typename ElementType, int DIM >
static ElementType
maximum(const ElementType v[DIM])
{
  ElementType result = std::numeric_limits< ElementType >::min();
  for(size_t i = 0; i < DIM; ++i)
  {
    if(v[i] > result)
    {
      result = v[i];
    }
  }
  return result;
}
template< typename ElementType >
static ElementType
maximum(const std::vector< ElementType > &v)
{
  ElementType result = std::numeric_limits< ElementType >::min();
  const size_t dim = v.size();
  for(size_t i = 0; i < dim; ++i)
  {
    if(v[i] > result)
    {
      result = v[i];
    }
  }
  return result;
}
///////////////////////////////////////////////////////////////////////////
template< typename ElementType, int DIM >
static ElementType
minimum(const ElementType v[DIM])
{
  ElementType result = std::numeric_limits< ElementType >::max();
  for(size_t i = 0; i < DIM; ++i)
  {
    if(v[i] < result)
    {
      result = v[i];
    }
  }
  return result;
}
template< typename ElementType >
static ElementType
minimum(const std::vector< ElementType > &v)
{
  ElementType result = std::numeric_limits< ElementType >::max();
  const size_t dim = v.size();
  for(size_t i = 0; i < dim; ++i)
  {
    if(v[i] < result)
    {
      result = v[i];
    }
  }
  return result;
}
///////////////////////////////////////////////////////////////////////////
template< typename ElementType, int DIM >
static ElementType
mean(const ElementType v[DIM])
{
  ElementType result = 0;
  for(size_t i = 0; i < DIM; ++i)
  {
    result += v[i];
  }
  if(DIM > 0)
    return result / static_cast< ElementType >(DIM);
  else
    return result;
}
template< typename ElementType >
static ElementType
mean(const std::vector< ElementType > &v)
{
  ElementType result = 0;
  const size_t dim = v.size();
  for(size_t i = 0; i < dim; ++i)
  {
    result += v[i];
  }
  if(dim > 0)
    return result / static_cast< ElementType >(dim);
  else
    return result;
}
///////////////////////////////////////////////////////////////////////////
template< typename ElementType >
static ElementType
mean2(const std::vector< ElementType > &v)
{
  ElementType result = 0;
  const size_t dim = v.size();
  for(size_t i = 0; i < dim; ++i)
  {
    result += v[i] * v[i];
  }
  if(dim > 0)
    return result / static_cast< ElementType >(dim);
  else
    return result;
}

template< typename ElementType >
static ElementType
mean4(const std::vector< ElementType > &v)
{
  ElementType result = 0;
  const size_t dim = v.size();
  for(size_t i = 0; i < dim; ++i)
  {
    result += v[i] * v[i] * v[i] * v[i];
  }
  if(dim > 0)
    return result / static_cast< ElementType >(dim);
  else
    return result;
}
/// Be careful, this function will not work if any of the values is negative.
template< typename ElementType >
static ElementType
mean_sqrt(const std::vector< ElementType > &v)
{
  ElementType result = 0;
  const size_t dim = v.size();
  for(size_t i = 0; i < dim; ++i)
  {
    result += sqrt(v[i]);
  }
  if(dim > 0)
    return result / static_cast< ElementType >(dim);
  else
    return result;
}
/// Be careful, this function will not work if any of the values is negative.
template< typename ElementType >
static ElementType
mean_sqrt_sqrt(const std::vector< ElementType > &v)
{
  ElementType result = 0;
  const size_t dim = v.size();
  for(size_t i = 0; i < dim; ++i)
  {
    result += sqrt(sqrt(v[i]));
  }
  if(dim > 0)
    return result / static_cast< ElementType >(dim);
  else
    return result;
}
///////////////////////////////////////////////////////////////////////////
template< typename ElementType, int DIM >
static ElementType
weighted_mean(const ElementType v[DIM], const ElementType weights[DIM])
{
  ElementType result = 0;
  ElementType sum_weights = 0;
  for(size_t i = 0; i < DIM; ++i)
  {
    result += v[i] * weights[i];
    sum_weights += weights[i];
  }
  if(fabs(sum_weights) > std::numeric_limits< ElementType >::epsilon())
    return result / sum_weights;
  else
    return result;
}
template< typename ElementType >
static ElementType
weighted_mean(const std::vector< ElementType > &v,
              const std::vector< ElementType > &weights)
{
  ElementType result = 0;
  ElementType sum_weights = 0;
  const size_t dim = ((v.size() <= weights.size()) ? v.size() : weights.size());
  for(size_t i = 0; i < dim; ++i)
  {
    result += v[i] * weights[i];
    sum_weights += weights[i];
  }
  if(fabs(sum_weights) > std::numeric_limits< ElementType >::epsilon())
    return result / sum_weights;
  else
    return result;
}

/// Compute the median: for the time being, very costly method because we use a
/// sorting algorithm (and thus we need to get a copy of the original vector...
/// If more efficient computation is needed, please update that piece of code by
/// using the quickselect algorithm
template< class ElementType >
static ElementType
median(std::vector< ElementType > v)
{
  std::sort(v.begin(), v.end());
  if(v.size() % 2 == 0)
  {
    return static_cast< ElementType >(
        (v[v.size() / 2] + v[(v.size() / 2) - 1]) /
        2.0f); // in case if someone use an int type for ElementType...
  }
  else
  {
    return v[v.size() / 2];
  }
}

template< class ElementType >
static ElementType
percentile(std::vector< ElementType > v, float k)
{
  std::sort(v.begin(), v.end());
  if(k > 1.f)
    k = 1.f;
  else if(k < 0.f)
    k = 0.f;
  return v[(int)floor(k * (float)v.size())];
}
///////////////////////////////////////////////////////////////////////////
template< class T >
struct IndexCmp
{
  T V;
  IndexCmp(const T &v) : V(v) {}
  bool operator()(const unsigned int a, const unsigned int b) const
  {
    return V[a] < V[b];
  }
};
/*
template <class T>
inline unsigned int* sortArrayIndices(T* A, size_t size) {
        vector<unsigned int> VI(size);
        for (size_t i = 0; i < size; i++) {
                VI[i] = i;
        }
        vector<T> V(A, A+size);
        std::sort( VI.begin(), VI.end(), IndexCmp<vector<T>&>(V) );
        unsigned int *I = new unsigned int[size];
        for (size_t i = 0; i < size; i++) {
                I[i] = VI[i];
        }
        return I;
}*/

template< class T >
inline std::vector< unsigned int >
sort_vector_indices(std::vector< T > &v)
{
  std::vector< unsigned int > vi(v.size());
  for(size_t i = 0; i < v.size(); i++)
  {
    vi[i] = i;
  }
  std::sort(vi.begin(), vi.end(), IndexCmp< std::vector< T > & >(v));
  return vi;
}
template< class ElementType >
static ElementType
weightedmedian(const std::vector< ElementType > &v,
               const std::vector< ElementType > &weights)
{
  std::vector< unsigned int > vi = sort_vector_indices(v);
  double sumw = 0.0, csumw = 0.0;
  size_t i, nb_e = weights.size();
  for(i = 0; i < nb_e; ++i)
  {
    sumw += (double)weights[i];
  }
  for(i = 0; i < nb_e; ++i)
  {
    csumw += (double)weights[vi[i]] / (double)sumw;
    if(csumw >= .5f)
      break;
  }
  return v[vi[i]];
}
///////////////////////////////////////////////////////////////////////////
template< class ElementType >
static ElementType
variance(const std::vector< ElementType > &v,
         ElementType mean_value,
         bool unbiased = true)
{
  ElementType result = 0;
  size_t nb_e = v.size();
  if(nb_e <= 1)
  {
    std::cerr << "Vector::Stats::variance: at least 2 elements are needed for "
                 "computing the unbiased variance."
              << std::endl;
  }
  else
  {
    for(size_t i = 0; i < nb_e; ++i)
    {
      result += (v[i] - mean_value) * (v[i] - mean_value);
    }
    if(unbiased)
      result /= (static_cast< ElementType >(nb_e) - 1);
    else
      result /= static_cast< ElementType >(nb_e);
  }
  return result;
}

template< class ElementType >
inline ElementType
skewness(const std::vector< ElementType > &v,
         ElementType mean_value,
         ElementType biased_variance_value)
{
  ElementType result = 0;
  size_t nb_e = v.size();
  for(size_t i = 0; i < nb_e; ++i)
  {
    result += (v[i] - mean_value) * (v[i] - mean_value) * (v[i] - mean_value);
  }
  result /= static_cast< ElementType >(nb_e);
  result /= sqrt(biased_variance_value * biased_variance_value *
                     biased_variance_value +
                 static_cast< float >(1e-30));
  return result;
}

template< class ElementType >
inline ElementType
kurtosis(const std::vector< ElementType > &v,
         ElementType mean_value,
         ElementType biased_variance_value)
{
  ElementType result = 0;
  size_t nb_e = v.size();
  for(size_t i = 0; i < nb_e; ++i)
  {
    result += (v[i] - mean_value) * (v[i] - mean_value) * (v[i] - mean_value) *
              (v[i] - mean_value);
  }
  result /= (nb_e * biased_variance_value * biased_variance_value +
             static_cast< float >(1e-30));
  return result;
}
} // namespace Stats
///////////////////////////////////////////////////////////////////////////
/// Compute v1 * v2 scalar product
template< typename ElementType, int DIM >
static ElementType
dot_product(const ElementType v1[DIM], const ElementType v2[DIM])
{
  ElementType s = 0;
  for(int i = 0; i < DIM; ++i)
    s += v1[i] * v2[i];
  return s;
}
/// Compute v1 * v2 scalar product
template< typename ElementType >
static ElementType
dot_product(const std::vector< ElementType > &v1,
            const std::vector< ElementType > &v2)
{
  ElementType s = 0;
  size_t nb_e = v1.size();
  if(nb_e > v2.size())
    nb_e = v2.size();

  for(size_t i = 0; i < nb_e; ++i)
    s += v1[i] * v2[i];
  return s;
}

template< typename GeometryTraits >
static typename GeometryTraits::Scalar
dot_product(const typename GeometryTraits::Vector &v1,
            const typename GeometryTraits::Vector &v2,
            const GeometryTraits &gt)
{
  return gt.dot_product(v1, v2);
}

/// Compute v1 x v2 cross product.
/// Usually for such a function one can expect to get a return value as a
/// "ElementType*" (dynamic allocation via new), however, with the purpose to
/// avoid as many as possible unfreed chunk of memory, we rather use
/// std::vector.
template< typename ElementType, int DIM >
static std::vector< ElementType >
cross_product(const ElementType v1[DIM], const ElementType v2[DIM])
{
  size_t nb_e = 3; // the bilinear cross product exists only in 3D and 7D (but
                   // we do not consider 7D for the time being)
  if(nb_e != DIM)
    assert(false);

  std::vector< ElementType > res(nb_e);

  for(size_t i = 0; i < nb_e; ++i)
    res[i] = v1[(1 + i) % nb_e] * v2[(2 + i) % nb_e] -
             v1[(2 + i) % nb_e] * v2[(1 + i) % nb_e];

  return res;
}

/// Compute v1 x v2 cross product
template< typename ElementType >
static std::vector< ElementType >
cross_product(const std::vector< ElementType > &v1,
              const std::vector< ElementType > &v2)
{
  size_t nb_e = 3; // the bilinear cross product exists only in 3D and 7D (but
                   // we do not consider 7D for the time being)
  if((nb_e != v1.size()) || (nb_e != v2.size()))
    assert(false);

  std::vector< ElementType > res(nb_e);

  for(size_t i = 0; i < nb_e; ++i)
    res[i] = v1[(1 + i) % nb_e] * v2[(2 + i) % nb_e] -
             v1[(2 + i) % nb_e] * v2[(1 + i) % nb_e];

  return res;
}

template< typename GeometryTraits >
static typename GeometryTraits::Vector
cross_product(const typename GeometryTraits::Vector &v1,
              const typename GeometryTraits::Vector &v2,
              const GeometryTraits &gt)
{
  return gt.cross_product(v1, v2);
}

/// Compute v1 - v2 (subtraction)
template< typename ElementType, int DIM >
static std::vector< ElementType >
sub(const ElementType v1[DIM], const ElementType v2[DIM])
{
  std::vector< ElementType > res(DIM);

  for(size_t i = 0; i < DIM; ++i)
    res[i] = v1[i] - v2[i];

  return res;
}

/// Compute v1 - v2 (subtraction)
template< typename ElementType >
static std::vector< ElementType >
sub(const std::vector< ElementType > &v1, const std::vector< ElementType > &v2)
{
  size_t nb_e = v1.size();
  if(nb_e > v2.size())
    nb_e = v2.size();

  std::vector< ElementType > res(nb_e);

  for(size_t i = 0; i < nb_e; ++i)
    res[i] = v1[i] - v2[i];

  return res;
}

/// Compute p1 - p2 (subtraction)
template< typename GeometryTraits >
static typename GeometryTraits::Vector
sub(const typename GeometryTraits::Point &p1,
    const typename GeometryTraits::Point &p2,
    const GeometryTraits &gt)
{
  return gt.sub(p1, p2);
}

/// Compute v1 + v2 (addition)
template< typename ElementType, int DIM >
static std::vector< ElementType >
add(const ElementType v1[DIM], const ElementType v2[DIM])
{
  std::vector< ElementType > res(DIM);

  for(size_t i = 0; i < DIM; ++i)
    res[i] = v1[i] + v2[i];

  return res;
}

/// Compute v1 + v2 (addition)
template< typename ElementType >
static std::vector< ElementType >
add(const std::vector< ElementType > &v1, const std::vector< ElementType > &v2)
{
  size_t nb_e = v1.size();
  if(nb_e > v2.size())
    nb_e = v2.size();

  std::vector< ElementType > res(nb_e);

  for(size_t i = 0; i < nb_e; ++i)
    res[i] = v1[i] + v2[i];

  return res;
}

/// Compute p + v (addition)
template< typename GeometryTraits >
static typename GeometryTraits::Point
add(const typename GeometryTraits::Point &p,
    const typename GeometryTraits::Vector &v,
    const GeometryTraits &gt)
{
  return gt.add_p(p, v);
}

/// Compute ||v1 - v2||_L2 norm (distance between v1 and v2)
template< typename ElementType, int DIM >
static double
l2_distance(const ElementType v1[DIM], const ElementType v2[DIM])
{


  std::vector< ElementType > sub_res(Vector::sub< ElementType, DIM >(v1, v2));
  return sqrt(dot_product< ElementType >(sub_res, sub_res));
}

/// Compute ||V||_L2 norm (distance between v1 and v2)
template< typename ElementType, int DIM >
static double
l2_distance(const ElementType v[DIM])
{
  return sqrt(dot_product< ElementType, DIM >(v, v));
}

/// Compute ||p1 - p2||_L2 norm (distance between p1 and p2)
template< typename GeometryTraits >
static double
l2_distance(const typename GeometryTraits::Point &p1,
            const typename GeometryTraits::Point &p2,
            const GeometryTraits &gt)
{
  return gt.length(p1, p2);
}

/// Compute ||v1 - v2||_L2 norm (distance between v1 and v2)
template< typename ElementType >
static double
l2_distance(const std::vector< ElementType > &v1,
            const std::vector< ElementType > &v2)
{
  std::vector< ElementType > sub_res(sub(v1, v2));
  return sqrt(dot_product< ElementType >(sub_res, sub_res));
}

/// Compute ||V||_L2 norm (distance between v1 and v2)
template< typename GeometryTraits >
static double
l2_distance(const typename GeometryTraits::Vector &v, const GeometryTraits &gt)
{
  return gt.length(v);
}

/// Compute ||V||_L2 norm (distance between v1 and v2)
template< typename ElementType >
static double
l2_distance(const std::vector< ElementType > &v)
{
  return sqrt(dot_product< ElementType >(v, v));
}

/// normalize a Vector: returns a pair composed of normalized input vector and
/// of input vector length
/// Input Vector V is not normalized
template< typename GeometryTraits >
static typename GeometryTraits::Vector
normalize(const typename GeometryTraits::Vector &v, const GeometryTraits &gt)
{
  return gt.normalize(v);
}

/// Tells if 2 nD vectors are collinear or not
template< typename ElementType >
static bool
are_collinear(const std::vector< ElementType > &v1, /// equals to P2 - P1
              const std::vector< ElementType > &v2  /// equals to P1 - P3
)
{ /*
         // follow algorithm presented here:
     http://math.stackexchange.com/questions/208577/find-if-three-points-in-3-dimensional-space-are-collinear
         // yet this algorithm does not seam to work properly (some non aligned
     points are detected aligned) ElementType A = 0, B=0, C=0 ; size_t
     nbE=V1.size() ; if(V2.size()<nbE) nbE=V2.size() ;

         for(size_t i=0; i<nbE; ++i)
         {
                 A += V1[i]*V1[i] ;
                 B += V2[i]*V1[i] ;
                 C += V2[i]*V2[i] ;
         }
         if ((fabs(B) < MACH_EPS_FLOAT)) // the 2 vectors are orthogonal
                 return false;
         B *= 2 ;
         return (fabs(B*B - 4 * A*C)< MACH_EPS_FLOAT) ;*/
  return Vector::l2_distance(Vector::cross_product< ElementType >(v1, v2)) <
         MACH_EPS_FLOAT;
}

template< typename GeometryTraits >
static bool
are_collinear(const typename GeometryTraits::Vector &v1, /// equals to P2 - P1
              const typename GeometryTraits::Vector &v2  /// equals to P1 - P3
)
{
  return are_collinear(
      std::vector< typename GeometryTraits::Scalar >{v1[0], v1[1], v1[2]},
      std::vector< typename GeometryTraits::Scalar >{v2[0], v2[1], v2[2]});
}

/// Tells if 3 nD points are aligned or not
template< typename ElementType >
static bool
are_aligned(const std::vector< ElementType > &p1, /// the central point
            const std::vector< ElementType > &p2,
            const std::vector< ElementType > &p3)
{
  return are_collinear< ElementType >(sub< ElementType >(p2, p1),
                                      sub< ElementType >(p1, p3));
}

/// Compute v1 * coef
template< typename ElementType >
static std::vector< ElementType >
scalar_mult(const std::vector< ElementType > &v1, ElementType coef)
{
  size_t nb_e = v1.size();

  std::vector< ElementType > res(nb_e);

  for(size_t i = 0; i < nb_e; ++i)
    res[i] = v1[i] * coef;

  return res;
}

template< typename GeometryTraits >
static typename GeometryTraits::Vector
scalar_mult(const typename GeometryTraits::Vector &v,
            double coef,
            const GeometryTraits &gt)
{
  return gt.scalar_mult(v, coef);
}
///////////////////////////////////////////////////////////////////////////
/// V1 and V2 must be unit vectors!
template< typename ElementType, int DIM >
static double
get_angle_from_unit_vectors(const ElementType v1[DIM],
                            const ElementType v2[DIM])
{
  double cosphi = Vector::dot_product< ElementType, DIM >(v1, v2);
  double sinphi =
      Vector::l2_distance(Vector::cross_product< ElementType, DIM >(v1, v2));

  return std::atan2(sinphi, cosphi);
}
template< typename ElementType, int DIM >
static double
get_angle_in_degree_from_unit_vectors(const ElementType v1[DIM],
                                      const ElementType v2[DIM])
{
  return rad2deg(
      (Vector::get_angle_from_unit_vectors< ElementType, DIM >(v1, v2)));
}

/// V1 and V2 must be unit vectors!
template< typename ElementType >
static double
get_angle_from_unit_vectors(const std::vector< ElementType > &v1,
                            const std::vector< ElementType > &v2)
{
  double cosphi = Vector::dot_product< ElementType >(v1, v2);
  double sinphi =
      Vector::l2_distance(Vector::cross_product< ElementType >(v1, v2));

  return std::atan2(sinphi, cosphi);
}
template< typename ElementType >
static double
get_angle_in_degree_from_unit_vectors(const std::vector< ElementType > &v1,
                                      const std::vector< ElementType > &v2)
{
  return rad2deg((Vector::get_angle_from_unit_vectors< ElementType >(v1, v2)));
}
///////////////////////////////////////////////////////////////////////////
/// V1 and V2 does not need to be unit vectors!
template< class ElementType >
static double
get_angle_from_non_unit_vectors(const std::vector< ElementType > &v1,
                                const std::vector< ElementType > &v2)
{
  std::vector< ElementType > copy_v1(v1);
  double lv1 = std::sqrt(dot_product(copy_v1, copy_v1));
  if(lv1 > MACH_EPS_FLOAT)
    copy_v1 = scalar_mult(copy_v1, 1.f / lv1);
  else
    copy_v1 = std::vector< ElementType >(copy_v1.size(), 0);

  std::vector< ElementType > copy_v2(v2);

  double lv2 = std::sqrt(dot_product(copy_v2, copy_v2));
  if(lv2 > MACH_EPS_FLOAT)
    copy_v2 = scalar_mult(copy_v2, 1.f / lv2);
  else
    copy_v2 = std::vector< ElementType >(copy_v2.size(), 0);

  return get_angle_from_unit_vectors(copy_v1, copy_v2);
}
template< typename ElementType >
static double
get_angle_in_degree_from_non_unit_vectors(const std::vector< ElementType > &v1,
                                          const std::vector< ElementType > &v2)
{
  return rad2deg(
      (Vector::get_angle_from_non_unit_vectors< ElementType >(v1, v2)));
}
template< class ElementType >
static double
get_angle_from3positions(
    const std::vector< ElementType > &p1,
    const std::vector< ElementType > &p2, // point onto which angle is computed
    const std::vector< ElementType > &p3)
{
  std::vector< ElementType > v1 = sub(p1, p2);
  std::vector< ElementType > v2 = sub(p3, p2);

  return get_angle_from_non_unit_vectors(v1, v2);
}
template< typename ElementType >
static double
get_angle_in_degree_from3positions(
    const std::vector< ElementType > &p1,
    const std::vector< ElementType > &p2, // point onto which angle is computed
    const std::vector< ElementType > &p3)
{
  return rad2deg((Vector::get_angle_from3positions< ElementType >(p1, p2, p3)));
}
///////////////////////////////////////////////////////////////////////////
// Rotate a coordinate system to be perpendicular to the given normal
// ElementType should be either float, double or long double
template< typename ElementType >
static void
rot_coord_sys(
    const std::vector< ElementType >
        &old_u, /// input old_u which is assumed to be unit-length and
                /// perpendicular to old_v
    const std::vector< ElementType >
        &old_v, /// input old_v which is assumed to be unit-length and
                /// perpendicular to old_u
    const std::vector< ElementType > &new_norm, /// input new normal
    std::vector< ElementType >
        &new_u, /// output new_u which is unit-length and perpendicular to new_v
    std::vector< ElementType >
        &new_v /// output new_v which is unit-length and perpendicular to new_u
)
{
  new_u = old_u;
  new_v = old_v;
  std::vector< ElementType > old_norm = cross_product(old_u, old_v);
  ElementType ndot =
      static_cast< ElementType >(dot_product(old_norm, new_norm));
  if(ndot <= -1.0f)
  {
    new_u = -new_u;
    new_v = -new_v;
    return;
  }
  std::vector< ElementType > perp_old =
      sub(new_norm, scalar_mult(old_norm, ndot));
  std::vector< ElementType > dperp =
      scalar_mult(add(old_norm, new_norm), 1.0f / (1 + ndot));
  new_u =
      sub(new_u,
          scalar_mult(
              dperp, static_cast< ElementType >(dot_product(new_u, perp_old))));
  new_v =
      sub(new_v,
          scalar_mult(
              dperp, static_cast< ElementType >(dot_product(new_v, perp_old))));
}
///////////////////////////////////////////////////////////////////////////
// Reproject a curvature tensor from the basis spanned by old_u and old_v
// (which are assumed to be unit-length and perpendicular) to the
// new_u, new_v basis.
// ElementType should be either float, double or long double
template< typename ElementType >
void
proj_curv(const std::vector< ElementType > &old_u, /// input
          const std::vector< ElementType > &old_v, /// input
          ElementType old_ku,                      /// input
          ElementType old_kuv,                     /// input
          ElementType old_kv,                      /// input
          const std::vector< ElementType > &new_u, /// input
          const std::vector< ElementType > &new_v, /// input
          ElementType &new_ku,                     /// output
          ElementType &new_kuv,                    /// output
          ElementType &new_kv                      /// output
)
{
  std::vector< ElementType > r_new_u, r_new_v;
  rot_coord_sys(new_u, new_v, cross_product(old_u, old_v), r_new_u, r_new_v);

  ElementType u1 = static_cast< ElementType >(dot_product(r_new_u, old_u));
  ElementType v1 = static_cast< ElementType >(dot_product(r_new_u, old_v));
  ElementType u2 = static_cast< ElementType >(dot_product(r_new_v, old_u));
  ElementType v2 = static_cast< ElementType >(dot_product(r_new_v, old_v));
  new_ku = old_ku * u1 * u1 + old_kuv * (2.0f * u1 * v1) + old_kv * v1 * v1;
  new_kuv = old_ku * u1 * u2 + old_kuv * (u1 * v2 + u2 * v1) + old_kv * v1 * v2;
  new_kv = old_ku * u2 * u2 + old_kuv * (2.0f * u2 * v2) + old_kv * v2 * v2;
}
} // namespace Vector

namespace Matrix {
namespace Square {
/// Compute the 3 x 3 covariance matrix of N 3D points that are stored into a
/// one-dimensional array pointCoords [P1x, P1y, P1z, ..., PNx, PNy, PNz].
template< typename CoordinateType, size_t N >
static std::vector< CoordinateType >
covar(const CoordinateType point_coords[3 * N])
{
  std::vector< CoordinateType > cd(9); // 3 x 3 matrix

  cd[0] = 0.0f;
  for(size_t r = 0; r < N; r++)
    cd[0] += point_coords[r * 3] * point_coords[r * 3];
  for(size_t j = 0; j <= 1; j++)
  {
    cd[1 + j * 3] = 0.0f;
    for(size_t r = 0; r < N; r++)
    {
      cd[1 + j * 3] += point_coords[1 + r * 3] * point_coords[j + r * 3];
    }
  }
  for(size_t j = 0; j <= 2; j++)
  {
    cd[2 + j * 3] = 0.0f;
    for(size_t r = 0; r < N; r++)
    {
      cd[2 + j * 3] += point_coords[2 + r * 3] * point_coords[j + r * 3];
    }
  }
  cd[3] = cd[1];
  cd[6] = cd[2];
  cd[7] = cd[5];

  return cd;
}

/// Compute the 3 x 3 covariance matrix of N 3D points that are stored into a
/// one-dimensional std::vector [P1x, P1y, P1z, ..., PNx, PNy, PNz].
template< typename CoordinateType >
static std::vector< CoordinateType >
covar(const std::vector< CoordinateType > &point_coords)
{
  std::vector< CoordinateType > cd(9); // 3 x 3 matrix
  const size_t n = point_coords.size() / 3;
  cd[0] = 0.0f;
  for(size_t r = 0; r < n; r++)
    cd[0] += point_coords[r * 3] * point_coords[r * 3];
  for(size_t j = 0; j <= 1; j++)
  {
    cd[1 + j * 3] = 0.0f;
    for(size_t r = 0; r < n; r++)
    {
      cd[1 + j * 3] += point_coords[1 + r * 3] * point_coords[j + r * 3];
    }
  }
  for(size_t j = 0; j <= 2; j++)
  {
    cd[2 + j * 3] = 0.0f;
    for(size_t r = 0; r < n; r++)
    {
      cd[2 + j * 3] += point_coords[2 + r * 3] * point_coords[j + r * 3];
    }
  }
  cd[3] = cd[1];
  cd[6] = cd[2];
  cd[7] = cd[5];

  return cd;
}
///////////////////////////////////////////////////////////////////////
/// Compute coef * V x V^t (V is a vertical DIM-D vector)
template< typename ElementType, int DIM = 3 >
static void
vector_times_transpose_mult(const ElementType *p_vector,
                            ElementType pp_matrix[][DIM],
                            ElementType coeff)
{
  int j;
  for(int i = 0; i < DIM; ++i)
    for(j = 0; j < DIM; ++j)
      pp_matrix[i][j] = coeff * p_vector[i] * p_vector[j];
}
///////////////////////////////////////////////////////////////////////
/// Compute the matrix tranformation of a vector.
/// Use can use it to rotate, scale or even translate your vertex (in the later
/// case homogeneous coordinates should be used).
template< typename ElementType, int DIM = 3 >
static void
transformation(
    const ElementType pp_matrix[DIM][DIM], /// represent our square matrix that
                                           /// contains the transformation.
    const ElementType in_vertex[DIM],      /// represent the one dimensional
                                           /// untransformed column vertex
    ElementType out_vertex[DIM] /// represent the one dimensional transformed
                                /// column vertex
)
{
  int j;
  for(int i = 0; i < DIM; ++i) // for each matrix row
  {
    out_vertex[i] = 0; // sum init
    for(j = 0; j < DIM; ++j)
      out_vertex[i] += pp_matrix[i][j] * in_vertex[j];
  }
}

/// Compute the matrix tranformation of N vectors [V1 V2 V3 ... VN]
/// (column-vectors). Use can use it to rotate, scale or even translate your
/// vertex (in the later case homogeneous coordinates should be used).
template< typename ElementType, int DIM, int N >
static void
transformation(
    const ElementType pp_matrix[DIM][DIM], /// represent our square matrix that
                                           /// contains the transformation.
    const ElementType
        in_vertices[DIM][N], /// represent the untransformed column vertices
    ElementType out_vertices[DIM]
                            [N] /// represent the transformed column vertices
)
{
  int j;
  for(int i = 0; i < DIM; ++i) // for each matrix row
  {
    for(int n = 0; n < N; ++n)
    {                         // for each column vector
      out_vertices[i][n] = 0; // sum init
      for(j = 0; j < DIM; ++j)
        out_vertices[i][n] += pp_matrix[i][j] * in_vertices[j][n];
    }
  }
}

/// add two matrices
template< typename ElementType, int DIM >
static void
add(const ElementType p_matrix[][DIM], ElementType p_matrix_sum[][DIM])
{
  int j;
  for(int i = 0; i < DIM; ++i)
    for(j = 0; j < DIM; ++j)
      p_matrix_sum[i][j] += p_matrix[i][j];
}

template< typename ElementType, int DIM >
static bool
is_diagonal(const ElementType p_matrix[DIM][DIM])
{
  int j;
  for(int i = 0; i < DIM; ++i)
    for(j = 0; j < DIM; ++j)
      if((i != j) && fabs(p_matrix[i][j]) > MACH_EPS_FLOAT)
        return false;
  return true;
}


///////////////////////////////////////////////////////////////////////
} // namespace Square
} // namespace Matrix

} // namespace Math
} // namespace FEVV

