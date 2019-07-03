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

#include "FEVV/Wrappings/Geometry_traits.h" // For forward convenience
#include "FEVV/Wrappings/Geometry_traits_operators.h"
#include <CGAL/Kernel/global_functions.h> // for CGAL::unit_normal, CGAL::scalar_product, CGAL::cross_product

namespace FEVV {

/**
 * \ingroup Geometry_traits_group
 * \brief Generic class used as a trick to circumvent C++ absence of
 *        inner member specialisation. This templated struct avoids testing
 *        the dimension argument of
 *        \link
 *          Geometry_traits_for_cgal_get_int_Point
 *          Geometry_traits_for_cgal::get< int >( Point) member
 *        \endlink
 *        at run time.
 */
template< typename Point, typename Scalar, int D >
struct CGAL_point_get
{
  static Scalar apply(const Point &p) { return 0; }
};

/*
 * \ingroup Geometry_traits_group
 * \brief Template specialisation of \ref CGAL_point_get
 */
template< typename Point, typename Scalar >
struct CGAL_point_get< Point, Scalar, 0 >
{
  static Scalar apply(const Point &p) { return p.x(); }
};

/*
 * \ingroup Geometry_traits_group
 * \brief Template specialisation of \ref CGAL_point_get
 */
template< typename Point, typename Scalar >
struct CGAL_point_get< Point, Scalar, 1 >
{
  static Scalar apply(const Point &p) { return p.y(); }
};

/*
 * \ingroup Geometry_traits_group
 * \brief Template specialisation of \ref CGAL_point_get
 */
template< typename Point, typename Scalar >
struct CGAL_point_get< Point, Scalar, 2 >
{
  static Scalar apply(const Point &p) { return p.z(); }
};

/**
 * \ingroup Geometry_traits_group
 * \anchor Geometry_traits_for_cgal
 * \brief  Many mesh representations of CGLA, e.g.
 *         CGAL::Exact_predicate_inexact_construction_kernel and
 *         CGAL::Cartesian<double> have compatible kernels (i.e. kernels
 *         exposing the same geometry related traits). Hence their respective
 *         \ref Geometry_traits template specializations will define exactly
 *         the same code. In order to avoid code replication this
 *         \ref Geometry_traits_for_cgal class "factorizes" the code that
 *         will be shared by CGAL meshes sharing such a common kernel
 *         structure. The various CGAL meshes needing simply inherit from this
 *         class like
 *          - \link
 *              Geometry_traits_specialization_SurfaceMesh
 *              Geometry_traits< Mesh , FEVV::Surface_mesh_kernel_generator<> >
 *            \endlink
 *          - \link
 *              Geometry_traits_specialization_cgal_exact_predicates_inexact
 *              Geometry_traits< Mesh,
 * CGAL::Exact_predicates_inexact_constructions_kernel > \endlink \see    The
 * \ref Geometry_traits generic class. For usage refer to \link
 * Geometry_traits_documentation_dummy Geometry traits documentation \endlink
 */
template< typename MeshT, typename KernelT >
class Geometry_traits_for_cgal : public KernelT
{
public:
  typedef Geometry_traits_for_cgal Self;
  typedef MeshT Mesh;
  typedef KernelT Kernel;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT Scalar;

  Geometry_traits_for_cgal(const Mesh &m) : m_mesh(const_cast< Mesh & >(m)) {}

  static Scalar get_x(const Point &p) { return p.x(); }

  static Scalar get_y(const Point &p) { return p.y(); }

  static Scalar get_z(const Point &p) { return p.z(); }

  /// \anchor Geometry_traits_for_cgal_get_int_Point
  template< int D >
  static Scalar get(const Point &p)
  {
    return CGAL_point_get< Point, Scalar, D >::apply(p);
  }

  static Vector unit_normal(const Point &p1, const Point &p2, const Point &p3)
  {
    return CGAL::unit_normal(p1, p2, p3);
  }

  static Vector normal(const Point &p1, const Point &p2, const Point &p3)
  {
    return typename Kernel::Construct_normal_3()(p1, p2, p3);
  }

  static Scalar dot_product(const Vector &v1, const Vector &v2)
  {
    return CGAL::scalar_product(v1, v2);
  }

  static Vector cross_product(const Vector &v1, const Vector &v2)
  {
    return CGAL::cross_product(v1, v2);
  }

  static Scalar length2(const Vector &v) { return v.squared_length(); }

  static Scalar length(const Vector &v) { return sqrt(length2(v)); }

  static Scalar length(const Point &p1, const Point &p2)
  {
    Vector v = p1 - p2;
    return length(v);
  }

  static Vector normalize(const Vector &v)
  {
    Scalar dist = length(v);
    Vector res;
    if(dist > 2e-7)
    {
      res = v * 1. / dist;
    }
    else
      res = v;
    return res;
  }

  static Vector add(const Vector &v1, const Vector &v2) { return v1 + v2; };

  static Point add_p(
      const Point &p,
      const Vector &v) // we need addP and add functions to have function names
                       // consistent with those of OpenMesh geometry trait
  {
    /*Point result(p1[0] + v[0],
                                 p1[1] + v[1],
                                 p1[2] + v[2] );
    return result;*/
    return p +
           v; // defined in
              // http://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html
  };

  static Point sub_p(const Point &p,
                     const Vector &v) // subP to be consistent with addP
  {
    /*Point result(p[0] - v[0],
                                 p[1] - v[1],
                                 p[2] - v[2] );
    return result;*/
    return p -
           v; // defined in
              // http://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html
  };

  static Vector sub(const Point &p1, const Point &p2) { return p1 - p2; };


  static Vector scalar_mult(const Vector &v, Scalar s) { return v * s; }

  static Vector scalar_mult(Scalar s, const Vector &v)
  {
    return scalar_mult(v, s);
  }

  static const Vector NULL_VECTOR;
  static const Point ORIGIN;
  
protected:
  MeshT &m_mesh;
};

/**
 * \ingroup Geometry_traits_group
 * \brief Initialisation of static member NULL_VECTOR of
 *        \ref Geometry_traits_for_cgal class.
 */
template< typename MeshT, typename KernelT >
const typename Geometry_traits_for_cgal< MeshT, KernelT >::Vector
    Geometry_traits_for_cgal< MeshT, KernelT >::NULL_VECTOR = CGAL::NULL_VECTOR;

/**
 * \ingroup Geometry_traits_group
 * \brief Initialisation of static member ORIGIN of
 *        \ref Geometry_traits_for_cgal class.
 */
template< typename MeshT, typename KernelT >
const typename Geometry_traits_for_cgal< MeshT, KernelT >::Point
    Geometry_traits_for_cgal< MeshT, KernelT >::ORIGIN = CGAL::ORIGIN;

} // namespace FEVV

