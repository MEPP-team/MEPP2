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

#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#if !defined(CGAL_USE_OM_POINTS)
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>
#endif
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/Geometry_traits_operators.h"

namespace FEVV {

/**
 * \ingroup Geometry_traits_group
 * \anchor  RetrieveKernel_specialization_OpenMesh
 * \brief   OpenMesh specialization. Refer to \ref RetrieveKernel for the
 *          documentation of the generic class.
 */
template< typename T >
struct RetrieveKernel< OpenMesh::PolyMesh_ArrayKernelT< T > >
{
  typedef T Kernel;
};

/**
 * \ingroup Geometry_traits_group
 * \anchor  Geometry_traits_specialization_OpenMesh
 * \brief   OpenMesh specialization of the \ref Geometry_traits generic class.
 *          For usage refer to
 *          \link Geometry_traits_documentation_dummy
 *                Geometry traits documentation
 *          \endlink

 */
template< typename T >
class Geometry_traits< OpenMesh::PolyMesh_ArrayKernelT< T >, T >
{
public:
  typedef OpenMesh::PolyMesh_ArrayKernelT< T > Mesh;
  typedef Geometry_traits< Mesh, T > Self;
  typedef typename Mesh::AttribKernel Kernel;
  typedef typename Kernel::Normal Vector;
  typedef typename Kernel::Scalar Scalar;

// When CGAL_USE_OM_POINTS is defined then OpenMesh CGAL wrappers (refer to
//    CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h) build on the fly
// of the put( vertex_point_t ) / get( vertex_point_t ) methods property
// maps that do not use OpenMesh kernel provided Point type. Instead the
// property maps use some hardwired CGAL Point type. In order for this
// Geometry_traits specialized for OpenMesh to be compatible with this
// CGAL choice for the boost/graph wrappers we adapt the Geometry Traits
// accordingly
#if defined(CGAL_USE_OM_POINTS)
  typedef typename Kernel::Point Point;
#else
  typedef CGAL::Exact_predicates_inexact_constructions_kernel CGALKernel;
  typedef CGALKernel::Point_3 Point;
  typedef CGALKernel::Vector_3 CGALVector;
#endif

  Geometry_traits(const Mesh &m) : m_mesh(const_cast< Mesh & >(m)) {}

  /// Deprecated syntax for coordinate access. Use the get_[x|y|z] versions
  template< int D >
  static Scalar get(const Point &p)
  {
    return static_cast< Scalar >(p[D]);
  }

  static Scalar get_x(const Point &p) { return static_cast< Scalar >(p[0]); }

  static Scalar get_y(const Point &p) { return static_cast< Scalar >(p[1]); }

  static Scalar get_z(const Point &p) { return static_cast< Scalar >(p[2]); }

  Vector unit_normal(const Point &p1, const Point &p2, const Point &p3) const
  {
#if defined(CGAL_USE_OM_POINTS)
    return m_mesh.calc_face_normal(p1, p2, p3).normalize();
#else
    CGALVector cga_lresult = CGAL::unit_normal(p1, p2, p3);
    return Vector(cga_lresult[0], cga_lresult[1], cga_lresult[2]);
#endif
  }

  Vector normal(const Point &p1, const Point &p2, const Point &p3) const
  {
#if defined(CGAL_USE_OM_POINTS)
    return m_mesh.calc_face_normal(p1, p2, p3);
#else
    CGALVector cga_lresult = CGALKernel::Construct_normal_3()(p1, p2, p3);
    return Vector(cga_lresult[0], cga_lresult[1], cga_lresult[2]);
#endif
  }

  static Scalar dot_product(const Vector &v1, const Vector &v2)
  {
    return v1 | v2;
  }

  static Vector cross_product(const Vector &v1, const Vector &v2)
  {
    return v1 % v2;
  }

  static Scalar length2(const Vector &v) { return v.sqrnorm(); }

  static Scalar length(const Vector &v) { return v.length(); }

  static Scalar length(const Point &p1, const Point &p2)
  {
    // here,  we have have native OpenMesh Points (VectorT) or native CGAL
    // Vector_3
    auto v = p1 - p2;
#if defined(CGAL_USE_OM_POINTS)
    return length(v);
#else
    Scalar dotprodloc = 0;
    for(int i = 0; i < v.dimension(); ++i)
      dotprodloc += v[i] * v[i];
    return sqrt(dotprodloc);
#endif
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

  static Point
  add_p(const Point &p,
        const Vector &v) // we need add_p and add functions, beacuse in OpenMesh
                         // Point and Vector have the same type
  { // here,  we have have native OpenMesh Points (VectorT) or native CGAL
    // Vector_3
#if defined(CGAL_USE_OM_POINTS)
    return p + v;
#else
    return GeometryTraits::add_p< Self >(p, v);
#endif
  };

  static Point sub_p(const Point &p,
                     const Vector &v) // sub_p to be consistent with add_p
  { // here,  we have have native OpenMesh Points (VectorT) or native CGAL
    // Vector_3
#if defined(CGAL_USE_OM_POINTS)
    return p - v;
#else
    return GeometryTraits::sub_p< Self >(p, v);
#endif
  };

  static Vector sub(const Point &p1, const Point &p2)
  { // here,  we have have native OpenMesh Points (VectorT) or native CGAL
    // Vector_3 (and there is not automatic convertion from Vector_3 to Vector
#if defined(CGAL_USE_OM_POINTS)
    return p1 - p2;
#else
    return GeometryTraits::sub< Self >(p1, p2);
#endif
  };

  static Vector scalar_mult(const Vector &v, Scalar s) { return v * s; }

  static Vector scalar_mult(Scalar s, const Vector &v)
  {
    return scalar_mult(v, s);
  }

  static const Vector NULL_VECTOR;
  static const Point ORIGIN;
  
protected:
  Mesh &m_mesh;
};

/// Initialisation of static member NULL_VECTOR of above OpenMesh specialisation
template< typename T >
const typename Geometry_traits< OpenMesh::PolyMesh_ArrayKernelT< T >,
                                T >::Vector
    Geometry_traits< OpenMesh::PolyMesh_ArrayKernelT< T >, T >::NULL_VECTOR =
        typename Geometry_traits< OpenMesh::PolyMesh_ArrayKernelT< T >,
                                  T >::Vector(0, 0, 0);
								  
/// Initialisation of static member ORIGIN of above OpenMesh specialisation
template< typename T >
const typename Geometry_traits< OpenMesh::PolyMesh_ArrayKernelT< T >,
                                T >::Point
    Geometry_traits< OpenMesh::PolyMesh_ArrayKernelT< T >, T >::ORIGIN =
        typename Geometry_traits< OpenMesh::PolyMesh_ArrayKernelT< T >,
                                  T >::Point(0, 0, 0);								  

} // namespace FEVV

