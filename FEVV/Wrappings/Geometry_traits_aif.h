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

#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/Geometry_traits_operators.h"

namespace FEVV {

/**
 * \ingroup Geometry_traits_group
 * \brief   AIF specialization. Refer to
 *          \link
 *          RetrieveKernel_specialization_SurfaceMesh
 *          for concrete usage
 *          \endlink
 * \note    Re-create an artificial kernel to be homogeneous
 *          with other mesh implementations (CGAL, OpenMesh...)
 *          that provide a real kernel.
 */
class AIF_mesh_kernel_generator
{
public:
  typedef AIF_mesh_kernel_generator Kernel;
  typedef DataStructures::AIF::AIFMesh::Point Point;
  typedef DataStructures::AIF::AIFMesh::Vector Vector;
  typedef DataStructures::AIF::AIFMesh::CoordinateType Scalar;

  // part dedicated to CGAL adaptation without including CGAL types
  typedef typename Point::CoordinateType FT;
  typedef Point Point_3;
  typedef Vector Vector_3;
  typedef typename Point::SuperClass::const_iterator Cartesian_const_iterator_3;
  typedef typename Point::SuperClass::iterator Cartesian_iterator_3;
  typedef Cartesian_const_iterator_3 Construct_cartesian_const_iterator_3;
};

/**
 * \ingroup Geometry_traits_group
 * \anchor  RetrieveKernel_specialization_AIF
 * \brief   AIF specialization. Refer to \ref RetrieveKernel for the
 *          documentation of the generic class.
 */
template<>
struct RetrieveKernel< DataStructures::AIF::AIFMesh >
{
  typedef AIF_mesh_kernel_generator Kernel;
};


/**
 * \ingroup Geometry_traits_group
 * \anchor  Geometry_traits_specialization_AIF
 * \brief   AIF specialization of the \ref Geometry_traits generic class.
 *          For usage refer to
 *          \link Geometry_traits_documentation_dummy
 *                Geometry traits documentation
 *          \endlink
 */
template< typename MeshT >
class Geometry_traits< MeshT, AIF_mesh_kernel_generator >
{
public:
  typedef MeshT Mesh;
  typedef typename RetrieveKernel< Mesh >::Kernel Kernel;
  typedef typename Kernel::Point Point;
  typedef typename Kernel::Vector Vector;
  typedef typename Kernel::Scalar Scalar;
  typedef typename FEVV::DataStructures::AIF::AIFPropertiesHelpers PropHelpers;

  Geometry_traits(const Mesh &m) : m_mesh(const_cast< Mesh & >(m)) {}

  template< int D >
  static Scalar get(const Point &p)
  {
    return p[D];
  }

  static Scalar get_x(const Point &p) { return get< 0 >(p); }

  static Scalar get_y(const Point &p) { return get< 1 >(p); }

  static Scalar get_z(const Point &p) { return get< 2 >(p); }

  static Vector unit_normal(const Point &p1, const Point &p2, const Point &p3)
  {
    return PropHelpers::compute_unit_normal(p1, p2, p3);
  }

  static Vector normal(const Point &p1, const Point &p2, const Point &p3)
  {
    return PropHelpers::compute_normal(p1, p2, p3);
  }

  static Scalar dot_product(const Vector &v1, const Vector &v2)
  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  }

  static Vector cross_product(const Vector &v1, const Vector &v2)
  {
    Vector res(0, 0, 0);
    size_t nb_e = 3;
    for(size_t i = 0; i < nb_e; ++i)
      res[i] = v1[(1 + i) % nb_e] * v2[(2 + i) % nb_e] -
               v1[(2 + i) % nb_e] * v2[(1 + i) % nb_e];

    return res;
  }

  static Scalar length2(const Vector &v)
  {
    return dot_product(v, v);
    ;
  }

  static Scalar length(const Vector &v) { return sqrt(length2(v)); }

  static Scalar length(const Point &p1, const Point &p2)
  {
    Vector v(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
    return length(v);
  }

  static Vector normalize(const Vector &v)
  {
    Scalar dist = length(v);
    Vector res;
    if(dist > 2e-7)
    {
      res = Vector(v[0] / dist, v[1] / dist, v[2] / dist);
    }
    else
      res = v;
    return res;
  }

  static Vector add(const Vector &v1, const Vector &v2)
  {
    Vector result(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
    return result;
  };

  static Point add_p(
      const Point &p,
      const Vector &v) // we need addP and add functions to have function names
                       // consistent with those of OpenMesh geometry trait
  {
    Point result(p[0] + v[0], p[1] + v[1], p[2] + v[2]);
    return result;
  };

  static Point sub_p(const Point &p1,
                     const Vector &v) // subP to be consistent with addP
  {
    Point result(p1[0] - v[0], p1[1] - v[1], p1[2] - v[2]);
    return result;
  };

  static Vector sub(const Point &p1, const Point &p2)
  {
    Vector result(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
    return result;
  };

  static Vector scalar_mult(const Vector &v, Scalar s)
  {
    Vector result(v[0] * s, v[1] * s, v[2] * s);
    return result;
  }

  static Vector scalar_mult(Scalar s, const Vector &v)
  {
    return scalar_mult(v, s);
  }

  static const Vector NULL_VECTOR;
  static const Point ORIGIN;  

protected:
  Mesh &m_mesh;
};

// Initialisation of static member NULL_VECTOR of above AIFMesh specialization
template< typename MeshT >
const typename Geometry_traits< MeshT, AIF_mesh_kernel_generator >::Vector
    Geometry_traits< MeshT, AIF_mesh_kernel_generator >::NULL_VECTOR =
        typename Geometry_traits< MeshT >::Vector(0.0, 0.0, 0.0);

// Initialisation of static member ORIGIN of above AIFMesh specialization
template< typename MeshT >
const typename Geometry_traits< MeshT, AIF_mesh_kernel_generator >::Point
    Geometry_traits< MeshT, AIF_mesh_kernel_generator >::ORIGIN =
        typename Geometry_traits< MeshT >::Point(0.0, 0.0, 0.0);

} // namespace FEVV

