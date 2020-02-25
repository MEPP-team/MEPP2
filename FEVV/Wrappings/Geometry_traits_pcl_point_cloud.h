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

#include "FEVV/DataStructures/DataStructures_pcl_point_cloud.h"
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/Geometry_traits_operators.h"

namespace FEVV {


/**
 * \ingroup Geometry_traits_group
 * \brief   PCL Point Cloud specialization. Refer to
 *          \link
 *          RetrieveKernel_specialization_SurfaceMesh
 *          for concrete usage
 *          \endlink
 * \note    Re-create an artificial kernel to be homogeneous
 *          with other mesh implementations (CGAL, OpenMesh...)
 *          that provide a real kernel.
 */
class PCLPointCloud_kernel_generator
{
public:
  typedef PCLPointCloud_kernel_generator  Kernel;
  typedef FEVV::PCLPoint                  Point;
  typedef FEVV::PCLVector                 Vector;
  typedef FEVV::PCLKernelType             Scalar;
};

/**
 * \ingroup Geometry_traits_group
 * \anchor  RetrieveKernel_specialization_PCLPointCloud
 * \brief   PCL Point Cloud specialization. Refer to \ref RetrieveKernel for the
 *          documentation of the generic class.
 */
template<>
struct RetrieveKernel< FEVV::PCLPointCloud >
{
  typedef PCLPointCloud_kernel_generator  Kernel;
};


/**
 * \ingroup Geometry_traits_group
 * \anchor  Geometry_traits_specialization_PCLPointCloud
 * \brief   PCLPointCloud specialization of the \ref Geometry_traits
 *          generic class.
 *          For usage refer to
 *          \link Geometry_traits_documentation_dummy
 *                Geometry traits documentation
 *          \endlink
 */
template< typename MeshT >
class Geometry_traits< MeshT, PCLPointCloud_kernel_generator >
{
public:
  typedef MeshT Mesh;
  typedef Geometry_traits< MeshT, PCLPointCloud_kernel_generator > Self;
  typedef typename RetrieveKernel< Mesh >::Kernel Kernel;
  typedef typename Kernel::Point  Point;
  typedef typename Kernel::Vector Vector;
  typedef typename Kernel::Scalar Scalar;

  Geometry_traits(const Mesh &m) : m_mesh(const_cast< Mesh & >(m)) {}

  static Scalar get_x(const Point &p) { return p[0]; }

  static Scalar get_y(const Point &p) { return p[1]; }

  static Scalar get_z(const Point &p) { return p[2]; }

  static Vector normal(const Point &p1, const Point &p2, const Point &p3)
  {
    // calculate two vectors from the three points
    Vector v1 = p1 - p2; // Eigen
    Vector v2 = p1 - p3; // Eigen

    // take the cross product of the two vectors to get the normal
    return cross_product(v1, v2);
  }

  static Vector unit_normal(const Point &p1, const Point &p2, const Point &p3)
  {
    return normal(p1, p2, p3).normalized(); // Eigen
  }

  static Scalar dot_product(const Vector &v1, const Vector &v2)
  {
    return v1.dot(v2); // Eigen
  }

  static Vector cross_product(const Vector &v1, const Vector &v2)
  {
    return v1.cross(v2); // Eigen
  }

  static Scalar length2(const Vector &v) { return v.squaredNorm(); } // Eigen

  static Scalar length(const Vector &v) { return v.norm(); } // Eigen

  static Scalar length(const Point &p1, const Point &p2)
  {
    Vector v = p1 - p2; // Eigen
    return length(v);
  }

  static Vector normalize(const Vector &v)
  {
    return v.normalized(); // Eigen
  }

  static Vector add(const Vector &v1, const Vector &v2)
  {
    return v1 + v2; // Eigen
  };

  // we need addP and add functions to have function names
  // consistent with those of OpenMesh geometry trait
  static Point add_p(const Point &p, const Vector &v)
  {
    return FEVV::GeometryTraits::add_p< Self >(p, v);
  };

  // subP to be consistent with addP
  static Point sub_p(const Point &p, const Vector &v)
  {
    return FEVV::GeometryTraits::sub_p< Self >(p, v);
  };

  static Vector sub(const Point &p1, const Point &p2)
  {
    return FEVV::GeometryTraits::sub< Self >(p1, p2);
  };


  static Vector scalar_mult(const Vector &v, Scalar s)
  {
    return s * v; // Eigen
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

/**
 * \ingroup Geometry_traits_group
 * \brief Initialisation of static member NULL_VECTOR
 */
template< typename MeshT >
const typename Geometry_traits< MeshT, PCLPointCloud_kernel_generator >::Vector
    Geometry_traits< MeshT, PCLPointCloud_kernel_generator >::NULL_VECTOR(0, 0, 0);

/**
 * \ingroup Geometry_traits_group
 * \brief Initialisation of static member ORIGIN
 */
template< typename MeshT >
const typename Geometry_traits< MeshT, PCLPointCloud_kernel_generator >::Point
    Geometry_traits< MeshT, PCLPointCloud_kernel_generator >::ORIGIN(0, 0, 0);

} // namespace FEVV

