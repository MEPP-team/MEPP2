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
  typedef Eigen::Vector3f                 Vector;
    // see pcl-1.8.1/include/pcl-1.8/pcl/impl/point_types.hpp
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
  typedef typename RetrieveKernel< Mesh >::Kernel Kernel;
  typedef typename Kernel::Point  Point;
  typedef typename Kernel::Vector Vector;
  typedef typename Kernel::Scalar Scalar;

  Geometry_traits(const Mesh &m) : m_mesh(const_cast< Mesh & >(m)) {}

  static Scalar get_x(const Point &p) { return p.x; }

  static Scalar get_y(const Point &p) { return p.y; }

  static Scalar get_z(const Point &p) { return p.z; }

protected:
  Mesh &m_mesh;
};


} // namespace FEVV

