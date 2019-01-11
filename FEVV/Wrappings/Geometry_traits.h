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

namespace FEVV {

/**
 * \defgroup Geometry_traits_group Geometry traits.
 * ## Original needs
 * A generic library relies on a set of unified interfaces to which all concrete
 * native external implementations (e.g. CGAL::Polyhedron_3 or OpenMesh),
 * that will be used for instantiation, must comply to.
 * The separation of concerns induces the separation of the topology related
 * interface from the geometry related one.
 * The main purpose of Geometry_trait is to provide such an interface i.e.
 *  - to expose in an unified manner the various underlying geometry related
 *    types of the considered native external implementations (for example
 *    `OpenMesh::PolyMesh_ArrayKernelT<>` underlying kernel (`AttribKernel`)
 *    uses a `Point` type for its vertex geometry whereas
 *    `CGAL::Polyhedron_3< CGAL::Cartesian<double>, ...>` kernel uses a
 *    `Point_3` type),
 *  - to provide common functionality for those type technically implemented
 *    as generic (for the geometry) algorithms
 *
 * Emerging such a geometry interface requires to "extend" all considered
 * native implementations in order to comply with the geometry interface that
 * we wish to promote and on top of we shall provide our generic algorithms.
 * Among the possible techniques for extending the native implementation
 * (wrappers, inheritance...) the design choice was led by the additional
 * need of following the Boost Graph Library syntatic approach of so called
 * "free functions". This aesthetic need is drawn from a wish to provide
 * notations that are homogeneous with the ones chosen for the topology related
 * interface for which we chose to use the CGAL provided
 * [boost graph wrappers](http://doc.cgal.org/latest/BGL/group__PkgBGL.html).
 * And when CGAL wrapped their native mesh data stratructures (as well as
 * OpenMesh), they not only adopted some
 * [BGL topology concepts](http://doc.cgal.org/latest/BGL/annotated.html),
 * but they also committed themselves to also follow BGL syntatic approach.
 *
 * ## Design notes
 * For the syntax of free functions to remain simple, we thus need to provide
 * (traditionaly as first argument) the concerned mesh instance (e.g.
 * `rescale( myMesh, 6.0)`).
 * Technically this means that every template specialization (for the various
 * native mesh implementations we wish to provide the geometry wrappers for) of
 * the free functions we propose must be distinguished from the others. This
 * distinction must be realized with the only information provided to the free
 * functions which is the (mesh) type of their first argument.
 *
 * FEVV geometry interfaces uses two mechanisms to achieve its goal:
 *   - \ref RetrieveKernel that "extracts" the geometry related information
 *     out of a Mesh type, (using the underlying kernel when available),
 *     and its specializations:
 *      - \link
 *           RetrieveKernel_specialization_OpenMesh
 *           RetrieveKernel< OpenMesh::PolyMesh_ArrayKernelT< > >
 *        \endlink
 *      - \link
 *           RetrieveKernel_specialization_SurfaceMesh
 *           RetrieveKernel< CGAL::Surface_mesh< > >
 *        \endlink
 *      - \link
 *           RetrieveKernel_specialization_cgal_Polyhedron
 *           RetrieveKernel< CGAL::Polyhedron_3< > >
 *        \endlink
 *   - \ref Geometry_traits a trait acting as unifying interface and that
 *     wraps the underlying mesh native implementations in order to
 *     blur their discrepancies.
 *      - \link
 *           Geometry_traits_specialization_OpenMesh
 *           Geometry_traits< Mesh,
 * OpenMesh::PolyMesh_ArrayKernelT<>::AttribKernel > \endlink
 *      - \link
 *           Geometry_traits_specialization_SurfaceMesh
 *           Geometry_traits< Mesh, FEVV::Surface_mesh_kernel_generator< Point>
 * > \endlink
 *      - \link
 *           Geometry_traits_specialization_cgal_cartesian
 *           Geometry_traits< Mesh, CGAL::Cartesian<double> >
 *        \endlink
 *      - \link
 *           Geometry_traits_specialization_cgal_exact_predicates_inexact
 *           Geometry_traits< Mesh,
 * CGAL::Exact_predicates_inexact_constructions_kernel > \endlink
 */

/**
 * \ingroup Geometry_traits_group
 * \brief A generic definition, that is template specialized for every
 *        supported native implementation, allowing for the retrieval (or
 *        extraction) of the geometry related information out of a Mesh
 *        type (possibly using the respective enderlying kernels of the native
 *        implementations when available).
 *        RetrieveKernel implementatinos (through template specializations):
 * \see   \link
 *           RetrieveKernel_specialization_OpenMesh
 *           RetrieveKernel< OpenMesh::PolyMesh_ArrayKernelT< T > >
 *        \endlink
 * \see   \link
 *           RetrieveKernel_specialization_SurfaceMesh
 *           RetrieveKernel< CGAL::Surface_mesh< T > >
 *        \endlink
 * \see   \link
 *           RetrieveKernel_specialization_cgal_Polyhedron
 *           RetrieveKernel< CGAL::Polyhedron_3< < >
 *        \endlink
 */
template< typename Mesh >
struct RetrieveKernel
{
};

/**
 * \ingroup Geometry_traits_group
 * \class  Geometry_traits
 * \brief  Refer to \ref Geometry_traits_documentation_dummy for further
 *         documentation on provided types and algorithms (if for example
 *         you are looking on what you need to provide when realizing a new
 *         specialization). If you are trying to understand what this
 *         Geometry_traits class is about, just consider Geometry_traits as a
 *         generic definition whose only purpose is to provide an anchoring
 * point for a unified geometry interface. The Geometry_traits template class is
 * used through template specializations. \tparam Mesh The Mesh out of which to
 * build the Geometry_traits. This template is only needed for mesh data
 * structures that have a loose relationship with the underlying geometry
 * related types (e.g. OpenMesh). For data structures having a kernel, this is
 * where the geometry is located. Hence, technically the Mesh template argument
 * is only needed for data structures that embed their geometry related types
 * and algorithmes (e.g. like OpenMesh). \tparam Kernel The geometric kernel
 * when available. This is defaulted to the RetrieveKernel<> utility template
 * whose purpose is to re-extract the kernel from mesh data structures providing
 * one. \see    \link Geometry_traits_specialization_OpenMesh Geometry_traits<
 * Mesh, OpenMesh::PolyMesh_ArrayKernelT<>::AttribKernel > \endlink \see \link
 *           Geometry_traits_specialization_SurfaceMesh
 *           Geometry_traits< Mesh, FEVV::Surface_mesh_kernel_generator< Point>
 * > \endlink \see    \link Geometry_traits_specialization_cgal_cartesian
 *           Geometry_traits< Mesh, CGAL::Cartesian<double> >
 *         \endlink
 * \see    \link
 *           Geometry_traits_specialization_cgal_exact_predicates_inexact
 *           Geometry_traits< Mesh,
 * CGAL::Exact_predicates_inexact_constructions_kernel > \endlink
 */

template< typename Mesh,
          typename Kernel = typename RetrieveKernel< Mesh >::Kernel >
class Geometry_traits
{
};

/**
 * \ingroup Geometry_traits_group
 * \anchor  Geometry_traits_documentation_dummy
 * \brief   Documentation (as a template dummy class) of \ref Geometry_traits
 * \tparam  MeshT   refer to \ref Geometry_traits template parameters
 *                  documentation
 * \tparam  KernelT refer to \ref Geometry_traits template parameters
 *                  documentation
 */
template< typename MeshT, typename KernelT >
class Geometry_traits_documentation_dummy
{
public:
  /// Unused internaly but nice to embed for unforseen caller usage
  typedef MeshT Mesh;
  /// Unused internaly but nice to embed for unforseen caller usage
  typedef KernelT Kernel;
  typedef typename Kernel::Point Point;
  typedef typename Kernel::Normal Vector;
  typedef typename Kernel::Scalar Scalar;

  /// Coordinate accessor
  template< int D >
  Scalar get(const Point &p)
  {
    return Scalar();
  }

  /**
   * \brief Returns the unit normal of the three argument points (implemented as
   *  \link GeometryTrait_free_function_unit_normal unit_normal() free
   *  function \endlink
   */
  Vector
  unit_normal(const Mesh &m, const Point &p1, const Point &p2, const Point &p3);

  /**
   * \brief Returns the normal of the three argument points (implemented as
   *  \link GeometryTrait_free_function_normal normal() free function \endlink
   */
  Vector
  normal(const Mesh &m, const Point &p1, const Point &p2, const Point &p3);

  /**
   * \brief Returns the normal of the three argument points (implemented as
   *  \link GeometryTrait_free_function_get_x get_x() free function \endlink
   */
  Scalar get_x(const Mesh &, const Point &p1);

  /**
   * \brief Returns the normal of the three argument points (implemented as
   *  \link GeometryTrait_free_function_get_y get_y() free function \endlink
   */
  Scalar get_y(const Mesh &, const Point &p1);

  /**
   * \brief Returns the normal of the three argument points (implemented as
   *  \link GeometryTrait_free_function_get_z get_z() free function \endlink
   */
  Scalar get_z(const Mesh &, const Point &p1);

  /**
   * \brief Returns the normal of the three argument points (implemented as
   *  \link GeometryTrait_free_function_length length() free function \endlink
   */
  Scalar length(const Mesh &, Vector &v);

  /**
   * \brief Returns the vector resulting from the addition of the two
   *        argument vectors.
   */
  Vector add_v(const Vector &v1, const Vector &v2);
};

} // namespace FEVV

