// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

/*
 * Specialization of graph traits extension
 * for CGAL::Polyhedron_3
 */

#include "Graph_traits_extension.h"
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

namespace FEVV {


/**
 * \brief  Real current number of vertices of the mesh.
 *         Specialization for CGAL::Polyhedron_3.
 *
 * \param  g   mesh
 *
 * \return  real number of vertices of the mesh
 */
template< typename T1, // class PolyhedronTraits_3
          typename T2, // class PolyhedronItems_3 =
                       //   CGAL::Polyhedron_items_3,
          template< class A,
                    class B,
                    class C > class T3, // class T_HDS = HalfedgeDS_default
          typename T4 >                 // class Alloc = CGAL_ALLOCATOR(int)
typename boost::graph_traits<
    CGAL::Polyhedron_3< T1, T2, T3, T4 > const >::vertices_size_type
size_of_vertices(const CGAL::Polyhedron_3< T1, T2, T3, T4 > &g)
{
  return g.size_of_vertices();
}


/**
 * \brief  Real current number of edges of the mesh.
 *         Specialization for CGAL::Polyhedron_3.
 *
 * \param  g   mesh
 *
 * \return  real number of edges of the mesh
 */
template< typename T1, // class PolyhedronTraits_3
          typename T2, // class PolyhedronItems_3 =
                       //   CGAL::Polyhedron_items_3,
          template< class A,
                    class B,
                    class C > class T3, // class T_HDS = HalfedgeDS_default
          typename T4 >                 // class Alloc = CGAL_ALLOCATOR(int)
typename boost::graph_traits<
    CGAL::Polyhedron_3< T1, T2, T3, T4 > const >::edges_size_type
size_of_edges(const CGAL::Polyhedron_3< T1, T2, T3, T4 > &g)
{
  return g.size_of_halfedges() / 2;
}


/**
 * \brief  Real current number of halfedges of the mesh.
 *         Specialization for CGAL::Polyhedron_3.
 *
 * \param  g   mesh
 *
 * \return  real number of halfedges of the mesh
 */
template< typename T1, // class PolyhedronTraits_3
          typename T2, // class PolyhedronItems_3 =
                       //   CGAL::Polyhedron_items_3,
          template< class A,
                    class B,
                    class C > class T3, // class T_HDS = HalfedgeDS_default
          typename T4 >                 // class Alloc = CGAL_ALLOCATOR(int)
typename boost::graph_traits<
    CGAL::Polyhedron_3< T1, T2, T3, T4 > const >::halfedges_size_type
size_of_halfedges(const CGAL::Polyhedron_3< T1, T2, T3, T4 > &g)
{
  return g.size_of_halfedges();
}


/**
 * \brief  Real current number of faces of the mesh.
 *         Specialization for CGAL::Polyhedron_3.
 *
 * \param  g   mesh
 *
 * \return  real number of faces of the mesh
 */
template< typename T1, // class PolyhedronTraits_3
          typename T2, // class PolyhedronItems_3 =
                       //   CGAL::Polyhedron_items_3,
          template< class A,
                    class B,
                    class C > class T3, // class T_HDS = HalfedgeDS_default
          typename T4 >                 // class Alloc = CGAL_ALLOCATOR(int)
typename boost::graph_traits<
    CGAL::Polyhedron_3< T1, T2, T3, T4 > const >::faces_size_type
size_of_faces(const CGAL::Polyhedron_3< T1, T2, T3, T4 > &g)
{
  return g.size_of_facets();
}


} // namespace FEVV


// ELO note: Clear mesh
// see
//   http://doc.cgal.org/latest/Polyhedron/classCGAL_1_1Polyhedron__3.html
//
//   void clear()
//   removes all vertices, halfedges, and facets.

