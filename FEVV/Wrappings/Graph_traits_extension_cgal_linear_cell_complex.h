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

/*
 * Specialization of graph traits extension
 * for CGAL::Linear_cell_complex
 */

#include "Graph_traits_extension.h"
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#define CGAL_LCC_TEMPLATE_ARGS                                                 \
  template< unsigned int d_,                                                   \
            unsigned int ambient_dim,                                          \
            class Traits_,                                                     \
            class Items_,                                                      \
            class Alloc_,                                                      \
            template< unsigned int, class, class, class, class > class CMap,   \
            class Storage_ >

#define CGAL_LCC_TYPE                                                          \
  CGAL::Linear_cell_complex_for_combinatorial_map< d_,                         \
                                                   ambient_dim,                \
                                                   Traits_,                    \
                                                   Items_,                     \
                                                   Alloc_,                     \
                                                   CMap,                       \
                                                   Storage_ >


namespace FEVV {


/**
 * \brief  Real current number of vertices of the mesh.
 *         Specialization for CGAL::Linear_cell_complex.
 *
 * \param  g   mesh
 *
 * \return  real number of vertices of the mesh
 */
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits< CGAL_LCC_TYPE >::vertices_size_type
size_of_vertices(const CGAL_LCC_TYPE &lcc)
{
  return lcc.template attributes< 0 >().size();
}


/**
 * \brief  Real current number of edges of the mesh.
 *         Specialization for CGAL::Linear_cell_complex.
 *
 * \param  g   mesh
 *
 * \return  real number of edges of the mesh
 */
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits< CGAL_LCC_TYPE >::edges_size_type
size_of_edges(const CGAL_LCC_TYPE &lcc)
{
  return lcc.number_of_darts() / 2;
}


/**
 * \brief  Real current number of halfedges of the mesh.
 *         Specialization for CGAL::Linear_cell_complex.
 *
 * \param  g   mesh
 *
 * \return  real number of halfedges  of the mesh
 */
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits< CGAL_LCC_TYPE >::halfedges_size_type
size_of_halfedges(const CGAL_LCC_TYPE &lcc)
{
  return lcc.number_of_darts();
}


/**
 * \brief  Real current number of faces of the mesh.
 *         Specialization for CGAL::Linear_cell_complex.
 *
 * \param  g   mesh
 *
 * \return  real number of faces of the mesh
 */
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits< CGAL_LCC_TYPE >::faces_size_type
size_of_faces(const CGAL_LCC_TYPE &lcc)
{
  return lcc.template attributes< 2 >().size();
}


} // namespace FEVV


// ELO note: Clear mesh
// see
//   https://doc.cgal.org/latest/Combinatorial_map/classGenericMap.html
//
//   void clear()
//        Deletes all the darts and all the attributes of the generic map.
//
// probably used by
//   void clear(FaceGraph& g)
//   see CGAL-4.x/include/CGAL/boost/graph/helpers.h
//


#undef CGAL_LCC_TEMPLATE_ARGS
#undef CGAL_LCC_TYPE


