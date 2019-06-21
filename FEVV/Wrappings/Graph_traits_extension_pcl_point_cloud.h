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
 * for FEVV::PCLPointCloud
 */

#include "Graph_traits_extension.h"

namespace FEVV {


/**
 * \brief  Real current number of vertices of the mesh.
 *
 * \param  g   mesh
 *
 * \return  real number of vertices of the mesh
 */
inline
typename boost::graph_traits< FEVV::PCLPointCloud >::vertices_size_type
size_of_vertices(const FEVV::PCLPointCloud &g)
{
  return num_vertices(g);
}


/**
 * \brief  Real current number of edges of the mesh.
 *
 * \param  g   mesh
 *
 * \return  real number of edges of the mesh
 */
inline
int
size_of_edges(const FEVV::PCLPointCloud &g)
{
  return 0;
}


/**
 * \brief  Real current number of halfedges of the mesh.
 *
 * \param  g   mesh
 *
 * \return  real number of halfedges of the mesh
 */
inline
int
size_of_halfedges(const FEVV::PCLPointCloud &g)
{
  return 0;
}


/**
 * \brief  Real current number of faces of the mesh.
 *
 * \param  g   mesh
 *
 * \return  real number of faces of the mesh
 */
inline
int
size_of_faces(const FEVV::PCLPointCloud &g)
{
  return 0;
}


} // namespace FEVV


// ELO note: Clear mesh
// see
//   http://doc.cgal.org/latest/Surface_mesh/classCGAL_1_1Surface__mesh.html#a247d4ad3e6b106ae22e5306203812642
//
//   void clear()
//   removes all vertices, halfedge, edges and faces. Collects garbage and
//   clears all properties.
// see
//   CGAL-4.9/include/CGAL/boost/graph/helpers.h
//


