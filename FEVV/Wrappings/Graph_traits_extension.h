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

#include <boost/graph/graph_traits.hpp>
#include <iterator> // for std::distance

namespace FEVV {


/**
 * \brief  Real current number of vertices of the mesh.
 *         Generic version.
 *
 * \param  g   mesh
 *
 * \return  real number of vertices of the mesh
 */
template< typename MeshT >
typename boost::graph_traits< MeshT >::vertices_size_type
size_of_vertices(const MeshT &g)
{
  return (typename boost::graph_traits< MeshT >::vertices_size_type)
      std::distance(vertices(g).first, vertices(g).second);
}


/**
 * \brief  Real current number of edges of the mesh.
 *         Generic version.
 *
 * \param  g   mesh
 *
 * \return  real number of edges of the mesh
 */
template< typename MeshT >
typename boost::graph_traits< MeshT >::edges_size_type
size_of_edges(const MeshT &g)
{
  return (typename boost::graph_traits< MeshT >::edges_size_type)std::distance(
      edges(g).first, edges(g).second);
}


/**
 * \brief  Real current number of halfedges of the mesh.
 *         Generic version.
 *
 * \param  g   mesh
 *
 * \return  real number of halfedges of the mesh
 */
template< typename MeshT >
typename boost::graph_traits< MeshT >::halfedges_size_type
size_of_halfedges(const MeshT &g)
{
  return (typename boost::graph_traits< MeshT >::halfedges_size_type)
      std::distance(halfedges(g).first, halfedges(g).second);
}


/**
 * \brief  Real current number of faces of the mesh.
 *         Generic version.
 *
 * \param  g   mesh
 *
 * \return  real number of faces of the mesh
 */
template< typename MeshT >
typename boost::graph_traits< MeshT >::faces_size_type
size_of_faces(const MeshT &g)
{
  return (typename boost::graph_traits< MeshT >::faces_size_type)std::distance(
      faces(g).first, faces(g).second);
}


} // namespace FEVV

