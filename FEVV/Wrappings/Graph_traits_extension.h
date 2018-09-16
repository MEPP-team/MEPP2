#ifndef FEVV_GRAPH_TRAITS_EXTENSION_H
#define FEVV_GRAPH_TRAITS_EXTENSION_H

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

#endif //  FEVV_GRAPH_TRAITS_EXTENSION_H
