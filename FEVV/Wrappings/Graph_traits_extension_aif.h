#ifndef FEVV_GRAPH_TRAITS_EXTENSION_AIF_H
#define FEVV_GRAPH_TRAITS_EXTENSION_AIF_H


/*
 * Specialization of graph traits extension
 * for AIF
 */

#include "Graph_traits_extension.h"
#include "Graph_traits_aif.h"
#include "FEVV/DataStructures/AIF/AIFMesh.hpp"

namespace FEVV {


/**
 * \brief  Real current number of vertices of the mesh.
 *         Specialization for AIF.
 *
 * \param  g   mesh
 *
 * \return  real number of vertices of the mesh
 */
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::vertices_size_type
size_of_vertices(const FEVV::DataStructures::AIF::AIFMesh &g)
{
  return static_cast< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::vertices_size_type >(
      g.GetNumberOfVertices());
}


/**
 * \brief  Real current number of edges of the mesh.
 *         Specialization for AIF.
 *
 * \param  g   mesh
 *
 * \return  real number of edges of the mesh
 */
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::edges_size_type
size_of_edges(const FEVV::DataStructures::AIF::AIFMesh &g)
{
  return static_cast< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::edges_size_type >(
      g.GetNumberOfEdges());
}


/**
 * \brief  Real current number of faces of the mesh.
 *         Specialization for AIF.
 *
 * \param  g   mesh
 *
 * \return  real number of faces of the mesh
 */
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::faces_size_type
size_of_faces(const FEVV::DataStructures::AIF::AIFMesh &g)
{
  return static_cast< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::faces_size_type >(
      g.GetNumberOfFaces());
}


} // namespace FEVV


#endif //  FEVV_GRAPH_TRAITS_EXTENSION_AIF_H
