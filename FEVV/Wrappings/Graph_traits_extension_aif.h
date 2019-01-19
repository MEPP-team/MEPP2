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

