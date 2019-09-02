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
#include <boost/graph/properties.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/Generic/Manifold/calculate_vertex_normal.hpp"

namespace FEVV {
namespace Filters {

/**
 * \brief Compute the respectice vertex normal for all the vertices
 *        of the considered mesh and populate the argument provided Face
 *        Normal Map with the resulting normals.
 * \sa    Consider using the other version of this function that doesn't
 *        require the Geometry_traits as argument.
 *
 * \param[in]  g    The considered mesh.
 * \param[in]  pm   The Point Map holding the vertices geometry
 * \param[in]  fnm  The Property Map providing the face normal vectors
 * \param[out] vnm  The Property Map that will be populated with the
 *                  resulting normal vectors.
 * \param[in]  gt   The Geometry Traits associated to the mesh
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename FaceNormalMap,
          typename VertexNormalMap,
          typename GeometryTraits >
void
calculate_vertex_normals(const HalfedgeGraph &g,
                         const PointMap &pm,
                         const FaceNormalMap &fnm, /// provided face normals
                         VertexNormalMap vnm,
                         const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;

  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb)
  {
    Vector vec = FEVV::Operators::calculate_vertex_normal< HalfedgeGraph,
                                                           PointMap,
                                                           FaceNormalMap,
                                                           GeometryTraits >(
        *vb, g, pm, fnm, gt);
    put(vnm, *vb, vec);
  }
}

/**
 * \brief Compute the respectice vertex normal for all the vertices
 *        of the considered mesh and populate the argument provided Face
 *        Normal Map with the resulting normals.
 * \note  This function is a syntactic sugar wrapper of the underlying
 *        version of the function that requires an additional Geometry_traits
 *        as argument.
 *
 * \param[in]  g    The considered mesh.
 * \param[in]  pm   The Point Map holding the vertices geometry
 * \param[in]  fnm  The Property Map providing the face normal vectors
 * \param[out] vnm  The Property Map that will be populated with the
 *                  resulting normal vectors.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename FaceNormalMap,
          typename VertexNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
calculate_vertex_normals(const HalfedgeGraph &g,
                         const PointMap &pm,
                         const FaceNormalMap &fnm, /// provided face normals
                         VertexNormalMap vnm)
{
  GeometryTraits gt(g);
  calculate_vertex_normals< HalfedgeGraph,
                            PointMap,
                            FaceNormalMap,
                            VertexNormalMap,
                            GeometryTraits >(g, pm, fnm, vnm, gt);
}

} // namespace Filters
} // namespace FEVV
