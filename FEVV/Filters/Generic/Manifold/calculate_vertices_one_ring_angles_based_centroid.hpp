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
#include "FEVV/Operators/Generic/Manifold/vertex_one_ring_angles_based_centroid.hpp"

namespace FEVV {
namespace Filters {

/**
 * \brief  Compute angle-based centroids for all vertices and store them in
 *         bangles_based_centroids_pm.
 *
 * \tparam FaceGraph a Mesh type that provides a Model of the
 *         FaceGraph Concept through a boost::graph_traits<> specialization.
 * \tparam PointMap A modifiable point map to manage vertex positions.
 * \tparam CentroidMap A modifiable point map to manage vertex' angle-based
 *         centroid.
 * \tparam GeometryTraits The geometric kernel when available. This is
 *         defaulted to FEVV::Geometry_traits<FaceGraph>.
 * \param  g The FaceGraph instance from which the vertices will be taken.
 * \param  pm The point map which associate a vertex descriptor with a vertex
 *         position.
 * \param  angles_based_centroids_pm  The point map to store the computed vertex
 *         centroids.
 * \param  SMOOTHING_FACTOR A (usually positive) factor used to compute the
 *         position of a stored centroid following
 *         V_pos + SMOOTHING_FACTORx(C_computed - V_pos).
 * \param  gt The geometry trait object.
 * \pre    g must be a triangular mesh.
 */
template< typename FaceGraph,
          typename PointMap,
          typename CentroidMap,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
void
calculate_vertices_one_ring_angles_based_centroid(
    const FaceGraph &g,
    const PointMap &pm, // const map
    CentroidMap
        angles_based_centroids_pm, // modifiable map (follows CGAL usage)
    const typename GeometryTraits::Scalar smoothing_factor,
    const GeometryTraits &gt)
{
  typedef boost::graph_traits< FaceGraph > GraphTraits;

  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits< PointMap >::value_type Point;
  typedef typename boost::property_traits< PointMap >::reference Reference;

  // this code works only for Polyhedron_3, Surface Mesh, OpenMesh and AIF
  auto iterator_pair = vertices(g); // vertices() returns a vertex_iterator pair

  // VertexMap barycenters_pm(get(boost::vertex_index, g));
  vertex_iterator vi = iterator_pair.first;
  vertex_iterator vi_end = iterator_pair.second;
  int cpt = 0;
  for(; vi != vi_end; ++vi)
  {
    cpt++;
    angles_based_centroids_pm[*vi] =
        Operators::vertex_one_ring_angles_based_centroid< FaceGraph,
                                                          PointMap,
                                                          GeometryTraits >(
            *vi, g, pm, smoothing_factor, gt);
  }
}

template< typename FaceGraph,
          typename PointMap,
          typename CentroidMap,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
void
calculate_vertices_one_ring_angles_based_centroid(
    const FaceGraph &g,
    const PointMap &pm, // const map (follows CGAL usage)
    CentroidMap
        angles_based_centroids_pm, // modifiable map (follows CGAL usage)
    const typename GeometryTraits::Scalar smoothing_factor = 0.1)

{
  GeometryTraits gt(g);
  return calculate_vertices_one_ring_angles_based_centroid< FaceGraph,
                                                           PointMap,
                                                           CentroidMap,
                                                           GeometryTraits >(
      g, pm, angles_based_centroids_pm, smoothing_factor, gt);
}

} // namespace Filters
} // namespace FEVV
