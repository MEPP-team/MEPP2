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
namespace Operators {

/**
 * \brief   Geometric fold detection, usually during iterations of an
 *          CGAL::Halfedge_around_target_circulator object cir_he.
 *
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \tparam  PointMap A modifiable point map to manage vertex positions.
 * \tparam  GeometryTraits The geometric kernel when available. This is
 *          defaulted to FEVV::Geometry_traits<FaceGraph>.
 * \param   h1 The first involved halfedge. May be
 *          prev(prev(opposite(*cir_he, g), g), g).
 * \param   h2 The second involved halfedge. May be next(*cir_he, g).
 * \param   g The FaceGraph instance from which the halfedges are taken.
 * \param   pm The mesh point map which associates vertex to positions.
 * \param   gt The geometry trait object.
 * \see     vertex_one_ring_angles_based_centroid filter.
 */
template< typename FaceGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
inline bool
is_geometrical_fold(
    const typename boost::graph_traits< FaceGraph >::halfedge_descriptor h1,
    const typename boost::graph_traits< FaceGraph >::halfedge_descriptor h2,
    const FaceGraph &g,
    const PointMap &pm,
    const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  // Variables declaration
  Vector next_h1, nnext, left_dir = gt.NULL_VECTOR, vg, next_h2,
                         right_dir = gt.NULL_VECTOR, vh;
  //////////////////////////////////////////////////////////////////
  next_h1 =
      gt.sub(get(pm, target(h1, g)),
             get(pm,
                 target(next(h1, g),
                        g))); // target(next(h1, g),g) instead of
                              // source(prev(h1, g) to handle polygonal cases

  // Compute unit normal
  nnext = FEVV::Operators::
      calculate_vertex_normal< FaceGraph, PointMap, GeometryTraits >(
          target(next(h1, g), g), g, pm, gt);

  if(next_h1 != gt.NULL_VECTOR)
    left_dir = gt.cross_product(next_h1, nnext); // LeftDir = next_h1 X Nnext
  else
    return true;
  //////////////////////////////////////////////////////////////////
  next_h2 = gt.sub(get(pm, target(h2, g)), get(pm, source(prev(h2, g), g)));

  // Compute unit normal
  nnext = FEVV::Operators::
      calculate_vertex_normal< FaceGraph, PointMap, GeometryTraits >(
          source(prev(h2, g), g), g, pm, gt);

  if(next_h2 != gt.NULL_VECTOR)
    right_dir = gt.cross_product(next_h2, nnext); // RightDir = next_h2 X Nnext
  else
    return true;

  vg = gt.sub(get(pm, target(h2, g)), get(pm, target(next(h1, g), g)));
  vh = gt.sub(get(pm, target(h1, g)), get(pm, source(prev(h2, g), g)));

  // fold to the right     fold to the left
  return ((gt.dot_product(left_dir, vg) <= 0.0) ||
          (gt.dot_product(right_dir, vh) <= 0.0));
}

} // namespace Operators
} // namespace FEVV
