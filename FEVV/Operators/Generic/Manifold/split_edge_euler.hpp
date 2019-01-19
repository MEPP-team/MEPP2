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

#include <CGAL/boost/graph/internal/helpers.h> // set_border
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/boost/graph/Euler_operations.h> // for add_face(vr,g)

#include <CGAL/boost/graph/helpers.h> // is_border + is_triangle_mesh

#include <cassert>

#include "FEVV/Wrappings/Geometry_traits.h"

namespace FEVV {
namespace Operators {

/**
 * \brief Split an edge of the graph.
 *        The edge to split is given as a halfedge.
 *
 * \tparam  MutableFaceGraph a Mesh type that provides a Model of the
 *          MutableFaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \tparam  PointMap A modifiable point map to manage vertex positions.
 * \param g The MutableFaceGraph instance from which the edge will be
 *          split.
 * \param pm The mesh point map which associates vertex to positions.
 * \param h The edge to be split designated through one of its
 *          halfedges. 
 */
template< typename MutableFaceGraph, typename PointMap >
void
split_edge_euler(
    MutableFaceGraph &g,
    PointMap pm,
    typename boost::graph_traits< MutableFaceGraph >::halfedge_descriptor &h)
{
  typedef FEVV::Geometry_traits< MutableFaceGraph > GeometryTraits;
  typedef typename GeometryTraits::Point Point;

  typedef boost::graph_traits< MutableFaceGraph > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  GeometryTraits gt(g);

  if(h == boost::graph_traits< MutableFaceGraph >::null_halfedge())
    return;

  halfedge_descriptor opp_h = opposite(h, g);
  // ensure the edge is incident to triangles (only)
  if(CGAL::is_border(h, g))
  {
    if(!CGAL::is_triangle(opp_h, g))
      return;
  }
  else if(CGAL::is_border(opp_h, g))
  {
    if(!CGAL::is_triangle(h, g))
      return;
  }
  else
  {
    if(!CGAL::is_triangle(h, g) || !CGAL::is_triangle(opp_h, g))
      return;
  }

  // deal with vertex geometry
  vertex_descriptor vs = source(h, g);
  vertex_descriptor vt = target(h, g);

  // use CGAL::Euler function
  std::cout
      << "Warning: use CGAL::Euler::split_edge() and CGAL::Euler::split_face"
      << "\n";
  halfedge_descriptor res;
  if(CGAL::is_border(h, g))
  {
    res = CGAL::Euler::split_edge(opp_h, g); // does not work with AIFMesh
    CGAL::Euler::split_face(res, next(opp_h, g), g);
  }
  else if(CGAL::is_border(opp_h, g))
  {
    res = CGAL::Euler::split_edge(h, g);
    CGAL::Euler::split_face(res, next(h, g), g);
  }
  else
  {
    res = CGAL::Euler::split_edge(h, g);
    CGAL::Euler::split_face(res, next(h, g), g);
    CGAL::Euler::split_face(opposite(h, g), next(opposite(res, g), g), g);
  }
  put(pm,
      target(res, g),
      Point((gt.get_x(get(pm, vs)) + gt.get_x(get(pm, vt))) * 0.5f,
            (gt.get_y(get(pm, vs)) + gt.get_y(get(pm, vt))) * 0.5f,
            (gt.get_z(get(pm, vs)) + gt.get_z(get(pm, vt))) * 0.5f));
}

} // namespace Operators
} // namespace FEVV
