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
 * \brief Flip an edge of the graph.
 *        The edge to flip is given as a halfedge.
 *
 * \tparam  MutableFaceGraph a Mesh type that provides a Model of the
 *          MutableFaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param g The MutableFaceGraph instance out of which the h edge will be
 *          flipped.
 * \param h The edge to be flipped designated through one of its
 *          halfedges.
 */
template< typename MutableFaceGraph >
void
flip_edge(
    MutableFaceGraph &g,
    typename boost::graph_traits< MutableFaceGraph >::halfedge_descriptor &h)
{
  typedef FEVV::Geometry_traits< MutableFaceGraph > GeometryTraits;

  typedef boost::graph_traits< MutableFaceGraph > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  GeometryTraits gt(g);

  if(h == boost::graph_traits< MutableFaceGraph >::null_halfedge())
    return;

  halfedge_descriptor opp_h = opposite(h, g);
  // ensure the edge is incident to 2 triangles (only)
  if(CGAL::is_border(h, g))
  {
    return;
  }
  else if(CGAL::is_border(opp_h, g))
  {
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

  typedef CGAL::Halfedge_around_face_iterator< MutableFaceGraph >
      Halfedge_around_face_iterator;
  vertex_descriptor midpoint_vertex;

  std::vector< vertex_descriptor > face_vertices;
  // here we assume a 2_manifold context
  vertex_descriptor opp_v1 = target(next(h, g), g);
  vertex_descriptor opp_v2 = target(next(opp_h, g), g);

  halfedge_descriptor h1 = next(h, g), h2 = prev(h, g);
  halfedge_descriptor h3 = next(opp_h, g), h4 = prev(opp_h, g);
  remove_face(face(h, g), g);
  CGAL::internal::set_border(h1, g);
  CGAL::internal::set_border(h2, g);

  remove_face(face(opp_h, g), g);
  CGAL::internal::set_border(h3, g);
  CGAL::internal::set_border(h4, g);

  face_vertices.push_back(vs);
  face_vertices.push_back(opp_v2);
  face_vertices.push_back(opp_v1);
  face_descriptor f1 = CGAL::Euler::add_face(face_vertices, g);
  face_vertices.clear();
  face_vertices.push_back(opp_v1);
  face_vertices.push_back(opp_v2);
  face_vertices.push_back(vt);
  face_descriptor f2 = CGAL::Euler::add_face(face_vertices, g);
  ///////////////////////////////////////////////////////////////////////////
  Halfedge_around_face_iterator hi, he;
  int cpt = 0;
  for(boost::tie(hi, he) = CGAL::halfedges_around_face(halfedge(f1, g), g);
      hi != he;
      ++hi)
  {
    if(CGAL::is_border_edge(*hi, g))
    {
      // std::cout << "border!" << std::endl;
      CGAL::internal::set_border(opposite(*hi, g), g);
    }
  }
  for(boost::tie(hi, he) = CGAL::halfedges_around_face(halfedge(f2, g), g);
      hi != he;
      ++hi)
  {
    if(CGAL::is_border_edge(*hi, g))
    {
      // std::cout << "border!" << std::endl;
      CGAL::internal::set_border(opposite(*hi, g), g);
    }
  }
}

} // namespace Operators
} // namespace FEVV
