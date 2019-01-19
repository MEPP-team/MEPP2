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
 * \brief   Split an edge of the graph.
 *          The edge to split is given as an halfedge.

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
split_edge(
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

  typedef CGAL::Halfedge_around_face_iterator< MutableFaceGraph >
      Halfedge_around_face_iterator;
  vertex_descriptor midpoint_vertex;
  midpoint_vertex = add_vertex(g);
  put(pm,
      midpoint_vertex,
      Point((gt.get_x(get(pm, vs)) + gt.get_x(get(pm, vt))) * 0.5f,
            (gt.get_y(get(pm, vs)) + gt.get_y(get(pm, vt))) * 0.5f,
            (gt.get_z(get(pm, vs)) + gt.get_z(get(pm, vt))) * 0.5f));

  std::vector< vertex_descriptor > face_vertices;
  // here we assume a 2_manifold context, we thus have 2 cases: border or no
  // border
  if(CGAL::is_border(h, g))
  {
    // opp_h is not a border
    vertex_descriptor opp_v = target(next(opp_h, g), g);

    halfedge_descriptor h1 = next(opp_h, g), h2 = prev(opp_h, g);
    remove_face(face(opp_h, g), g);
    // if (halfedge(edge(opp_h, g),g) !=
    // boost::graph_traits<MutableFaceGraph>::null_halfedge())
    //  remove_edge(edge(opp_h, g), g);
    CGAL::internal::set_border(h1, g);
    CGAL::internal::set_border(h2, g);

    face_vertices.push_back(opp_v);
    face_vertices.push_back(vt);
    face_vertices.push_back(midpoint_vertex);
    CGAL::Euler::add_face(face_vertices, g);

    face_vertices.clear();
    face_vertices.push_back(opp_v);
    face_vertices.push_back(midpoint_vertex);
    face_vertices.push_back(vs);
    CGAL::Euler::add_face(face_vertices, g);
  }
  else if(CGAL::is_border(opp_h, g))
  {
    // h is not a border
    vertex_descriptor opp_v = target(next(h, g), g);

    halfedge_descriptor h1 = next(h, g), h2 = prev(h, g);
    remove_face(face(h, g), g);
    // if (halfedge(edge(h, g),g) !=
    // boost::graph_traits<MutableFaceGraph>::null_halfedge())
    //  remove_edge(edge(h, g), g);
    CGAL::internal::set_border(h1, g);
    CGAL::internal::set_border(h2, g);

    face_vertices.push_back(opp_v);
    face_vertices.push_back(vs);
    face_vertices.push_back(midpoint_vertex);
    CGAL::Euler::add_face(face_vertices, g);
    face_vertices.clear();
    face_vertices.push_back(opp_v);
    face_vertices.push_back(midpoint_vertex);
    face_vertices.push_back(vt);
    CGAL::Euler::add_face(face_vertices, g);
  }
  else
  { // no border
    vertex_descriptor opp_v1 = target(next(h, g), g);
    vertex_descriptor opp_v2 = target(next(opp_h, g), g);
    halfedge_descriptor h1 = next(h, g), h2 = prev(h, g);
    halfedge_descriptor h3 = next(opp_h, g), h4 = prev(opp_h, g);
    remove_face(face(h, g), g);
    CGAL::internal::set_border(h1, g);
    CGAL::internal::set_border(h2,
                               g); // CGAL::internal::set_border(prev(h1,g), g);

    remove_face(face(opp_h, g), g);
    CGAL::internal::set_border(h3, g);
    CGAL::internal::set_border(h4, g);

    face_vertices.push_back(vs);
    face_vertices.push_back(midpoint_vertex);
    face_vertices.push_back(
        opp_v1); // prefer new point in the middle due to CGAL::add_face
    face_descriptor f1 = CGAL::Euler::add_face(face_vertices, g);
    face_vertices.clear();
    face_vertices.push_back(opp_v1);
    face_vertices.push_back(midpoint_vertex);
    face_vertices.push_back(vt);
    face_descriptor f2 = CGAL::Euler::add_face(face_vertices, g);
    ///////////////////////////////////////////////////////////////////////////
    Halfedge_around_face_iterator hi, he;
    // int cpt=0;
    for(boost::tie(hi, he) = CGAL::halfedges_around_face(halfedge(f1, g), g);
        hi != he;
        ++hi)
    {
      if(CGAL::is_border_edge(*hi, g))
      {
        // std::cout << "border!" << std::endl;
        CGAL::internal::set_border(opposite(*hi, g), g);
        auto h_if_exist = halfedge(vs, opp_v2, g);
        if(h_if_exist.second && (target(*hi, g) == midpoint_vertex) &&
           (next(opposite(*hi, g), g) != h_if_exist.first))
          set_next(opposite(*hi, g), h_if_exist.first, g);
        //++cpt;
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
        auto h_if_exist = halfedge(opp_v2, vt, g);
        if(h_if_exist.second && (source(*hi, g) == midpoint_vertex) &&
           (next(h_if_exist.first, g) != opposite(*hi, g)))
          set_next(h_if_exist.first, opposite(*hi, g), g);
        //++cpt;
      }
    }
    // assert(cpt==2);
    ///////////////////////////////////////////////////////////////////////////
    face_vertices.clear();
    face_vertices.push_back(vt);
    face_vertices.push_back(midpoint_vertex);
    face_vertices.push_back(
        opp_v2); // prefer new point in the middle due to CGAL::add_face
    f1 = CGAL::Euler::add_face(face_vertices, g);
    for(boost::tie(hi, he) = CGAL::halfedges_around_face(halfedge(f1, g), g);
        hi != he;
        ++hi)
    {
      if(CGAL::is_border_edge(*hi, g))
      {
        // std::cout << "border!" << std::endl;
        CGAL::internal::set_border(opposite(*hi, g), g);
        //++cpt;
      }
    }
    // assert(cpt==3);

    face_vertices.clear();
    face_vertices.push_back(opp_v2);
    face_vertices.push_back(midpoint_vertex);
    face_vertices.push_back(vs);
    f1 = CGAL::Euler::add_face(face_vertices, g);
    for(boost::tie(hi, he) = CGAL::halfedges_around_face(halfedge(f1, g), g);
        hi != he;
        ++hi)
    {
      if(CGAL::is_border_edge(*hi, g))
      {
        std::cout << "border!" << std::endl;
        //++cpt;
      }
    }
    // assert(cpt==3);
  }
}

} // namespace Operators
} // namespace FEVV
