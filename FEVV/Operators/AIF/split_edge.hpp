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

//#include <CGAL/boost/graph/internal/helpers.h> // set_border
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/boost/graph/Euler_operations.h> // for add_face(vr,g)
//#include <CGAL/boost/graph/helpers.h> // is_border + is_triangle_mesh

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Operators/AIF/topology_predicates.hpp"

#include <cassert>

namespace FEVV {
namespace Operators {

/**
 * \brief Split an edge of the graph.
 *
 * \tparam  MutableFaceIncidentGraph a Mesh type that provides a Model of the
 *          MutableFaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param g The MutableFaceIncidentGraph instance from which the e edge will be
 *          split.
 * \param e The edge to be split.
 * \param midpoint_vertex The new vertex used for the splitting.
 * \return The newly inserted edge incident to inserted vertex (target vertex is the new one).
 */
template< typename MutableFaceIncidentGraph // similar concept to
                                             // MutableCellIncidentGraph, but
                                             // limited to 0, 1 and 2
                                             // dimensional cells
         >
typename boost::graph_traits<MutableFaceIncidentGraph>::edge_descriptor	
split_edge(
    MutableFaceIncidentGraph &g,
    typename boost::graph_traits< MutableFaceIncidentGraph >::edge_descriptor
        e, 
    typename boost::graph_traits<MutableFaceIncidentGraph>::vertex_descriptor 
	    midpoint_vertex)
{
  typedef FEVV::Geometry_traits< MutableFaceIncidentGraph > GeometryTraits;
  typedef typename GeometryTraits::Point Point;

  typedef boost::graph_traits< MutableFaceIncidentGraph > GraphTraits;
  typedef typename GraphTraits::edge_descriptor edge_descriptor;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  GeometryTraits gt(g);

  if(e == boost::graph_traits< MutableFaceIncidentGraph >::null_edge())
    return boost::graph_traits<MutableFaceIncidentGraph>::null_edge();

  // ensure the edge is incident to triangles (only)
  if(!FEVV::Operators::has_only_incident_triangular_faces(e, g))
  {
    return boost::graph_traits<MutableFaceIncidentGraph>::null_edge();
  }

  // deal with vertex geometry
  vertex_descriptor vs = source(e, g);
  vertex_descriptor vt = target(e, g);

  // when a topological modication is done over an edge, better to remove that
  // edge than to reuse it because we have to make sure the final user will
  // update edge properties accordingly.
  if (degree(e, g) == 0)
  {// one dangling edge is tranformed into 2 successive dangling edges
#ifdef USE_ADD_EDGE_AIF // conflict with function in Euler_operations.h (CGAL)
    std::pair<edge_descriptor, bool> p = add_edge(vs, midpoint_vertex, g);
    add_edge(midpoint_vertex, vt, g);
#else	
    edge_descriptor p = CGAL::Euler::add_edge(vs, midpoint_vertex, g);
    CGAL::Euler::add_edge(midpoint_vertex, vt, g);
#endif
    remove_edge(e, g); // remove at the end, to make sure no vertex is removed
#ifdef USE_ADD_EDGE_AIF
    return p.first;
#else
    return p;
#endif	
  }
  else
  {
    std::vector< vertex_descriptor > face1_vertices, face2_vertices;
    auto face_range_pair = in_edges(e, g); // get incident faces
    std::vector< face_descriptor > copy_f(
        face_range_pair.first,
        face_range_pair
            .second); // copy to not invalidate iterator while iterating
    typename std::vector< face_descriptor >::iterator iter_f(copy_f.begin()),
        iter_fe(copy_f.end());
    for(; iter_f != iter_fe; ++iter_f)
    {
      auto edge_range_pair = out_edges(*iter_f, g);
      auto iter_e = edge_range_pair.first;
      vertex_descriptor third_v;
      /////////////////////////////////////////////////////////////////////////
      while(*iter_e != e)
        ++iter_e;
      /////////////////////////////////////////////////////////////////////////
      edge_descriptor prev_e = *iter_e;
      ++iter_e;
      if(iter_e == edge_range_pair.second)
        iter_e = edge_range_pair.first;
      /////////////////////////////////////////////////////////////////////////
      if(source(prev_e, g) == source(*iter_e, g))
      {
        face1_vertices.push_back(target(*iter_e, g));
        face1_vertices.push_back(target(prev_e, g));
        face1_vertices.push_back(midpoint_vertex);

        face2_vertices.push_back(source(prev_e, g));
        face2_vertices.push_back(target(*iter_e, g));
        face2_vertices.push_back(midpoint_vertex);
      }
      else if(target(prev_e, g) == source(*iter_e, g))
      {
        face1_vertices.push_back(target(*iter_e, g));
        face1_vertices.push_back(source(prev_e, g));
        face1_vertices.push_back(midpoint_vertex);

        face2_vertices.push_back(target(prev_e, g));
        face2_vertices.push_back(target(*iter_e, g));
        face2_vertices.push_back(midpoint_vertex);
      }
      else if(source(prev_e, g) == target(*iter_e, g))
      {
        face1_vertices.push_back(source(*iter_e, g));
        face1_vertices.push_back(target(prev_e, g));
        face1_vertices.push_back(midpoint_vertex);

        face2_vertices.push_back(source(prev_e, g));
        face2_vertices.push_back(source(*iter_e, g));
        face2_vertices.push_back(midpoint_vertex);
      }
      else if(target(prev_e, g) == target(*iter_e, g))
      {
        face1_vertices.push_back(source(*iter_e, g));
        face1_vertices.push_back(source(prev_e, g));
        face1_vertices.push_back(midpoint_vertex);

        face2_vertices.push_back(target(prev_e, g));
        face2_vertices.push_back(source(*iter_e, g));
        face2_vertices.push_back(midpoint_vertex);
      }

      CGAL::Euler::add_face(face1_vertices, g);
      CGAL::Euler::add_face(face2_vertices, g);
      face1_vertices.clear();
      face2_vertices.clear();

      remove_face(*iter_f,
                  g); // remove at the end, to make sure no vertex or edge is
                      // removed while still needing by new geometries
    }
	
	  std::pair<edge_descriptor, bool> p2 = edge(vs, midpoint_vertex, g);
      return p2.first;
  }
}

/**
 * \brief Split an edge of the graph.
 *
 * \tparam  MutableFaceIncidentGraph a Mesh type that provides a Model of the
 *          MutableFaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param g The MutableFaceIncidentGraph instance from which the e edge will be
 *          split.
 * \param e The edge to be split.
 * \return The newly inserted edge incident to inserted vertex (target vertex is the new one).
 *         In this version, the new vertex is created internally via add_vertex(g).
 */
template< typename MutableFaceIncidentGraph // similar concept to MutableCellIncidentGraph, but limited to 0, 1 and 2 dimensional cells
               >
typename boost::graph_traits<MutableFaceIncidentGraph>::edge_descriptor  
	split_edge(MutableFaceIncidentGraph& g, typename boost::graph_traits<MutableFaceIncidentGraph>::edge_descriptor e)
{
  typedef boost::graph_traits<MutableFaceIncidentGraph>       GraphTraits;
  typedef typename GraphTraits::vertex_descriptor             vertex_descriptor;

  vertex_descriptor midpoint_vertex;
  midpoint_vertex = add_vertex(g);

  return split_edge<MutableFaceIncidentGraph>(g, e, midpoint_vertex);
}

/**
 * \brief Split an edge of the graph.
 *
 * \tparam  MutableFaceIncidentGraph a Mesh type that provides a Model of the
 *          MutableFaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \tparam  PointMap A modifiable point map to manage vertex positions.
 * \tparam  GeometryTraits  The geometric kernel when available. This is
 *                          defaulted to FEVV::Geometry_traits<MutableFaceIncidentGraph>. 
 * \param g The MutableFaceIncidentGraph instance from which the e edge will be
 *          split.
 * \param pm The mesh point map which associates vertex to positions.
 * \param e The edge to be split.
 * \return The newly inserted edge incident to inserted vertex (target vertex is the new one).
 *         In this version, the new vertex is created internally via add_vertex(g).
 */
template< typename MutableFaceIncidentGraph, // similar concept to MutableCellIncidentGraph, but limited to 0, 1 and 2 dimensional cells
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits<MutableFaceIncidentGraph> 
         >
typename boost::graph_traits<MutableFaceIncidentGraph>::edge_descriptor 	  
	split_edge(MutableFaceIncidentGraph& g, PointMap pm, typename boost::graph_traits<MutableFaceIncidentGraph>::edge_descriptor e)
{
  typedef typename GeometryTraits::Point                      Point;

  typedef boost::graph_traits<MutableFaceIncidentGraph>       GraphTraits;
  typedef typename GraphTraits::edge_descriptor               edge_descriptor;
  typedef typename GraphTraits::vertex_descriptor             vertex_descriptor;
  GeometryTraits gt(g);

  vertex_descriptor vs = source(e, g);
  vertex_descriptor vt = target(e, g);
  edge_descriptor res = split_edge<MutableFaceIncidentGraph>(g, e);

  if (res == boost::graph_traits<MutableFaceIncidentGraph>::null_edge())
    return boost::graph_traits<MutableFaceIncidentGraph>::null_edge();

  // deal with vertex geometry
  vertex_descriptor midpoint_vertex = target(res, g);

  put(pm, midpoint_vertex, Point( (gt.get_x(get(pm, vs)) + gt.get_x(get(pm, vt)))*0.5f,
                                  (gt.get_y(get(pm, vs)) + gt.get_y(get(pm, vt)))*0.5f,
                                  (gt.get_z(get(pm, vs)) + gt.get_z(get(pm, vt)))*0.5f));
  return res ;
}

} // namespace Operators
} // namespace FEVV
