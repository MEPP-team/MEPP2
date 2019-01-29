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
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/boost/graph/Euler_operations.h> // for add_face(vr,g)

#include "FEVV/Operators/AIF/similarity.hpp"
#include "FEVV/DataStructures/AIF/AIFMesh.hpp"

#include <algorithm> // std::swap
#include <cassert>

namespace FEVV {
namespace Operators {


typedef FEVV::DataStructures::AIF::AIFTopologyHelpers AIFHelpers;

/**
 * \brief Collapse an edge of the graph.
 *        The edge to collapse is given as a non-oriented edge.
 *        The edge source vertex is kept, and the edge
 *        target vertex is removed from the graph.
 *
 * \tparam  MutableFaceIncidentGraph a Mesh type that provides a Model of the
 *          MutableFaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param g The MutableFaceIncidentGraph instance from which the e edge will be
 *          deleted.
 * \param e The edge to be deleted.
 */
template< typename MutableFaceIncidentGraph > 
void
collapse_edge_keep_source(MutableFaceIncidentGraph &g,
                          typename boost::graph_traits<
                                                MutableFaceIncidentGraph >::
                                                edge_descriptor &e)
{
  typedef boost::graph_traits< MutableFaceIncidentGraph > GraphTraits;
  typedef typename GraphTraits::edge_descriptor edge_descriptor;
  /////////////////////////////////////////////////////////////////////////////
  if (e == GraphTraits::null_edge())
    return;

  auto pair_e_bool = edge(source(e, g), target(e, g), g);
  if(!pair_e_bool.second)
    return;
  /////////////////////////////////////////////////////////////////////////////
  edge_descriptor retrieved_edge = pair_e_bool.first;

  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  /////////////////////////////////////////////////////////////////////////////
  //                REMOVE INCIDENT TRIANGULAR FACES
  /////////////////////////////////////////////////////////////////////////////
  // 1 remove retrieved_edge from all its incident faces
  std::vector< face_descriptor > faces_to_remove;
  auto faces_range_pair = in_edges(e, g); // get incident faces
  auto iter_f = faces_range_pair.first;
  for(; iter_f != faces_range_pair.second; ++iter_f)
  {
    remove_out_edge(
        *iter_f, retrieved_edge); // remove edge retrieved_edge from face *iterF
    if(degree(*iter_f, g) < 3)
      faces_to_remove.push_back(*iter_f);
  }
  /////////////////////////////////////////////////////////////////////////////
  // 2 remove all faces that must be removed from mesh
  typename std::vector< face_descriptor >::iterator it(faces_to_remove.begin()),
      ite(faces_to_remove.end());
  for(; it != ite; ++it)
    remove_face(*it, g);
  faces_to_remove.clear();
  /////////////////////////////////////////////////////////////////////////////
  //                REMOVE THE EDGE AND MERGE ITS TWO VERTICES
  /////////////////////////////////////////////////////////////////////////////
  // 3 remove retrieved_edge
  vertex_descriptor t = target(retrieved_edge, g);
  vertex_descriptor s = source(retrieved_edge, g);

  remove_edge(retrieved_edge, g);
  /////////////////////////////////////////////////////////////////////////////
  // cannot handle multiple incident similar edges (assumed by Florian Caillaud
  // => easier encoding for mesh compression)
  FEVV::Operators::resolve_similar_incident_edges(s, g);
  FEVV::Operators::resolve_similar_incident_edges(t, g);
  // cannot handle multiple incident similar faces (assumed by Florian Caillaud
  // => easier encoding for mesh compression) Else a particular attention would
  // be needed to distinguish between similar faces added intentionaly (or read
  // from file) against those obtained via the edge collapse that must be
  // resolved
  FEVV::Operators::resolve_similar_incident_faces(s, g, false);
  FEVV::Operators::resolve_similar_incident_faces(t, g, false);
  /////////////////////////////////////////////////////////////////////////////
  // 4 merge retrieved_edge' vertices (source and target): the source is kept!

  // all incident edges of target vertex are put onto source vertex
  auto src_incident_edges_range = in_edges(t, g);
  auto iter_e = src_incident_edges_range.first;
  for(; iter_e != src_incident_edges_range.second; ++iter_e)
  {
    AIFHelpers::add_vertex(*iter_e, s, AIFHelpers::vertex_position(*iter_e, t));
  }
  /////////////////////////////////////////////////////////////////////////////
  // merge the 2 vertices
  // for each incident edge of the target vertex t
  //    if that edge was not incident to source vertex s, add it to s
  //     else add current edge to a list of edges to be removed from mesh
  //          replace this edge to be removed by the similar edge incident to s
  // remove all edges in the list of edges to be removed
  std::vector< edge_descriptor > edges_to_remove;
  src_incident_edges_range = in_edges(t, g);
  iter_e = src_incident_edges_range.first;
  for(; iter_e != src_incident_edges_range.second; ++iter_e)
  {
    edge_descriptor edge = AIFHelpers::common_edge(
        (*iter_e)->get_first_vertex(), (*iter_e)->get_second_vertex());
    if(edge == AIFHelpers::null_edge())
      AIFHelpers::add_edge(s, *iter_e);
    else
    {
      edges_to_remove.push_back(*iter_e);
      FEVV::Operators::merge_similar_edges(edge, *iter_e, g);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  typename std::vector< edge_descriptor >::iterator it_e(
      edges_to_remove.begin()),
      it_ee(edges_to_remove.end());
  for(; it_e != it_ee; ++it_e)
  {
    // AIFHelpers::add_vertex(*itE, t, AIFHelpers::vertex_position(*itE, s)); //
    // to unsure that t is correctly updated during the next removal
    remove_edge(*it_e, g); // REDUNDANT EDGES ARE REMOVED (can also remove t
                           // vertex if above ligne is uncommented)
  }
  edges_to_remove.clear();
  /////////////////////////////////////////////////////////////////////////////
  AIFHelpers::clear_vertex(t);
  /////////////////////////////////////////////////////////////////////////////
  // 5 unsure that the incident faces of s does not contain t vertex as incident
  // vertex
  //   then remove similar faces
  FEVV::Operators::resolve_similar_incident_faces(s, g, false);
  /////////////////////////////////////////////////////////////////////////////
  // facesRange = AIFHelpers::incident_faces(s);
  // std::cout << "nb incident faces = " << facesRange.size() << std::endl;
  assert(!FEVV::Operators::contains_similar_incident_edges(s, g));
  /////////////////////////////////////////////////////////////////////////////
  // 6 remove t vertex
  remove_vertex(t, g);
}


/**
 * \brief Collapse an edge of the graph.
 *        The edge to collapse is given as a non-oriented edge.
 *        The edge target vertex is kept, and the edge
 *        source vertex is removed from the graph.
 *
 * \tparam  MutableFaceIncidentGraph a Mesh type that provides a Model of the
 *          MutableFaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param g The MutableFaceIncidentGraph instance from which the e edge will be
 *          deleted.
 * \param e The edge to be deleted.
 */
template< typename MutableFaceIncidentGraph >
void collapse_edge_keep_target(MutableFaceIncidentGraph &g,
                               typename boost::graph_traits<
                                                MutableFaceIncidentGraph >::
                                                edge_descriptor &e)
{
  typedef boost::graph_traits< MutableFaceIncidentGraph > GraphTraits;
  typedef typename GraphTraits::edge_descriptor edge_descriptor;

  auto pair_e_bool = edge(target(e, g), source(e, g), g);
  if(pair_e_bool.second)
  {
    edge_descriptor retrieved_edge = pair_e_bool.first;
    collapse_edge_keep_source(g, retrieved_edge);
  }
}

} // namespace Operators
} // namespace FEVV
