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

namespace FEVV {
namespace Operators {

/*!
 * \brief  Function determining if the arguments are similar edges (parallel
 *         edges).
 *
 * \tparam  FaceIncidentGraph a Mesh type that provides a Model of the
 *          FaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param	  edge1	The first involving edge
 * \param	  edge2	The second involving edge
 * \param   g The FaceIncidentGraph instance from which the edges are taken.
 * \return	true if the arguments have the same vertices (not necessarily in
 *          the same position.)
 */
template< typename FaceIncidentGraph >
static bool
are_similar_edges(
    const typename boost::graph_traits< FaceIncidentGraph >::edge_descriptor
        edge1,
    const typename boost::graph_traits< FaceIncidentGraph >::edge_descriptor
        edge2,
    const FaceIncidentGraph &g)
{
  if(edge1 == edge2)
    return true;
  return (source(edge1, g) == source(edge2, g) &&
          target(edge1, g) == target(edge2, g)) ||
         (source(edge1, g) == target(edge2, g) &&
          target(edge1, g) == source(edge2, g));
}

/*!
 * \brief   Function determining if a pair of similar/parallel edges exists
 *          among the incident edges of vertex v_to_keep.
 *
 * \tparam  FaceIncidentGraph a Mesh type that provides a Model of the
 *          FaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param	  v_to_keep	The involving vertex
 * \param   g The FaceIncidentGraph instance from which the vertex is taken.
 * \return	true if the vertex has at least one pair of incident similar
 *          edges.
 */
template< typename FaceIncidentGraph >
static bool
contains_similar_incident_edges(
    const typename boost::graph_traits< FaceIncidentGraph >::vertex_descriptor
        v_to_keep,
    const FaceIncidentGraph &g)
{
  auto edges_range_pair = in_edges(v_to_keep, g);
  auto iter_e = edges_range_pair.first;
  for(; iter_e != edges_range_pair.second; ++iter_e)
  {
    auto iter_e_bis = iter_e;
    ++iter_e_bis;
    for(; iter_e_bis != edges_range_pair.second; ++iter_e_bis)
    {
      if(Operators::are_similar_edges(*iter_e, *iter_e_bis, g))
        return true;
    }
  }
  return false;
}

/*!
 * \brief Function determining if the arguments are similar faces (parallel
 *        faces).
 *
 * \tparam  FaceIncidentGraph a Mesh type that provides a Model of the
 *          FaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param	  f1	The first involving face
 * \param	  f2	The second involving face
 * \param   g The FaceIncidentGraph instance from which the faces are taken.
 * \param   take_into_account_face_rientation A boolean to tell if the similar
 *          faces must have the same orientation.
 * \return	true if the arguments are similar faces given
 *          take_into_account_face_rientation.
 */
template< typename FaceIncidentGraph >
static bool
are_similar_faces(
    const typename boost::graph_traits< FaceIncidentGraph >::face_descriptor f1,
    const typename boost::graph_traits< FaceIncidentGraph >::face_descriptor f2,
    const FaceIncidentGraph &g,
    bool take_into_account_face_rientation =
        false /// false means two faces sharing the same vertices not with the
              /// same orientation are similar, else they must have the same
              /// orientation
)
{
  if(degree(f1, g) == degree(f2, g))
  {
    auto edges_range_pair = out_edges(f1, g); // incident edges of f1
    auto iter_e = edges_range_pair.first;
    for(; iter_e != edges_range_pair.second; ++iter_e)
    { // for each incident edge of the face f1
      bool is_similar = false;
      auto edges_range_pair2 = out_edges(f2, g); // incident edges of f2
      auto iter_e2 = edges_range_pair2.first;
      for(; iter_e2 != edges_range_pair2.second; ++iter_e2)
      { // for each incident edge of the face f2
        if(Operators::are_similar_edges(*iter_e, *iter_e2, g))
        {
          if(take_into_account_face_rientation)
          { // the first time *iterE and *iterE2 are similar edges, we need to
            // continue the comparison up to come back
            auto iter_ee = iter_e;
            do
            {
              ++iter_e;
              if(iter_e == edges_range_pair.second)
                iter_e = edges_range_pair.first;
              ++iter_e2;
              if(iter_e2 == edges_range_pair2.second)
                iter_e2 = edges_range_pair2.first;

              if(!Operators::are_similar_edges(
                     *iter_e, *iter_e2, g))
                return false;

            } while(iter_e != iter_ee);

            return true;
          }
          else
            is_similar = true;
          break;
        }
      }
      if(!is_similar)
        return false;
    }
    return true;
  }

  return false;
}

/*!
 * \brief Function determining if the incident faces to e are all triangles.
 *
 * \tparam  IncidentMutableFaceGraph a Mesh type that provides a Model of the
 *          IncidentMutableFaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param	  e	The involving edge.
 * \param   g The IncidentMutableFaceGraph instance from which the edge e is
 *          taken.
 * \return  true if all incident faces to e are triangles.
 */
template< typename IncidentMutableFaceGraph >
inline bool
has_only_incident_triangular_faces(
    const typename boost::graph_traits<
        IncidentMutableFaceGraph >::edge_descriptor e,
    const IncidentMutableFaceGraph &g)
{
  auto face_range_pair = in_edges(e, g);
  auto iter = face_range_pair.first;
  for(; iter != face_range_pair.second; ++iter)
  {
    if(degree(*iter, g) != 3)
      return false;
  }
  return true;
}

} // namespace Operators
} // namespace FEVV

