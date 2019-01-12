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

namespace FEVV {
namespace Operators {

/**
 * \brief   Delete the edge provided as argument when it is parallel to
 *          another edge of the graph (this is an if and only if condition)
 *          i.e. when the edge (as opposed to half edges) provided as argument
 *          and (at least) another edge are adjacent to the same two vertices.
 * \tparam  MutableFaceGraph a Mesh type that provides a Model of the
 *          MutableFaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param g The MutableFaceGraph instance out of which the h edge will be
 *          deleted.
 * \param h The edge to be deleted designated through one of its
 *          halfedges.
 * \post    When present, the face adjacent to the argument halfedge is
 *          deleted form the graph too
 * \note    The argument is an halfedge as opposed to an edge because
 *          it enables to designate the (possible) adjacent face that will
 *          be removed.
 */
template< typename MutableFaceGraph >
void
delete_edge_if_parallel(
    MutableFaceGraph &g,
    typename boost::graph_traits< MutableFaceGraph >::halfedge_descriptor &h)
{
  typedef boost::graph_traits< MutableFaceGraph > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  if(h == boost::graph_traits< MutableFaceGraph >::null_halfedge())
    return;

  halfedge_descriptor hn = next(h, g);
  if(next(hn, g) != h)
    return; // no parallel edge, nothing to do

  // change edges connections
  halfedge_descriptor ho = opposite(h, g);
  halfedge_descriptor hon = next(ho, g);
  halfedge_descriptor hop = prev(ho, g);
  set_next(hop, hn, g);
  set_next(hn, hon, g);

  // change adjacent face
  face_descriptor fc = face(ho, g);
  if(fc != boost::graph_traits< MutableFaceGraph >::null_face())
  {
    set_face(hn, fc, g);
    set_halfedge(fc, hn, g);
  }

  // remove ancient adjacent face, if any
  face_descriptor fd = face(h, g);
  if(fd != boost::graph_traits< MutableFaceGraph >::null_face())
    remove_face(fd, g);

  // remove edge
  remove_edge(edge(h, g), g);
}

} // namespace Operators
} // namespace FEVV
