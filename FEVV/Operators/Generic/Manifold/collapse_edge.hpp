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
#include "FEVV/Operators/Generic/Manifold/delete_edge_if_parallel.hpp"

namespace FEVV {
namespace Operators {

/**
 * \brief Collapse an edge of the graph.
 *        The edge to collapse is given as a halfedge.
 *        The halfedge source vertex is kept, and the halfedge
 *        target vertex is removed from the graph.
 *
 * \tparam  MutableFaceGraph a Mesh type that provides a Model of the
 *          MutableFaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param g The MutableFaceGraph instance out of which the h edge will be
 *          deleted.
 * \param h The edge to be deleted designated through one of its
 *          halfedges.
 */
template< typename MutableFaceGraph >
void
collapse_edge_keep_source(
    MutableFaceGraph &g,
    typename boost::graph_traits< MutableFaceGraph >::halfedge_descriptor &h)
{
  typedef boost::graph_traits< MutableFaceGraph > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef CGAL::Halfedge_around_target_iterator< MutableFaceGraph >
      halfedge_around_target_iterator;

  if(h == boost::graph_traits< MutableFaceGraph >::null_halfedge())
    return;

  // deal with vertex geometry
  vertex_descriptor vs = source(h, g);
  vertex_descriptor vt = target(h, g);

  halfedge_around_target_iterator hi, he;
  for(boost::tie(hi, he) = CGAL::halfedges_around_target(h, g); hi != he; ++hi)
    set_target(*hi, vs, g);
  remove_vertex(vt, g);

  // change edge connectivity
  halfedge_descriptor hn = next(h, g);
  halfedge_descriptor hp = prev(h, g);
  halfedge_descriptor ho = opposite(h, g);
  halfedge_descriptor hon = next(ho, g);
  halfedge_descriptor hop = prev(ho, g);
  set_next(hp, hn, g);
  set_next(hop, hon, g);
  remove_edge(edge(h, g), g);

  // remove parallel edges
  FEVV::Operators::delete_edge_if_parallel(g, hn);
  FEVV::Operators::delete_edge_if_parallel(g, hop);
}

/**
 * \brief Collapse an edge of the graph.
 *        The edge to collapse is given as a halfedge.
 *        The halfedge target vertex is kept, and the halfedge
 *        source vertex is removed from the graph.
 *
 * \tparam  MutableFaceGraph a Mesh type that provides a Model of the
 *          MutableFaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param g The MutableFaceGraph instance out of which the h edge will be
 *          deleted.
 * \param h The edge to be deleted designated through one of its
 *          halfedges.
 */
template< typename MutableFaceGraph >
void
collapse_edge_keep_target(
    MutableFaceGraph &g,
    typename boost::graph_traits< MutableFaceGraph >::halfedge_descriptor &h)
{
  typedef boost::graph_traits< MutableFaceGraph > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor ho = opposite(h, g);
  collapse_edge_keep_source(g, ho);
  // TODO
  // calling collapse_edge_keep_source(g, opposite(h, g))
  // doesn't compile. WHY ?
}

} // namespace Operators
} // namespace FEVV

