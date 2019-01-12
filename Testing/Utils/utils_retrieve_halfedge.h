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

#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>
#include <boost/next_prior.hpp> // for boost::next(T x, Distance n)
#include <boost/graph/graph_traits.hpp>

/*
 * \brief Utility functions that retrieves the halfedge out
 *        of its source and destination indexes.
 * \tparam MutableFaceGraph Any mesh type having a boost::graph_traits
 *                    specialisation.
 * \param g           The graph which which to search for.
 * \param sourceIndex The index of the source vertex.
 * \param targetIndex The index of the target vertex.
 * \return            When found, the halfedge going from source to target.
 */
template< typename MutableFaceGraph >
typename boost::graph_traits< MutableFaceGraph >::halfedge_descriptor
retrieve_halfedge(MutableFaceGraph &g, int source_index, int target_index)
{
  typedef boost::graph_traits< MutableFaceGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;

  vertex_iterator v =
      vertices(g).first; // g.vertices_begin() cannot be used here, because it
                         // is not defined in any concept

  vertex_descriptor src = *boost::next(v, source_index);
  vertex_descriptor dst = *boost::next(v, target_index);

  halfedge_descriptor result;
  bool found = false;
  boost::tie(result, found) = halfedge(src, dst, g);
  if(!found)
  {
    return GraphTraits::null_halfedge();
  }
  return result;
}
