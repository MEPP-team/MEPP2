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

/*
 * \brief Utility functions that retrieves the vertex out
 *        of its index.
 * \tparam MutableFaceGraph Any mesh type having a boost::graph_traits
 *                    specialisation.
 * \param g           The graph which which to search for.
 * \param target_index The index of the target vertex.
 * \return            When found, the vertex with target index.
 */   
template<typename MutableFaceGraph>
typename boost::graph_traits<MutableFaceGraph>::vertex_descriptor
retrieve_vertex( MutableFaceGraph& g,
                 int target_index )
{
  typedef boost::graph_traits<MutableFaceGraph>     GraphTraits;
  typedef typename GraphTraits::vertex_iterator     vertex_iterator;
  typedef typename GraphTraits::vertex_descriptor   vertex_descriptor;

  vertex_iterator v = vertices(g).first; // g.vertices_begin() cannot be used here, 
                                         // because it is not defined in any concept

  vertex_descriptor dst = *boost::next(v, target_index);

  return dst;
}
