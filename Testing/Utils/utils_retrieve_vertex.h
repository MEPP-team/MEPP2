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
