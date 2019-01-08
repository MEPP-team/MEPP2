#pragma once

#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>

/*
 * \brief Utility functions that retrieves the face out
 *        of its index.
 * \tparam MutableFaceGraph Any mesh type having a boost::graph_traits
 *                    specialisation.
 * \param g           The graph which which to search for.
 * \param target_index The index of the target face.
 * \return            When found, the face with target index.
 */   
template<typename MutableFaceGraph>
typename boost::graph_traits<MutableFaceGraph>::face_descriptor
retrieve_face( MutableFaceGraph& g,
               int target_index )
{
  typedef boost::graph_traits<MutableFaceGraph>     GraphTraits;
  typedef typename GraphTraits::face_iterator       face_iterator;
  typedef typename GraphTraits::face_descriptor     face_descriptor;

  face_iterator f = faces(g).first;

  face_descriptor dst = *boost::next(f, target_index);

  return dst;
}
