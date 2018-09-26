#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>

/*
 * \brief Utility functions that retrieves the vertex out
 *        of its index.
 * \tparam MutableFaceGraph Any mesh type having a boost::graph_traits
 *                    specialisation.
 * \param g           The graph which which to search for.
 * \param targetIndex The index of the target vertex.
 * \return            When found, the vertex with target index.
 */   
template<typename MutableFaceGraph>
typename boost::graph_traits<MutableFaceGraph>::face_descriptor
retrieve_face( MutableFaceGraph& g,
               int targetIndex )
{
  typedef boost::graph_traits<MutableFaceGraph>     GraphTraits;
  typedef typename GraphTraits::face_iterator       face_iterator;
  typedef typename GraphTraits::face_descriptor     face_descriptor;

  face_iterator f = faces(g).first; // g.vertices_begin() cannot be used here, because it is not defined in any concept

  face_descriptor dst = *boost::next(f, targetIndex);

  return dst;
}
