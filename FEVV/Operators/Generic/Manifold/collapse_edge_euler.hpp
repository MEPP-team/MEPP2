#ifndef COLLAPSE_EDGE_EULER_H
#define COLLAPSE_EDGE_EULER_H

#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>

namespace FEVV {
namespace Operators {

/**
 * \brief Collapse an edge of the graph.
 *        The edge to collapse is given as a halfedge.
 *        The halfedge source vertex is kept, and the halfedge
 *        target vertex is removed from the graph.
 *
 * Template parameters:
 *           MutableFaceGraph: a Mesh type that provides
 *                 a Model of the MutableFaceGraph Concept
 *                 through a boost::graph_traits<> specialization.
 */
template< typename MutableFaceGraph >
void
collapse_edge_keep_source_euler(
    MutableFaceGraph &g,
    typename boost::graph_traits< MutableFaceGraph >::halfedge_descriptor &h)
{
  // use CGAL::Euler function
  std::cout << "Warning: using CGAL::Euler::collapse_edge()"
            << "\n";
  // CGAL's documentation doesn't specify which one of the
  // source and the target vertex is removed.
  // We experimentaly found that the target vertex is removed,
  // so we need to use "opposite()" here to remove the source
  // vertex.
  CGAL::Euler::collapse_edge(edge(opposite(h, g), g), g);
}

/**
 * \brief Collapse an edge of the graph.
 *        The edge to collapse is given as a halfedge.
 *        The halfedge target vertex is kept, and the halfedge
 *        source vertex is removed from the graph.
 *
 * Template parameters:
 *           MutableFaceGraph: a Mesh type that provides
 *                 a Model of the MutableFaceGraph Concept
 *                 through a boost::graph_traits<> specialization.
 */
template< typename MutableFaceGraph >
void
collapse_edge_keep_target_euler(
    MutableFaceGraph &g,
    typename boost::graph_traits< MutableFaceGraph >::halfedge_descriptor &h)
{
  typedef boost::graph_traits< MutableFaceGraph > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor ho = opposite(h, g);
  collapse_edge_keep_source_euler(g, ho);
  // todo
  // calling collapse_edge_keep_source(g, opposite(h, g))
  // doesn't compile. WHY ?
}

} // namespace Operators
} // namespace FEVV

#endif // COLLAPSE_EDGE_EULER_H
