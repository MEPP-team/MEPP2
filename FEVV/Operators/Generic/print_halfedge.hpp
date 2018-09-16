#ifndef PRINT_HALFEDGES_H
#define PRINT_HALFEDGES_H

#include <iostream>
#include <boost/graph/graph_traits.hpp>

namespace FEVV {
namespace Operators {


/**
 * \brief  Print the source and the target coordinates of one halfedge.
 *
 * \param  g   mesh
 * \param  h   the halfedge
 * \param  pm  point map of the mesh
 *
 * \ingroup  GenericManifoldFilters GenericNonManifoldFilters
 */
template< typename HalfedgeGraph, typename PointMap >
void
print_halfedge(
    const typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor &h,
    PointMap pm,
    const HalfedgeGraph &g)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;

  vertex_descriptor vs = source(h, g);
  vertex_descriptor vt = target(h, g);
  std::cout << "Halfedge: source = " << pm[vs] << " target = " << pm[vt]
            << "\n";
}


} // namespace Operators
} // namespace FEVV

#endif // PRINT_HALFEDGES_H
