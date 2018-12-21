#pragma once

#include <iostream>
#include <boost/graph/graph_traits.hpp>

namespace FEVV {
namespace Filters {


/**
 * \brief  Print all vertices coordinates.
 *
 * \param  g   mesh
 * \param  pm  point map of the mesh
 *
 * \ingroup  GenericManifoldFilters GenericNonManifoldFilters
 */
template< typename HalfedgeGraph, typename PointMap >
void
print_points(const HalfedgeGraph &g, PointMap pm)
{
  std::cout << "Points:" << std::endl;
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  for(vertex_iterator it = vertices(g).begin(); it != vertices(g).end(); ++it)
  {
    std::cout << pm[*it] << std::endl;
  }
}


} // namespace Filters
} // namespace FEVV
