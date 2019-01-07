#pragma once

#include <iostream>
#include <boost/graph/graph_traits.hpp>

namespace FEVV {
namespace Filters {


/**
 * \brief  Print all faces normals.
 *
 * \param  g   mesh
 * \param  nm  normal map of the mesh
 *
 * \ingroup  GenericManifoldFilters GenericNonManifoldFilters
 */
template< typename HalfedgeGraph, typename NormalMap >
void
print_face_normals(const HalfedgeGraph &g, NormalMap nm)
{
  std::cout << "Face normals:" << std::endl;
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::face_iterator face_iterator;
  for(face_iterator it = faces(g).begin(); it != faces(g).end(); ++it)
  {
    // Facet_iterator is a face_descriptor, so we can use it as the
    // key here
    std::cout << nm[*it] << std::endl;
  }
}


} // namespace Filters
} // namespace FEVV

