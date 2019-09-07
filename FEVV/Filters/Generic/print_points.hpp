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

#include <iostream>
#include <boost/graph/graph_traits.hpp>

namespace FEVV {
namespace Filters {


/**
 * \brief  Print all vertices coordinates.
 *
 * \param  g   mesh
 * \param  pm  point map of the mesh
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
