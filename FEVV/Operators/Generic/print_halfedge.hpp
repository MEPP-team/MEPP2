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
namespace Operators {


/**
 * \brief  Print the source and the target coordinates of one halfedge.
 *
 * \param  g   mesh
 * \param  h   the halfedge
 * \param  pm  point map of the mesh
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

