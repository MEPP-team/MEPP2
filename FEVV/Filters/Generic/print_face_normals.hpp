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
 * \brief  Print all faces normals.
 *
 * \param  g   mesh
 * \param  nm  normal map of the mesh
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
