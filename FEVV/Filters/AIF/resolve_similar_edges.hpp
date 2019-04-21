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

#include <boost/graph/graph_traits.hpp>
#include "FEVV/Operators/AIF/topology_predicates.hpp"

namespace FEVV {
namespace Filters {

/**
 * \brief Remove/resolve similar edges for the given mesh g.
 *
 * \tparam  MutableFaceIncidentGraph a Mesh type that provides a Model of the
 *          MutableFaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param   g The MutableFaceIncidentGraph instance.
 */
template< typename MutableFaceIncidentGraph >
static void
resolve_similar_edges(MutableFaceIncidentGraph &g)
{
  typedef
      typename boost::graph_traits< MutableFaceIncidentGraph >::edge_descriptor
          edge_descriptor;
  std::vector< edge_descriptor > edges_to_remove;
  /////////////////////////////////////////////////////////////////////////////
  auto edges_range_pair = edges(g);
  auto iter_e = edges_range_pair.first;
  for(; iter_e != edges_range_pair.second; ++iter_e)
  {
    if(std::find(edges_to_remove.begin(), edges_to_remove.end(), *iter_e) !=
       edges_to_remove.end())
      continue;
    auto iter_e_bis = iter_e;
    ++iter_e_bis;
    for(; iter_e_bis != edges_range_pair.second; ++iter_e_bis)
    {
      if(Operators::are_similar_edges(*iter_e, *iter_e_bis, g))
      {
        edges_to_remove.push_back(*iter_e_bis);
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  typename std::vector< edge_descriptor >::iterator it_e(
      edges_to_remove.begin()),
      it_ee(edges_to_remove.end());
  for(; it_e != it_ee; ++it_e)
    remove_edge(*it_e, g); // REDUNDANT EDGES ARE REMOVED
  edges_to_remove.clear();
}

} // namespace Filters
} // namespace FEVV
