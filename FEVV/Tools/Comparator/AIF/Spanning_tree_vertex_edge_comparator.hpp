// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
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

#include "FEVV/Tools/Comparator/Spanning_tree_vertex_edge_comparator.hpp"

#include <limits>
#include <iterator>
#include <map>
#include <set>

#include <stack>

namespace FEVV { 
namespace Comparator 
{
  template<>
  std::list< typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::vertex_descriptor>
	  get_adjacent_vertices(typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::vertex_descriptor v,
		  const FEVV::DataStructures::AIF::AIFMesh& /*g*/)
  {
	  typedef FEVV::DataStructures::AIF::AIFTopologyHelpers Helpers;
	  auto res = Helpers::adjacent_vertices(v);
	  return std::list< typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::vertex_descriptor>(res.begin(), res.end());
  }

  template<>
  std::list< typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::vertex_descriptor>
    get_not_processed_adjacent_vertices(typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::vertex_descriptor v,
      const FEVV::DataStructures::AIF::AIFMesh& g,
      std::map<typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::vertex_descriptor, bool>& processed_vertices,
      typename boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh>::edge_descriptor min_e)
  {
    typedef boost::graph_traits<FEVV::DataStructures::AIF::AIFMesh> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor                 vertex_descriptor;

    std::list<vertex_descriptor> adjacent_vertices;

    if (min_e != GraphTraits::null_edge())
    { // works in 2-manifold (thus to avoid when locally non-manifold)
      typedef typename GraphTraits::halfedge_descriptor             halfedge_descriptor;

      // find the (second) edge that separates old region and new region
      halfedge_descriptor h = halfedge(min_e, g);
      if (target(h, g) != v)
        h = opposite(h, g);

      size_t cpt = 0; // needed because sometimes the condition !(h == h_end) may
                   // rise an issue when h_end value is never met again (due to
                   // fact that we may have a locally non-manifold situation)
      halfedge_descriptor h_end = h;
      do
      {
        vertex_descriptor v_adj = source(h, g);
        auto it_res = processed_vertices.find(v_adj);
        if ((it_res == processed_vertices.end()) || !it_res->second)
        { // not yet processed
          if (it_res != processed_vertices.end())
            it_res->second = true;
          else
            processed_vertices[v_adj] = true;
          adjacent_vertices.push_back(source(h, g));
        }

        h = opposite(next(h, g), g);
        ++cpt;
      } while (!(h == h_end) && (cpt<= degree(v, g)));
    }
    else
    { // works in 2-manifold and non-manifold
      typedef FEVV::DataStructures::AIF::AIFTopologyHelpers           Helpers;

      auto res = Helpers::adjacent_vertices(v);

      auto iter_v = res.begin(), iter_v_e = res.end();
      for (; iter_v != iter_v_e; ++iter_v) {
        auto it_res = processed_vertices.find(*iter_v);
        if ((it_res != processed_vertices.end()) && it_res->second)
          continue;
        else {
          if (it_res != processed_vertices.end())
            it_res->second = true;
          else
            processed_vertices[*iter_v] = true;
          adjacent_vertices.push_back(*iter_v);
        }
      }
    }
	  return adjacent_vertices;
  }

} // namespace Comparator
} // namespace FEVV

