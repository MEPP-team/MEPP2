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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include "FEVV/Tools/Comparator/Spanning_tree_vertex_edge_comparator.hpp"

#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Memory_comparator.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/apply_color.h"

namespace FEVV {
namespace Filters {

/** 
  * \brief Encode vertex to split' one-ring edges to expand. 
  * For a vertex to split (represented by a Collapse_info object), will encode  
  * the edges to expand in a binary stream, 1 being an edge that should be 
  * expanded and 0 being an edge that should not. Takes as arguments a graph g,
  * a pointmap, a Collapse_info object (again, representing a vertex and many
  * necessary refinement info), the first edge to process (see the 
  * encode_connectivity_bitmask function), and the edge/connectivity bitmask 
  * (which will receive the connectivity bitstream).
**/
template<
    typename HalfedgeGraph,
    typename PointMap,
    typename halfedge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor,
    typename edge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::edge_descriptor,
    typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
void
connectivity_encoding(const HalfedgeGraph &g, /// in
                      const PointMap &/*pm*/, /// in
                      Collapse_info< HalfedgeGraph, PointMap > &mem, /// in/out
                      edge_descriptor starting_edge, /// in min edge (edge with the lowest index in the spanning tree)
                      std::list< bool > &edge_bitmask /// out
                     )
{
  // convert edge into halfedge
  halfedge_descriptor starting_halfedge = halfedge(starting_edge, g);
  // the current edge starts as the "starting edge" and will turn around its
  // target using the topology of the mesh
  halfedge_descriptor current_edge = starting_halfedge;

  bool first = true;
  int count = 0;
  do
  {
    if(count > 1)
		break; // not needed to add the last false
    if(source(current_edge, g) == mem.get_v3() || source(current_edge, g) == mem.get_v4())
    {
      // if v3 goes after v4, we have to reverse our residuals in order to
      // ensure the proper decompression of our mesh
      if(source(current_edge, g) == mem.get_v4() && first == true)
      {
        mem.record_reverse(true);
      }
      edge_bitmask.push_back(true);
      count++;
      first = false;
    }
    else
    {
      edge_bitmask.push_back(false);
    }
	
    current_edge = opposite(next(current_edge, g), g);
  } while(current_edge != starting_halfedge);
  if(mem.get_v4() != boost::graph_traits< HalfedgeGraph >::null_vertex() &&
     mem.get_v3() == boost::graph_traits< HalfedgeGraph >::null_vertex())
  {
    mem.record_reverse(false);
  }
  if(mem.get_v3() != boost::graph_traits< HalfedgeGraph >::null_vertex() &&
     mem.get_v4() == boost::graph_traits< HalfedgeGraph >::null_vertex())
  {
    mem.record_reverse(true);
  }
}

/// \brief Encode vertex to split' one-ring edges to expand, for 
///        all collapsed edges in the order of the sorted list list_memory.
template<
    typename HalfedgeGraph,
    typename PointMap,
    typename halfedge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor,
    typename edge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::edge_descriptor,
    typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
void
encode_connectivity_bitmask(
    const HalfedgeGraph &g,
    const PointMap &pm,
    std::list< Collapse_info< HalfedgeGraph, PointMap > > &list_memory, /// sorted list (according to st)
                                                                        /// but non-const reference
                                                                        /// due to Collapse_info updates
    std::list< bool > &edge_bitmask, /// out
    const FEVV::Comparator::Spanning_tree_vertex_edge_comparator< HalfedgeGraph,
                                                        PointMap,
                                                        Geometry >& spanningtree /// in computed st (needed
                                                                                 /// to find the min edges)
    )
{
  auto mem_it = list_memory.begin(), mem_it_e = list_memory.end();
  for( ; mem_it != mem_it_e; ++mem_it)
  {
    auto min_edge = spanningtree.get_spanning_tree_min_incident_edge((*mem_it).get_vkept());
    connectivity_encoding(g,
                         pm,
                         *mem_it,
                         min_edge,
                         edge_bitmask
                         );
  }
}

} // namespace Filters
} // namespace FEVV
