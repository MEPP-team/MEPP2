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

#pragma push_macro("ERROR")
#undef ERROR

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Collapse_info.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Connectivity_encoding.h"

#include <list>
#include <vector>
#include <map>

namespace FEVV {
namespace Filters {

/** A Refinement_info Object stores the necessary info to refine a batch.
  * Each batch has its corresponding Refinement_info object.
  * Input: A sorted list of Collapse_info objects, a halfedgeGraph corresponding
  *        to the current level of detail.
  * Output: Corresponding bitmasks and residuals arrays.
  **/
template<
    typename HalfedgeGraph,
    typename PointMap,
    typename Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector,
    typename Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point,
    typename vertex_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
    typename vertex_iterator =
        typename boost::graph_traits< HalfedgeGraph >::vertex_iterator >
class Refinement_info
{ 
public:
  Refinement_info(
      HalfedgeGraph &g, 
      const std::list< Collapse_info< HalfedgeGraph, PointMap > > &list_memory)
      : _g(g), _list_memory(list_memory)
  {}

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int get_num_vertices() const
  {
    auto iterator_pair = vertices(_g);
    vertex_iterator vi = iterator_pair.first;
    vertex_iterator vi_end = iterator_pair.second;
    auto nb_vertices = std::distance(iterator_pair.first, iterator_pair.second);

    return static_cast< int >(nb_vertices);
  }

  /////////////////////////////////////////////////////////////////////////////
  /// Creates a vertex bitmask according to a vertex spanning tree and 
  /// a list of Collapse_info objects (initialized at construction).
  void set_bitMask(
                   const FEVV::Comparator::SpanningTreeVertexEdgeComparator< HalfedgeGraph,
                                             PointMap > &st /// spanning tree of current LOD
                   )
  {
	const std::list< vertex_descriptor > &spanning_tree = st.get_spanning_tree_vertices(); /// st vertices ordered in the same ordered
                                                                                           /// that the internal memory list.
    // Bit optimization stuff (to not encode predictable zeros).
	bool last_was_1 = false;	
    std::list< vertex_descriptor > remaining_adjacent_vertices_of_last_v, tmp;
    std::map<vertex_descriptor, bool> processed_vertices;  
	  
    typename std::list< vertex_descriptor >::const_iterator it =
        spanning_tree.begin(), it_e = spanning_tree.end();
    auto it_list = _list_memory.begin(), it_list_e = _list_memory.end();
    for( ; (it != it_e) && (it_list != it_list_e); ++it)
    {
      // Bit optimization stuff. 
      processed_vertices[*it] = true; 
	  bool is_adjacent_to_former_1 = false;
      if(last_was_1)
      {
        auto it_to_remove = std::find(remaining_adjacent_vertices_of_last_v.begin(), 
                                      remaining_adjacent_vertices_of_last_v.end(), *it) ;
        is_adjacent_to_former_1 = (it_to_remove != remaining_adjacent_vertices_of_last_v.end());
        if(is_adjacent_to_former_1)
          remaining_adjacent_vertices_of_last_v.erase(it_to_remove);
        if(remaining_adjacent_vertices_of_last_v.empty())
          last_was_1 = false;
      }
	  
      if(*it == (*it_list).get_vkept()) // Each split vertex must be present in 
                                        // the bit mask.
      { // Code by 1 vertex that belongs to the list.
        _bitMask.push_back(1);
		
        _other_info_bits.push_back(std::get< 0 >((*it_list).get_midpoint_rounds()));
        _other_info_bits.push_back(std::get< 1 >((*it_list).get_midpoint_rounds()));
        _other_info_bits.push_back(std::get< 2 >((*it_list).get_midpoint_rounds()));
        _other_info_bits.push_back(std::get< 3 >((*it_list).get_midpoint_rounds()));
		
        ++it_list;
		
        // Bit optimization stuff.
        last_was_1 = true; 
        tmp = FEVV::Comparator::get_not_processed_adjacent_vertices(*it, _g, processed_vertices, st.get_spanning_tree_min_incident_edge(*it));
        remaining_adjacent_vertices_of_last_v.insert(remaining_adjacent_vertices_of_last_v.end(), tmp.begin(), tmp.end());
      }
      else if(!is_adjacent_to_former_1)   // Each non-split vertex is present in 
                                          // the bit mask only when not adjacent 
                                          // to a previously met split vertex.
      {
        _bitMask.push_back(0);
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Creates the edge bitmask
  void set_connectivity_topology(
      const FEVV::Comparator::SpanningTreeVertexEdgeComparator< HalfedgeGraph,
                                                          PointMap > &st, /// spanning tree of current LOD
                                PointMap &pm, /// pointmap of current LOD
      std::list< Collapse_info< HalfedgeGraph, PointMap > > &list_memory /// sorted list of Collapse_info objects
	                                                                    /// (same list that _list_memory, but
																		/// _list_memory use the original sorted
																		/// list while here list_memory can have 
																		/// its object' reverse field updated)
  )
  {
    encode_connectivity_bitmask(_g, 
                              pm, 
                              list_memory, 
                              _connectivity, 
                              st 
                              );
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void set_error_prediction()
  {
    if(!_list_memory.empty())
    {
      auto it_list = _list_memory.begin(), ite = _list_memory.end();
      for( ; it_list != ite; ++it_list)
      {
        _error_prediction.push_back((*it_list).get_error_prediction());
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void set_reverse_bool()
  {
    if(!_list_memory.empty())
    {
      auto it_list = _list_memory.begin(), ite = _list_memory.end();
      for( ; it_list != ite; ++it_list)
        _reverse_bool.push_back((*it_list).get_reverse());
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  const std::list< bool >& get_bitmask() const { return _bitMask; }
  const std::list< bool >& get_connectivity() const { return _connectivity; }
  const std::list< bool >& get_other_info() const { return _other_info_bits; }
  const std::list< std::vector< Vector > >& get_error_prediction() const
  {
    return _error_prediction;
  }

  const std::list< bool >& get_reverse_bool() const { return _reverse_bool; }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

private:
  HalfedgeGraph &_g;
  const std::list< Collapse_info< HalfedgeGraph, PointMap > > &_list_memory; /// sorted list of refinement info
                                                                            /// (vertices to encode) according
                                                                            /// to the vertex spanning traversal																		
  std::list< bool > _bitMask; /// vertex bit mask (which vertex to split during)
                              /// decoding
  std::list< bool > _connectivity; /// edge bit mask (which incident edge to a 
                                   /// vertex to split is a pivot edge to expand)
  std::list< bool > _other_info_bits;
  std::list< std::vector< Vector > > _error_prediction;
  std::list< bool > _reverse_bool;
};


} // namespace Filters
} // namespace FEVV
