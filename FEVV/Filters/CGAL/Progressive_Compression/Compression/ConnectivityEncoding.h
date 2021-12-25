#pragma once

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include "FEVV/Tools/Comparator/SpanningTreeVertexEdgeComparator.hpp"

#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/MemoryComparator.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/ApplyColor.h"

namespace FEVV {
namespace Filters {

/**
  * For a vertex (represented by a CollapseInfo object), will encode the edges 
  * to expand in a binary stream, 1 being an edge that should be expanded and 0
  * being an edge that should not. Takes as arguments a graph g, a pointmap,
  * [color maps for edges and vertices], a CollapseInfo object (again, 
  * representing a vertex and many necessary refinement info), the first edge 
  * to process (see the EncodeConnectivityBitmask function), the 
  * edge/connectivity bitmask (which will receive the connectivity bitstream), 
  * and a neighbours bitmask: used for debugging purposes, in order to separate
  * every connectivity bitmask (each vertex will have its connectivity written 
  * separately) 
**/
template<
    typename HalfedgeGraph,
    typename PointMap,
    //typename EdgeColorMap,
    //typename VertexColorMap,
    typename vertex_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
    typename halfedge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor,
    typename edge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::edge_descriptor,
    typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
void
ConnectivityEncoding(const HalfedgeGraph &g,
                     const PointMap &/*pm*/,
                     //EdgeColorMap &/*ecm*/,
                     //VertexColorMap &/*vcm*/,
                     CollapseInfo< HalfedgeGraph, PointMap > &mem,
                     edge_descriptor starting_edge, /// min edge (edge with the lowest index in the spanning tree)
                     std::list< bool > &edge_bitmask
                     //, std::vector< std::vector< bool > > &neighbours_bitmask /// for debug; \todo remove it (costly)
                     )
{
  // convert edge into halfedge
  halfedge_descriptor starting_halfedge = halfedge(starting_edge, g);
  // the current edge starts as the "starting edge" and will turn around its
  // target using the topology of the mesh
  halfedge_descriptor current_edge = starting_halfedge;

  // for *debug purposes*, color the edges in a known order in order to see any
  // issue in the path taken to go around the target vertex
  //typedef typename boost::property_traits< VertexColorMap >::value_type Color;
  //std::vector< Color > colors = {Color(1, 0, 0),
  //                               Color(0, 1, 1),
  //                               Color(1, 1, 0),
  //                               Color(0, 0, 1),
  //                               Color(0.4, 0.4, 0.4),
  //                               Color(1, 1, 1),
  //                               Color(0.9, 0.5, 0),
  //                               Color(0.4, 0, 0.4)};
  //int i = 0;

  //std::vector< bool > current;
  bool first = true;
  int count = 0;
  do
  {
    if( count > 1)
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
      //current.push_back(true);
      count++;
      first = false;
    }
    else
    {
      edge_bitmask.push_back(false);
      //current.push_back(false);
    }

    // ColorEdge(g,pm,ecm,vcm,edge(current_edge,g),colors[i%8]);
    //i++;
	
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
  //neighbours_bitmask.push_back(current);
}

template<
    typename HalfedgeGraph,
    typename PointMap,
    //typename EdgeColorMap,
    //typename VertexColorMap,
    typename vertex_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
    typename halfedge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor,
    typename edge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::edge_descriptor,
    typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >
void
EncodeConnectivityBitmask(
    const HalfedgeGraph &g,
    const PointMap &pm,
    //EdgeColorMap &ecm,
    //VertexColorMap &vcm,
    std::list< CollapseInfo< HalfedgeGraph, PointMap > > &list_memory, /// sorted list (according to st)
	                                                                   /// but non-const reference
																	   /// due to CollapseInfo updates
    std::list< bool > &edge_bitmask,
    const FEVV::Comparator::SpanningTreeVertexEdgeComparator< HalfedgeGraph,
                                                        PointMap,
                                                        Geometry >& spanningtree /// computed st
    //, std::vector< std::vector< bool > > &neighbours_bitmask /// for debug
    )
{
  auto mem_it = list_memory.begin(), mem_it_e = list_memory.end();
  //typedef typename boost::property_traits< VertexColorMap >::value_type Color;
  for(; mem_it != mem_it_e; ++mem_it)
  {
    auto min_edge = spanningtree.get_spanning_tree_min_incident_edge((*mem_it).get_vkept());
    ConnectivityEncoding(g,
                         pm,
                         //ecm,
                         //vcm,
                         *mem_it,
                         min_edge,
                         edge_bitmask
                         //, neighbours_bitmask
                         );
  }
}

} // namespace Filters
} // namespace FEVV
