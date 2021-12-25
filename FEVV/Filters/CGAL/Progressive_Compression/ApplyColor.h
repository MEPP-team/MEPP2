#pragma once
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/Euler_operations.h>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"
#include <iostream>
#include <vector>
#include <list>

namespace FEVV {
namespace Filters {
	
template<
	typename HalfedgeGraph,
	typename PointMap,
	typename EdgeColorMap,
	typename VertexColorMap,
	typename halfedge_descriptor =
	typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor >
void ColorStencil(HalfedgeGraph &g,
                  PointMap &pm,
                  EdgeColorMap &e_cm,
                  VertexColorMap &v_cm,
                  std::vector<halfedge_descriptor> &edges_to_color, 
                  typename boost::property_traits< VertexColorMap >::value_type color)
{
    size_t nbe = edges_to_color.size();
	for(size_t i = 0; i < nbe; i++)
	{
		put(e_cm, edge(edges_to_color[i], g), color); // colors every edge in color
	}
}

template<
	typename HalfedgeGraph,
	typename PointMap,
	typename EdgeColorMap,
	typename VertexColorMap,
	typename edge_descriptor =
	typename boost::graph_traits< HalfedgeGraph >::edge_descriptor >
void
ColorEdge(HalfedgeGraph &/*g*/,PointMap &/*pm*/,
          EdgeColorMap &e_cm,
          VertexColorMap &/*v_cm*/,
          edge_descriptor edge_to_color, 
          typename boost::property_traits< VertexColorMap >::value_type color)
{
	put(e_cm, edge_to_color, color);// colors every edge in color
}

template<
	typename HalfedgeGraph,
	typename PointMap,
	typename VertexColorMap,
	typename vertex_descriptor =
	typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor >
void
ColorCenter(HalfedgeGraph &/*g*/,
            PointMap &/*pm*/,
            VertexColorMap &v_cm,
            vertex_descriptor center,
            typename boost::property_traits<VertexColorMap>::value_type color )
{
	put(v_cm,center,color);
}
}
}
