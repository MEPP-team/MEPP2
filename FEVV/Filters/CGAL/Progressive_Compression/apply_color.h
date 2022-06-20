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
#include <CGAL/boost/graph/Euler_operations.h>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"

#include <vector>

namespace FEVV {
namespace Filters {
/**
 *  \brief Colors an edge by putting the specified color into
 *         the edge color map at the "edge" position.
 *         Used for debug purpose.
 */
template<
	typename HalfedgeGraph,
	typename PointMap,
	typename EdgeColorMap,
	typename edge_descriptor =
	typename boost::graph_traits< HalfedgeGraph >::edge_descriptor >
void
color_edge(HalfedgeGraph &/*g*/,PointMap &/*pm*/,
           EdgeColorMap &e_cm,
           edge_descriptor edge_to_color, 
           typename boost::property_traits< EdgeColorMap >::value_type color)
{
	put(e_cm, edge_to_color, color);
}

/**
 *  \brief Colors a vertex by putting the specified color into
 *         the vertex color map at the "vertex" position.
 *         Used for debug purpose. 
 */
template<
	typename HalfedgeGraph,
	typename PointMap,
	typename VertexColorMap,
	typename vertex_descriptor =
	typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor >
void
color_vertex(HalfedgeGraph &/*g*/,
             PointMap &/*pm*/,
             VertexColorMap &v_cm,
             vertex_descriptor vertex_to_color,
             typename boost::property_traits<VertexColorMap>::value_type color)
{
	put(v_cm, vertex_to_color, color);
}

} // namespace Filters
} // namespace FEVV
