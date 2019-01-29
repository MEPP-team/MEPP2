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
#include <boost/foreach.hpp>
#include <limits>

namespace FEVV {
namespace Filters {

/**
 * \brief  Generic call to compute min max from a Descriptor
 *		  Get the value of the property map for the given descriptor and
 *update minMetric and maxMetric accordingly
 *
 * \param  d                 Descriptor used as key for the property map
 * \param  prop_map			Property map containing numerical values
 * \param  minMetric         minimum value of the property map
 * \param  maxMetric         maximum value of the property map
 *
 */
template< typename Descriptor,
          typename PropertyMap,
          typename MapType = typename PropertyMap::value_type >
void
compute_min_max_descriptor(const Descriptor &d,
                           const PropertyMap &prop_map,
                           MapType &min_metric,
                           MapType &max_metric)
{
  MapType l_metric = get(prop_map, d);
  if(min_metric > l_metric)
    min_metric = l_metric;
  if(max_metric < l_metric)
    max_metric = l_metric;
}


/**
 * \brief  Compute min and max value of a numerical property map for the faces
 * of a mesh
 *
 * \param  g                 mesh with the property map
 * \param  prop_map			Property map containing numerical values
 * \param  minMetric         minimum value of the property map
 * \param  maxMetric         maximum value of the property map
 *
 */
template< typename HalfedgeGraph,
          typename PropertyMap,
          typename MapType = typename PropertyMap::value_type >
void
compute_min_max_faces(const HalfedgeGraph &g,
                      const PropertyMap &prop_map,
                      MapType &min_metric,
                      MapType &max_metric)
{
  typedef typename boost::graph_traits< HalfedgeGraph >::face_descriptor
      Face_Descriptor;
  min_metric = (std::numeric_limits< MapType >::max)();
  max_metric = (std::numeric_limits< MapType >::min)();
  BOOST_FOREACH(Face_Descriptor f, faces(g))
  {
    compute_min_max_descriptor(f, prop_map, min_metric, max_metric);
  }
}


/**
 * \brief  Compute min and max value of a numerical property map for the
 * vertices of a mesh
 *
 * \param  g                 mesh with the property map
 * \param  prop_map			Property map containing numerical values
 * \param  minMetric         minimum value of the property map
 * \param  maxMetric         maximum value of the property map
 *
 */
template< typename HalfedgeGraph,
          typename PropertyMap,
          typename MapType = typename PropertyMap::value_type >
void
compute_min_max_vertices(const HalfedgeGraph &g,
                         const PropertyMap &prop_map,
                         MapType &min_metric,
                         MapType &max_metric)
{
  typedef typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor
      Vertex_Descriptor;

  min_metric = (std::numeric_limits< MapType >::max)();
  max_metric = (std::numeric_limits< MapType >::min)();
  BOOST_FOREACH(Vertex_Descriptor v, vertices(g))
  {
    compute_min_max_descriptor(v, prop_map, min_metric, max_metric);
  }
}


/**
 * \brief  Compute min and max value of a numerical property map for the
 * halfedges of a mesh
 *
 * \param  g                 mesh with the property map
 * \param  prop_map			Property map containing numerical values
 * \param  minMetric         minimum value of the property map
 * \param  maxMetric         maximum value of the property map
 *
 */
template< typename HalfedgeGraph,
          typename PropertyMap,
          typename MapType = typename PropertyMap::value_type >
void
compute_min_max_halfedges(const HalfedgeGraph &g,
                          const PropertyMap &prop_map,
                          MapType &min_metric,
                          MapType &max_metric)
{
  typedef typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor
      Halfedge_Descriptor;
  min_metric = (std::numeric_limits< MapType >::max)();
  max_metric = (std::numeric_limits< MapType >::min)();

  BOOST_FOREACH(Halfedge_Descriptor h, halfedges(g))
  {
    compute_min_max_descriptor(h, prop_map, min_metric, max_metric);
  }
}
} // namespace Filters
} // namespace FEVV
