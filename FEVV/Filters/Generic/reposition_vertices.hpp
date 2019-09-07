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
#include <boost/graph/properties.hpp>

#include "FEVV/Wrappings/Geometry_traits.h"

namespace FEVV {
namespace Filters {
/**
 * \brief Apply new vertex positions to vertex position property map.
 *
 * \tparam PropertyGraph  a Mesh type that provides a Model of the
 *                        PropertyGraph Concept through a boost::graph_traits<>
 *                        specialization.
 * \tparam PointMap       A modifiable point map to manage vertex positions.
 * \tparam NewPointMap    A modifiable point map to manage vertex' position
 *                        proposals.
 * \tparam GeometryTraits  The geometric kernel when available. This is
 *                         defaulted to FEVV::Geometry_traits<PropertyGraph>.
 * \param  g              The PropertyGraph instance from which the vertices
 *                        will be taken.
 * \param  pm             The point map which associate a vertex descriptor
 *                        with a vertex position.
 * \param  smoothed_pm    The point map used to replace the vertex positions
 *                        in pm.
 * \param  gt             The geometry trait object.
 */
template< typename PropertyGraph,
          typename PointMap,
          typename NewPointMap,
          typename GeometryTraits = FEVV::Geometry_traits< PropertyGraph > >
void
reposition_vertices(
    const PropertyGraph &g,
    PointMap pm,                    // modifiable map (follows CGAL usage)
    const NewPointMap &smoothed_pm, // const map (follows CGAL usage)
    const GeometryTraits &gt)
{
  typedef boost::graph_traits< PropertyGraph > GraphTraits;

  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename boost::property_traits< PointMap >::value_type Point;
  typedef typename boost::property_traits< PointMap >::reference Reference;

  // this code works only for Polyhedron_3, Surface Mesh, OpenMesh and AIF
  auto iterator_pair = vertices(g); // vertices() returns a vertex_iterator pair
  // auto smoothed_pm = get_barycenter_pm<PropertyGraph, PointMap,
  // GeometryTraits >(g, pm, gt);
  vertex_iterator vi = iterator_pair.first;
  vertex_iterator vi_end = iterator_pair.second;

  for(; vi != vi_end; ++vi)
  {
    Point p = pm[*vi];
    put(pm, *vi, smoothed_pm[*vi]);
  }
}

/*
 * \brief The only purpose of this function is to instantiate the last argument
 *        in order to simplify caller's syntax. This is to circumvent an
 *        apparent impossibility of providing a default argument out of a
 *        type drilled down from a default template parameter i.e. we cannot
 *        define something like
 *          template<,, typename Geometry=FEVV::Geometry_traits<PropertyGraph> >
 *          apply_Smoothing(,, const Geometry& gt = Geomtry( g ) )
 */
template< typename PropertyGraph,
          typename PointMap,
          typename NewPointMap,
          typename GeometryTraits = FEVV::Geometry_traits< PropertyGraph > >
void
reposition_vertices(
    const PropertyGraph &g,
    PointMap pm,                   // modifiable map (follows CGAL usage)
    const NewPointMap &smoothed_pm // const map (follows CGAL usage)
)
{
  GeometryTraits gt(g);
  reposition_vertices< PropertyGraph, PointMap, NewPointMap, GeometryTraits >(
      g, pm, smoothed_pm, gt);
}

} // namespace Filters
} // namespace FEVV
