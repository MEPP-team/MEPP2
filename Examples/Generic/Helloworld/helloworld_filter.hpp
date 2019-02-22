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

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"

/**
 * \ingroup GenericManifoldFilters
 * \brief Refer \ref ExampleComprehensive "here" for a detailed presentation
 *        of that filter example, what it does and how to use it.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename VertexNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
helloworld_filter(const HalfedgeGraph &g,
                   PointMap &pm,
                   VertexColorMap &v_cm,
                   VertexNormalMap &v_nm,
                   const GeometryTraits &gt)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  typedef typename GraphTraits::face_iterator face_iterator;
  typedef typename boost::property_traits< PointMap >::value_type Point;
  typedef typename boost::property_traits< VertexColorMap >::value_type Color;
  typedef typename boost::property_traits< VertexNormalMap >::value_type Normal;

  // notes:
  //  - a property map stores a single property value for each
  //    vertex/halfedge/face it is related to ;
  //  - the property map content can only be accessed with a
  //    vertex/halfedge/face identifier (called a 'descriptor') ;
  //  - in order to display the whole property map content one must loop
  //    over all vertices/halfedges/faces

  // loop over vertices
  int counter = 0;
  auto iterator_pair = vertices(g);
  // vertices() returns a vertex_iterator pair ;
  // the same for faces(g) and edges(g)
  vertex_iterator vi = iterator_pair.first;
  vertex_iterator vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    counter++;
    std::cout << "processing vertex #" << counter << std::endl;

    // read property map data
    Point p = get(pm, *vi);    // gt.get_x(p), gt.get_y(p), gt.get_z(p)
    Color c = get(v_cm, *vi);  // c[0], c[1], c[2]
    Normal n = get(v_nm, *vi); // n[0], n[1], n[2]

    // write property map data
    Point newpoint(gt.get_x(p) + 1, gt.get_y(p) + 2, gt.get_z(p) + 1);
    put(pm, *vi, newpoint);
    Color newcolor(0.8, 0.2, 0.1);
    put(v_cm, *vi, newcolor);
    Normal newnormal(1.0, 2.0, 3.0);
    put(v_nm, *vi, newnormal);
  }

  // create and populate a filter-specific property map related to faces

  // create a property map
  typename FEVV::Face_pmap< HalfedgeGraph, float > my_property_map;
  my_property_map = FEVV::make_face_property_map< HalfedgeGraph, float >(g);

  // populate the property map
  // loop over the faces
  auto face_iterator_pair = faces(g);
  face_iterator fi = face_iterator_pair.first;
  face_iterator fi_end = face_iterator_pair.second;
  for(; fi != fi_end; ++fi)
  {
    static float value = 42.1f;
    put(my_property_map, *fi, value);
    value += 1.0f;
  }

  // display the property map
  // loop over the faces
  counter = 0;
  fi = face_iterator_pair.first;
  for(; fi != fi_end; ++fi)
  {
    counter++;
    float value = get(my_property_map, *fi);
    std::cout << "my_property_map[face#" << counter << "] = " << value
              << std::endl;
  }
}

// Helper function to simplify (syntactic sugar) the call to helloworld_filter
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename VertexNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
helloworld_filter(const HalfedgeGraph &g,
                   PointMap &pm,
                   VertexColorMap &v_cm,
                   VertexNormalMap &v_nm)

{
  GeometryTraits gt(g);
  helloworld_filter(g, pm, v_cm, v_nm, gt);
}
