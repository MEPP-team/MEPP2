// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include "FEVV/Wrappings/properties.h"
#include "FEVV/Wrappings/Wrappings_cgal_point_set.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_point_set.h"

#include <string>


namespace FEVV {

// _PMap_traits specialization for FEVV::CGALPointSet
// use CGAL::Point_set_3 internal property maps
// (which are the same as CGAL::Surface_mesh internal
// property maps)


// define default property map type (associative property map,
// vector property map...) for each cell type (vertex, face...) ;

// a counter used to generate a unique property map identifier
static size_t cgal_point_set_prop_map_counter = 0;

// default vertex property map
template< typename ValueT >
struct Vertex_pmap_traits< FEVV::CGALPointSet, ValueT >
{
  // use CGAL::Point_set_3 internal property map
  typedef typename FEVV::CGALPointSet::Property_map< ValueT > pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    std::string name = std::to_string(cgal_point_set_prop_map_counter);
    cgal_point_set_prop_map_counter++;
    return const_cast< FEVV::CGALPointSet & >(m)
        .add_property_map< ValueT >(name)
        .first;
  }
};

#if 0 // useless for point cloud

// default face property map
template< typename ValueT >
struct Face_pmap_traits< FEVV::CGALPointSet, ValueT >
{
  typedef typename boost::property_map< FEVV::CGALPointSet,
                                        boost::face_index_t >::const_type
      face_index_map_type;
  typedef typename boost::vector_property_map< ValueT, face_index_map_type >
      pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default edge property map
template< typename ValueT >
struct Edge_pmap_traits< FEVV::CGALPointSet, ValueT >
{
  typedef typename boost::property_map< FEVV::CGALPointSet,
                                        boost::edge_index_t >::const_type
      edge_index_map_type;
  typedef typename boost::vector_property_map< ValueT, edge_index_map_type >
      pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default halfedge property map
template< typename ValueT >
struct Halfedge_pmap_traits< FEVV::CGALPointSet, ValueT >
{
  typedef typename boost::property_map< FEVV::CGALPointSet,
                                        boost::halfedge_index_t >::const_type
      halfedge_index_map_type;
  typedef typename boost::vector_property_map< ValueT, halfedge_index_map_type >
      pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

#endif

// define standard property map types for (property,cell) pair,
// for example vertex-normal

// specialize the property maps traits for vertex-normal
template<>
struct _PMap_traits< FEVV::CGALPointSet, FEVV::vertex_normal_t >
{
  typedef FEVV::CGALPointSet::Vector_map   pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    return const_cast< FEVV::CGALPointSet& >(m).add_normal_map().first;
  }
};

#if 0
// specialize the property maps traits for vertex-tangent
template<>
struct _PMap_traits< FEVV::CGALPointSet, FEVV::vertex_tangent_t >
{
  typedef FEVV::CGALPointSet::Vector_map   pmap_type;
  typedef FEVV::CGALPointSet::Vector_3     Vector_3;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    return const_cast< FEVV::CGALPointSet & >(m)
        .add_property_map< Vector_3 >("vertex_tangent")
        .first;
  }
};

// specialize the property maps traits for vertex-texcoord
template<>
struct _PMap_traits< FEVV::CGALPointSet, FEVV::vertex_texcoord_t >
{
  typedef FEVV::CGALPointSet::Vector_map   pmap_type;
  typedef FEVV::CGALPointSet::Vector_3     Vector_3;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    return const_cast< FEVV::CGALPointSet & >(m)
        .add_property_map< Vector_3 >("vertex_texcoord")
        .first;
  }
};
#endif

// specialize the standard property map for vertex-color
template<>
struct _PMap_traits< FEVV::CGALPointSet, FEVV::vertex_color_t >
{
  typedef FEVV::CGALPointSetColor          Color;
  typedef FEVV::CGALPointSet::Property_map< Color > pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    return const_cast< FEVV::CGALPointSet & >(m)
        .add_property_map< Color >("vertex_color")
        .first;
  }
};


} // namespace FEVV
