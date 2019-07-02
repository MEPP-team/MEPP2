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


namespace FEVV {

// _PMap_traits specialization for FEVV::CGALPointSet
// use boost::vector_property_map


// define default property map type (associative property map,
// vector property map...) for each cell type (vertex, face...) ;

// default vertex property map is a vector property map
template< typename ValueT >
struct Vertex_pmap_traits< FEVV::CGALPointSet, ValueT >
{
  typedef typename boost::property_map< FEVV::CGALPointSet,
                                        boost::vertex_index_t >::const_type
      vertex_index_map_type;
  typedef typename boost::vector_property_map< ValueT, vertex_index_map_type >
      pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
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
  typedef typename FEVV::Geometry_traits< FEVV::CGALPointSet >::Vector
      value_type;
  typedef typename Vertex_pmap_traits< FEVV::CGALPointSet,
                                       value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for vertex-tangent
template<>
struct _PMap_traits< FEVV::CGALPointSet, FEVV::vertex_tangent_t >
{
  typedef typename FEVV::Geometry_traits< FEVV::CGALPointSet >::Vector
      value_type;
  typedef typename Vertex_pmap_traits< FEVV::CGALPointSet,
                                       value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for vertex-texcoord
template<>
struct _PMap_traits< FEVV::CGALPointSet, FEVV::vertex_texcoord_t >
{
  typedef typename FEVV::Geometry_traits< FEVV::CGALPointSet >::Vector
      value_type;
  typedef typename Vertex_pmap_traits< FEVV::CGALPointSet,
                                       value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the standard property map for vertex-color
template<>
struct _PMap_traits< FEVV::CGALPointSet, FEVV::vertex_color_t >
{
  typedef typename FEVV::Geometry_traits< FEVV::CGALPointSet >::Vector
      value_type;
  typedef typename Vertex_pmap_traits< FEVV::CGALPointSet,
                                       value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

#if 0 // useless for point cloud

// specialize the property maps traits for halfedge-normal
struct _PMap_traits< FEVV::CGALPointSet, FEVV::halfedge_normal_t >
{
  typedef typename FEVV::Geometry_traits< FEVV::CGALPointSet >::Vector
      value_type;
  typedef typename Halfedge_pmap_traits< FEVV::CGALPointSet,
                                         value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for halfedge-texcoord
struct _PMap_traits< FEVV::CGALPointSet, FEVV::halfedge_texcoord_t >
{
  typedef typename FEVV::Geometry_traits< FEVV::CGALPointSet >::Vector
      value_type;
  typedef typename Halfedge_pmap_traits< FEVV::CGALPointSet,
                                         value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for edge-color
struct _PMap_traits< FEVV::CGALPointSet, FEVV::edge_color_t >
{
  typedef typename FEVV::Geometry_traits< FEVV::CGALPointSet >::Vector
      value_type;
  typedef typename Edge_pmap_traits< FEVV::CGALPointSet,
                                     value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-normal
struct _PMap_traits< FEVV::CGALPointSet, FEVV::face_normal_t >
{
  typedef typename FEVV::Geometry_traits< FEVV::CGALPointSet >::Vector
      value_type;
  typedef typename Face_pmap_traits< FEVV::CGALPointSet,
                                     value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-color
struct _PMap_traits< FEVV::CGALPointSet, FEVV::face_color_t >
{
  typedef typename FEVV::Geometry_traits< FEVV::CGALPointSet >::Vector
      value_type;
  typedef typename Face_pmap_traits< FEVV::CGALPointSet,
                                     value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-material
struct _PMap_traits< FEVV::CGALPointSet, FEVV::face_material_t >
{
  typedef size_t value_type;
  typedef typename Face_pmap_traits< FEVV::CGALPointSet,
                                     value_type >::pmap_type pmap_type;

  static pmap_type create(const FEVV::CGALPointSet &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

#endif

} // namespace FEVV
