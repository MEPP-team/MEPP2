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
#include "FEVV/Wrappings/Wrappings_cgal_surface_mesh.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"

namespace FEVV {

// _PMap_traits specialization for CGAL::Surface_mesh
// use boost::vector_property_map


// define default property map type (associative property map,
// vector property map...) for each cell type (vertex, face...) ;

// default vertex property map
template< typename ValueT, typename PointT >
struct Vertex_pmap_traits< CGAL::Surface_mesh< PointT >, ValueT >
{
  typedef typename boost::property_map< CGAL::Surface_mesh< PointT >,
                                        boost::vertex_index_t >::const_type
      vertex_index_map_type;
  typedef typename boost::vector_property_map< ValueT, vertex_index_map_type >
      pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};


// default face property map
template< typename ValueT, typename PointT >
struct Face_pmap_traits< CGAL::Surface_mesh< PointT >, ValueT >
{
  typedef typename boost::property_map< CGAL::Surface_mesh< PointT >,
                                        boost::face_index_t >::const_type
      face_index_map_type;
  typedef typename boost::vector_property_map< ValueT, face_index_map_type >
      pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default edge property map
template< typename ValueT, typename PointT >
struct Edge_pmap_traits< CGAL::Surface_mesh< PointT >, ValueT >
{
  typedef typename boost::property_map< CGAL::Surface_mesh< PointT >,
                                        boost::edge_index_t >::const_type
      edge_index_map_type;
  typedef typename boost::vector_property_map< ValueT, edge_index_map_type >
      pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default halfedge property map
template< typename ValueT, typename PointT >
struct Halfedge_pmap_traits< CGAL::Surface_mesh< PointT >, ValueT >
{
  typedef typename boost::property_map< CGAL::Surface_mesh< PointT >,
                                        boost::halfedge_index_t >::const_type
      halfedge_index_map_type;
  typedef typename boost::vector_property_map< ValueT, halfedge_index_map_type >
      pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};


// define standard property map types for (property,cell) pair,
// for example vertex-normal

// specialize the property maps traits for vertex-normal
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::vertex_normal_t >
{
  typedef typename FEVV::Geometry_traits< CGAL::Surface_mesh< PointT > >::Vector
      value_type;
  typedef typename Vertex_pmap_traits< CGAL::Surface_mesh< PointT >,
                                       value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for vertex-tangent
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::vertex_tangent_t >
{
  typedef typename FEVV::Geometry_traits< CGAL::Surface_mesh< PointT > >::Vector
      value_type;
  typedef typename Vertex_pmap_traits< CGAL::Surface_mesh< PointT >,
                                       value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for vertex-texcoord
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::vertex_texcoord_t >
{
  typedef typename FEVV::Geometry_traits< CGAL::Surface_mesh< PointT > >::Vector
      value_type;
  typedef typename Vertex_pmap_traits< CGAL::Surface_mesh< PointT >,
                                       value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the standard property map for vertex-color
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::vertex_color_t >
{
  typedef typename FEVV::Geometry_traits< CGAL::Surface_mesh< PointT > >::Vector
      value_type;
  typedef typename Vertex_pmap_traits< CGAL::Surface_mesh< PointT >,
                                       value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for halfedge-normal
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::halfedge_normal_t >
{
  typedef typename FEVV::Geometry_traits< CGAL::Surface_mesh< PointT > >::Vector
      value_type;
  typedef typename Halfedge_pmap_traits< CGAL::Surface_mesh< PointT >,
                                         value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for halfedge-texcoord
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::halfedge_texcoord_t >
{
  typedef typename FEVV::Geometry_traits< CGAL::Surface_mesh< PointT > >::Vector
      value_type;
  typedef typename Halfedge_pmap_traits< CGAL::Surface_mesh< PointT >,
                                         value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for edge-color
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::edge_color_t >
{
  typedef typename FEVV::Geometry_traits< CGAL::Surface_mesh< PointT > >::Vector
      value_type;
  typedef typename Edge_pmap_traits< CGAL::Surface_mesh< PointT >,
                                     value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-normal
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::face_normal_t >
{
  typedef typename FEVV::Geometry_traits< CGAL::Surface_mesh< PointT > >::Vector
      value_type;
  typedef typename Face_pmap_traits< CGAL::Surface_mesh< PointT >,
                                     value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-color
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::face_color_t >
{
  typedef typename FEVV::Geometry_traits< CGAL::Surface_mesh< PointT > >::Vector
      value_type;
  typedef typename Face_pmap_traits< CGAL::Surface_mesh< PointT >,
                                     value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-material
template< typename PointT >
struct _PMap_traits< CGAL::Surface_mesh< PointT >, FEVV::face_material_t >
{
  typedef size_t value_type;
  typedef typename Face_pmap_traits< CGAL::Surface_mesh< PointT >,
                                     value_type >::pmap_type pmap_type;

  static pmap_type create(const CGAL::Surface_mesh< PointT > &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};


} // namespace FEVV

