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

#include "FEVV/Wrappings/properties.h"
#include "FEVV/Wrappings/Wrappings_cgal_linear_cell_complex.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_linear_cell_complex.h"


#define CGAL_LCC_TEMPLATE_ARGS                                                 \
  unsigned int d_, unsigned int ambient_dim, class Traits_, class Items_,      \
      class Alloc_,                                                            \
      template< unsigned int, class, class, class, class >                     \
      class CMap, class Storage_

#define CGAL_LCC_TYPE                                                          \
  CGAL::Linear_cell_complex_for_combinatorial_map< d_,                         \
                                                   ambient_dim,                \
                                                   Traits_,                    \
                                                   Items_,                     \
                                                   Alloc_,                     \
                                                   CMap,                       \
                                                   Storage_ >


namespace FEVV {


// _PMap_traits specialization for CGAL::Linear_cell_complex
// use boost::vector_property_map


// define default property map type (associative property map,
// vector property map...) for each cell type (vertex, face...) ;

// default vertex property map
template< typename ValueT, CGAL_LCC_TEMPLATE_ARGS >
struct Vertex_pmap_traits< CGAL_LCC_TYPE, ValueT >
{
  typedef typename boost::property_map< CGAL_LCC_TYPE,
                                        boost::vertex_index_t >::const_type
      vertex_index_map_type;
  typedef typename boost::vector_property_map< ValueT, vertex_index_map_type >
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};


// default face property map
template< typename ValueT, CGAL_LCC_TEMPLATE_ARGS >
struct Face_pmap_traits< CGAL_LCC_TYPE, ValueT >
{
  typedef typename boost::property_map< CGAL_LCC_TYPE,
                                        boost::face_index_t >::const_type
      face_index_map_type;
  typedef typename boost::vector_property_map< ValueT, face_index_map_type >
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default edge property map
template< typename ValueT, CGAL_LCC_TEMPLATE_ARGS >
struct Edge_pmap_traits< CGAL_LCC_TYPE, ValueT >
{
  typedef typename boost::property_map< CGAL_LCC_TYPE,
                                        boost::edge_index_t >::const_type
      edge_index_map_type;
  typedef typename boost::vector_property_map< ValueT, edge_index_map_type >
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default halfedge property map
template< typename ValueT, CGAL_LCC_TEMPLATE_ARGS >
struct Halfedge_pmap_traits< CGAL_LCC_TYPE, ValueT >
{
  typedef typename boost::property_map< CGAL_LCC_TYPE,
                                        boost::halfedge_index_t >::const_type
      halfedge_index_map_type;
  typedef typename boost::vector_property_map< ValueT, halfedge_index_map_type >
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};


// define standard property map types for (property,cell) pair,
// for example vertex-normal

// specialize the property maps traits for vertex-normal
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::vertex_normal_t >
{
  typedef typename FEVV::Geometry_traits< CGAL_LCC_TYPE >::Vector value_type;
  typedef typename Vertex_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for vertex-tangent
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::vertex_tangent_t >
{
  typedef typename FEVV::Geometry_traits< CGAL_LCC_TYPE >::Vector value_type;
  typedef typename Vertex_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for vertex-texcoord
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::vertex_texcoord_t >
{
  typedef typename FEVV::Geometry_traits< CGAL_LCC_TYPE >::Vector value_type;
  typedef typename Vertex_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the standard property map for vertex-color
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::vertex_color_t >
{
  typedef typename FEVV::Geometry_traits< CGAL_LCC_TYPE >::Vector value_type;
  typedef typename Vertex_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for halfedge-normal
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::halfedge_normal_t >
{
  typedef typename FEVV::Geometry_traits< CGAL_LCC_TYPE >::Vector value_type;
  typedef typename Halfedge_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for halfedge-texcoord
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::halfedge_texcoord_t >
{
  typedef typename FEVV::Geometry_traits< CGAL_LCC_TYPE >::Vector value_type;
  typedef typename Halfedge_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for edge-color
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::edge_color_t >
{
  typedef typename FEVV::Geometry_traits< CGAL_LCC_TYPE >::Vector value_type;
  typedef typename Edge_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-normal
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::face_normal_t >
{
  typedef typename FEVV::Geometry_traits< CGAL_LCC_TYPE >::Vector value_type;
  typedef typename Face_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-color
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::face_color_t >
{
  typedef typename FEVV::Geometry_traits< CGAL_LCC_TYPE >::Vector value_type;
  typedef typename Face_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-material
template< CGAL_LCC_TEMPLATE_ARGS >
struct _PMap_traits< CGAL_LCC_TYPE, FEVV::face_material_t >
{
  typedef size_t value_type;
  typedef typename Face_pmap_traits< CGAL_LCC_TYPE, value_type >::pmap_type
      pmap_type;

  static pmap_type create(const CGAL_LCC_TYPE &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

} // namespace FEVV


#undef CGAL_LCC_TEMPLATE_ARGS
#undef CGAL_LCC_TYPE


