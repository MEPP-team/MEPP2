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
#include "FEVV/DataStructures/DataStructures_openmesh.h"
#include "FEVV/Wrappings/Geometry_traits_openmesh.h"

namespace FEVV {

// _PMap_traits specialization for Openmesh
// use boost::vector_property_map

// template<typename OMTraits>
// struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT<OMTraits> >

// define default property map type (associative property map,
// vector property map...) for each cell type (vertex, face...) ;

// default vertex property map
template< typename ValueT, typename OMTraits >
struct Vertex_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >, ValueT >
{
  typedef
      typename boost::property_map< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                    boost::vertex_index_t >::const_type
          vertex_index_map_type;
  typedef typename boost::vector_property_map< ValueT, vertex_index_map_type >
      pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default face property map
template< typename ValueT, typename OMTraits >
struct Face_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >, ValueT >
{
  typedef
      typename boost::property_map< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                    boost::face_index_t >::const_type
          face_index_map_type;
  typedef typename boost::vector_property_map< ValueT, face_index_map_type >
      pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default edge property map
template< typename ValueT, typename OMTraits >
struct Edge_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >, ValueT >
{
  typedef
      typename boost::property_map< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                    boost::edge_index_t >::const_type
          edge_index_map_type;
  typedef typename boost::vector_property_map< ValueT, edge_index_map_type >
      pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default halfedge property map
template< typename ValueT, typename OMTraits >
struct Halfedge_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                             ValueT >
{
  typedef
      typename boost::property_map< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                    boost::halfedge_index_t >::const_type
          halfedge_index_map_type;
  typedef typename boost::vector_property_map< ValueT, halfedge_index_map_type >
      pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};


// define standard property map types for (property,cell) pair,
// for example vertex-normal

// specialize the property maps traits for vertex-normal

template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::vertex_normal_t >
{
  typedef typename FEVV::Geometry_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits > >::Vector value_type;
  typedef
      typename Vertex_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                   value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for vertex-tangent
template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::vertex_tangent_t >
{
  typedef typename FEVV::Geometry_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits > >::Vector value_type;
  typedef
      typename Vertex_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                   value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for vertex-texcoord
template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::vertex_texcoord_t >
{
  typedef typename FEVV::Geometry_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits > >::Vector value_type;
  typedef
      typename Vertex_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                   value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the standard property map for vertex-color
template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::vertex_color_t >
{
  typedef typename FEVV::Geometry_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits > >::Vector value_type;
  typedef
      typename Vertex_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                   value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for halfedge-normal
template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::halfedge_normal_t >
{
  typedef typename FEVV::Geometry_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits > >::Vector value_type;
  typedef typename Halfedge_pmap_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
      value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for halfedge-texcoord
template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::halfedge_texcoord_t >
{
  typedef typename FEVV::Geometry_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits > >::Vector value_type;
  typedef typename Halfedge_pmap_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
      value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::halfedge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for edge-color
template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::edge_color_t >
{
  typedef typename FEVV::Geometry_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits > >::Vector value_type;
  typedef
      typename Edge_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                 value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-normal
template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::face_normal_t >
{
  typedef typename FEVV::Geometry_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits > >::Vector value_type;
  typedef
      typename Face_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                 value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-color
template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::face_color_t >
{
  typedef typename FEVV::Geometry_traits<
      OpenMesh::PolyMesh_ArrayKernelT< OMTraits > >::Vector value_type;
  typedef
      typename Face_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                 value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// specialize the property maps traits for face-material
template< typename OMTraits >
struct _PMap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                     FEVV::face_material_t >
{
  typedef size_t value_type;
  typedef
      typename Face_pmap_traits< OpenMesh::PolyMesh_ArrayKernelT< OMTraits >,
                                 value_type >::pmap_type pmap_type;

  static pmap_type create(const OpenMesh::PolyMesh_ArrayKernelT< OMTraits > &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};


} // namespace FEVV

