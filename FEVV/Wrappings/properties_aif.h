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
#include "FEVV/Wrappings/Wrappings_aif.h"
#include "FEVV/Wrappings/Geometry_traits_aif.h"

//#include "FEVV/DataStructures/DataStructures_aif.h"

#define USE_NATIVE_PM
#ifdef USE_NATIVE_PM
#include "FEVV/DataStructures/AIF/AIFProperties.h"
#endif

namespace FEVV {

// Here is the place for some specialization for AIF

// default vertex property map
template< typename ValueT >
struct Vertex_pmap_traits< FEVV::DataStructures::AIF::AIFMesh, ValueT >
{
  typedef typename boost::property_map< FEVV::DataStructures::AIF::AIFMesh,
                                        boost::vertex_index_t >::const_type
      vertex_index_map_type;
  typedef typename boost::vector_property_map< ValueT, vertex_index_map_type >
      pmap_type;

  static pmap_type create(const FEVV::DataStructures::AIF::AIFMesh &m)
  {
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

// default face property map
template< typename ValueT >
struct Face_pmap_traits< AIFMesh, ValueT >
{
  typedef
      typename boost::property_map< AIFMesh, boost::face_index_t >::const_type
          face_index_map_type;
  typedef typename boost::vector_property_map< ValueT, face_index_map_type >
      pmap_type;

  static pmap_type create(const AIFMesh &m)
  {
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};


// default edge property map
template< typename ValueT >
struct Edge_pmap_traits< AIFMesh, ValueT >
{
  typedef
      typename boost::property_map< AIFMesh, boost::edge_index_t >::const_type
          edge_index_map_type;
  typedef typename boost::vector_property_map< ValueT, edge_index_map_type >
      pmap_type;

  static pmap_type create(const AIFMesh &m)
  {
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
  }
};

///////////////////////////////////////////////////////////////////////////////
// define standard property map types for (property,cell) pair,
// for example vertex-normal

// specialize the property maps traits for vertex-normal
template<>
struct _PMap_traits< AIFMesh, FEVV::vertex_normal_t >
{
  typedef typename FEVV::Geometry_traits< AIFMesh >::Vector value_type;
#ifdef USE_NATIVE_PM
  typedef typename FEVV::DataStructures::AIF::PropertyMap< value_type >
      pmap_type;
#else
  typedef
      typename Vertex_pmap_traits< AIFMesh, value_type >::pmap_type pmap_type;
#endif
  static pmap_type create(const AIFMesh &m)
  {
#ifdef USE_NATIVE_PM
    pmap_type *pmap =
        (const_cast< AIFMesh * >(&m))
            ->AddPropertyMap< AIFMesh::vertex_type::ptr, value_type >(
                "v:normal");
    return *pmap;
#else
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
#endif
  }
};

// specialize the property maps traits for vertex-tangent
template<>
struct _PMap_traits< AIFMesh, FEVV::vertex_tangent_t >
{
  typedef typename FEVV::Geometry_traits< AIFMesh >::Vector value_type;
#ifdef USE_NATIVE_PM
  typedef typename FEVV::DataStructures::AIF::PropertyMap< value_type >
      pmap_type;
#else
  typedef
      typename Vertex_pmap_traits< AIFMesh, value_type >::pmap_type pmap_type;
#endif
  static pmap_type create(const AIFMesh &m)
  {
#ifdef USE_NATIVE_PM
    pmap_type *pmap =
        (const_cast< AIFMesh * >(&m))
            ->AddPropertyMap< AIFMesh::vertex_type::ptr, value_type >(
                "v:tangent");
    return *pmap;
#else
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
#endif
  }
};

// specialize the property maps traits for vertex-texcoord
template<>
struct _PMap_traits< AIFMesh, FEVV::vertex_texcoord_t >
{
  typedef typename FEVV::Geometry_traits< AIFMesh >::Vector value_type;
#ifdef USE_NATIVE_PM
  typedef typename FEVV::DataStructures::AIF::PropertyMap< value_type >
      pmap_type;
#else
  typedef
      typename Vertex_pmap_traits< AIFMesh, value_type >::pmap_type pmap_type;
#endif
  static pmap_type create(const AIFMesh &m)
  {
#ifdef USE_NATIVE_PM
    pmap_type *pmap =
        (const_cast< AIFMesh * >(&m))
            ->AddPropertyMap< AIFMesh::vertex_type::ptr, value_type >(
                "v:texcoord");
    return *pmap;
#else
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
#endif
  }
};

// specialize the standard property map for vertex-color
template<>
struct _PMap_traits< AIFMesh, FEVV::vertex_color_t >
{
  typedef typename FEVV::Geometry_traits< AIFMesh >::Vector value_type;
#ifdef USE_NATIVE_PM
  typedef typename FEVV::DataStructures::AIF::PropertyMap< value_type >
      pmap_type;
#else
  typedef
      typename Vertex_pmap_traits< AIFMesh, value_type >::pmap_type pmap_type;
#endif
  static pmap_type create(const AIFMesh &m)
  {
#ifdef USE_NATIVE_PM
    pmap_type *pmap =
        (const_cast< AIFMesh * >(&m))
            ->AddPropertyMap< AIFMesh::vertex_type::ptr, value_type >(
                "v:color");
    return *pmap;
#else
    auto index_map = get(boost::vertex_index, m);
    pmap_type pmap(index_map);
    return pmap;
#endif
  }
};

// specialize the property maps traits for edge-color
template<>
struct _PMap_traits< AIFMesh, FEVV::edge_color_t >
{
  typedef typename FEVV::Geometry_traits< AIFMesh >::Vector value_type;
#ifdef USE_NATIVE_PM
  typedef typename FEVV::DataStructures::AIF::PropertyMap< value_type >
      pmap_type;
#else
  typedef typename Edge_pmap_traits< AIFMesh, value_type >::pmap_type pmap_type;
#endif
  static pmap_type create(const AIFMesh &m)
  {
#ifdef USE_NATIVE_PM
    pmap_type *pmap =
        (const_cast< AIFMesh * >(&m))
            ->AddPropertyMap< AIFMesh::edge_type::ptr, value_type >("e:color");
    return *pmap;
#else
    auto index_map = get(boost::edge_index, m);
    pmap_type pmap(index_map);
    return pmap;
#endif
  }
};


// specialize the property maps traits for face-normal
template<>
struct _PMap_traits< AIFMesh, FEVV::face_normal_t >
{
  typedef typename FEVV::Geometry_traits< AIFMesh >::Vector value_type;
#ifdef USE_NATIVE_PM
  typedef typename FEVV::DataStructures::AIF::PropertyMap< value_type >
      pmap_type;
#else
  typedef typename Face_pmap_traits< AIFMesh, value_type >::pmap_type pmap_type;
#endif
  static pmap_type create(const AIFMesh &m)
  {
#ifdef USE_NATIVE_PM
    pmap_type *pmap =
        (const_cast< AIFMesh * >(&m))
            ->AddPropertyMap< AIFMesh::face_type::ptr, value_type >("f:normal");
    return *pmap;
#else
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
#endif
  }
};

// specialize the property maps traits for face-color
template<>
struct _PMap_traits< AIFMesh, FEVV::face_color_t >
{
  typedef typename FEVV::Geometry_traits< AIFMesh >::Vector value_type;
#ifdef USE_NATIVE_PM
  typedef typename FEVV::DataStructures::AIF::PropertyMap< value_type >
      pmap_type;
#else
  typedef typename Face_pmap_traits< AIFMesh, value_type >::pmap_type pmap_type;
#endif
  static pmap_type create(const AIFMesh &m)
  {
#ifdef USE_NATIVE_PM
    pmap_type *pmap =
        (const_cast< AIFMesh * >(&m))
            ->AddPropertyMap< AIFMesh::face_type::ptr, value_type >("f:color");
    return *pmap;
#else
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
#endif
  }
};

// specialize the property maps traits for face-material
template<>
struct _PMap_traits< AIFMesh, FEVV::face_material_t >
{
  typedef size_t value_type;
#ifdef USE_NATIVE_PM
  typedef typename FEVV::DataStructures::AIF::PropertyMap< value_type >
      pmap_type;
#else
  typedef typename Face_pmap_traits< AIFMesh, value_type >::pmap_type pmap_type;
#endif
  static pmap_type create(const AIFMesh &m)
  {
#ifdef USE_NATIVE_PM
    pmap_type *pmap =
        (const_cast< AIFMesh * >(&m))
            ->AddPropertyMap< AIFMesh::face_type::ptr, value_type >(
                "f:material");
    return *pmap;
#else
    auto index_map = get(boost::face_index, m);
    pmap_type pmap(index_map);
    return pmap;
#endif
  }
};

} // namespace FEVV

