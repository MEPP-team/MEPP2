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

#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/Wrappings/Graph_traits_aif.h"

#include <CGAL/boost/graph/properties.h> // for boost::vertex_point_t


typedef FEVV::DataStructures::AIF::AIFMesh AIFMesh;
typedef FEVV::DataStructures::AIF::AIFVertex AIFVertex;
typedef AIFMesh::Point Point;
typedef AIFMesh::Vector Vector;


/*
 * A class that associate a vertex/edge/face descriptor to an index,
 * that will be used as a property map key.
 */
template< typename VEF >
class AIF_index_pmap
    : public boost::put_get_helper< unsigned int, AIF_index_pmap< VEF > >
{
public:
  typedef boost::readable_property_map_tag category;
  // typedef boost::read_write_property_map_tag  category;
  // typedef boost::lvalue_property_map_tag    category;
  typedef unsigned int value_type;
  typedef unsigned int reference;
  typedef VEF key_type;

  value_type operator[](const key_type &vd) const
  {
    return static_cast< value_type >(vd->GetIndex());
  }
};


namespace boost {

// FACE
template<>
struct property_map< FEVV::DataStructures::AIF::AIFMesh, boost::face_index_t >
{
  typedef AIF_index_pmap< boost::graph_traits< AIFMesh >::face_descriptor >
      type;
  typedef AIF_index_pmap< boost::graph_traits< AIFMesh >::face_descriptor >
      const_type;
};
// FACE

// VERTEX
template<>
struct property_map< FEVV::DataStructures::AIF::AIFMesh, boost::vertex_index_t >
{
  typedef AIF_index_pmap< boost::graph_traits< AIFMesh >::vertex_descriptor >
      type;
  typedef AIF_index_pmap< boost::graph_traits< AIFMesh >::vertex_descriptor >
      const_type;
};
template<>
struct property_map< FEVV::DataStructures::AIF::AIFMesh, boost::vertex_point_t >
{
  typedef FEVV::DataStructures::AIF::PropertyMap< Point > type;
  typedef FEVV::DataStructures::AIF::PropertyMap< Point > const_type;
};
// VERTEX

// HALFEDGE
template<>
struct property_map< FEVV::DataStructures::AIF::AIFMesh,
                     boost::halfedge_index_t >
{
  typedef AIF_index_pmap< boost::graph_traits< AIFMesh >::halfedge_descriptor >
      type;
  typedef AIF_index_pmap< boost::graph_traits< AIFMesh >::halfedge_descriptor >
      const_type;
};
// HALFEDGE

// EDGE
template<>
struct property_map< FEVV::DataStructures::AIF::AIFMesh, boost::edge_index_t >
{
  typedef AIF_index_pmap< boost::graph_traits< AIFMesh >::edge_descriptor >
      type;
  typedef AIF_index_pmap< boost::graph_traits< AIFMesh >::edge_descriptor >
      const_type;
};
// EDGE

/*****************************************************************************/

/**
 * \brief  Returns the vertex index property map of the mesh.
 */
inline typename boost::property_map< FEVV::DataStructures::AIF::AIFMesh,
                                     boost::vertex_index_t >::const_type
get(const boost::vertex_index_t &, const FEVV::DataStructures::AIF::AIFMesh &)
{
  return AIF_index_pmap< typename boost::graph_traits<
      const FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor >();
}

/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 */
inline typename FEVV::DataStructures::AIF::PropertyMap< Point > &
get(const boost::vertex_point_t, FEVV::DataStructures::AIF::AIFMesh &m)
{
  FEVV::DataStructures::AIF::PropertyMap< Point > *pm =
      m.GetPropertyMap< AIFVertex::ptr, Point >("v:point");
  return *pm;
}

/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 *         Const version.
 */
inline const typename FEVV::DataStructures::AIF::PropertyMap< Point > &
get(const boost::vertex_point_t, const FEVV::DataStructures::AIF::AIFMesh &m)
{
  FEVV::DataStructures::AIF::PropertyMap< Point > *pm =
      const_cast< FEVV::DataStructures::AIF::AIFMesh & >(m)
          .GetPropertyMap< AIFVertex::ptr, Point >("v:point");
  return *pm;
}

/*****************************************************************************/

/**
 * \brief  Returns the face index property map of the mesh.
 */
inline typename boost::property_map< FEVV::DataStructures::AIF::AIFMesh,
                                     boost::face_index_t >::const_type
get(const boost::face_index_t &, const FEVV::DataStructures::AIF::AIFMesh &)
{
  return AIF_index_pmap< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::face_descriptor >();
}

/*****************************************************************************/

/**
 * \brief  Returns the edge index property map of the mesh.
 */
inline typename boost::property_map< FEVV::DataStructures::AIF::AIFMesh,
                                     boost::edge_index_t >::const_type
get(const boost::edge_index_t &, const FEVV::DataStructures::AIF::AIFMesh &)
{
  return AIF_index_pmap< typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::edge_descriptor >();
}

/*****************************************************************************/
} // namespace boost


namespace FEVV {
namespace DataStructures {
namespace AIF {

//! Specialization of get(pmap, key) for AIF
template< typename ValueT, typename KeyT >
inline ValueT &
get(const typename FEVV::DataStructures::AIF::PropertyMap< ValueT > &pm, KeyT k)
{
  return pm[k];
}

/*****************************************************************************/

//! \brief  Specialization of get(pmap, mesh, key) for AIF.
inline typename boost::property_traits< typename boost::property_map<
    FEVV::DataStructures::AIF::AIFMesh,
    boost::vertex_index_t >::const_type >::value_type
get(const boost::vertex_index_t &prop,
    const FEVV::DataStructures::AIF::AIFMesh &sm,
    const typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor k)
{
  auto pm_index = boost::get(prop, sm);
  return get(pm_index, k);
}

/*****************************************************************************/

//! \brief  Specialization of get(pmap, mesh, key) for AIF.
inline Point
get(const boost::vertex_point_t &prop,
    const FEVV::DataStructures::AIF::AIFMesh &sm,
    typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor k)
{
  typename FEVV::DataStructures::AIF::PropertyMap< Point > pm_pos =
      boost::get(prop, sm);
  return get(pm_pos, k);
}

/*****************************************************************************/

//! \brief  Specialization of put(pmap, key, value) for AIF
template< typename ValueT, typename KeyT >
inline void
put(const typename FEVV::DataStructures::AIF::PropertyMap< ValueT > &pm,
    KeyT k,
    const ValueT &v)
{
  pm[k] = v;
}

/*****************************************************************************/

//! \brief  Specialization of put(pmap, mesh, key, val) for AIF
template< typename ValueT >
inline void
put(const boost::vertex_index_t &prop,
    const FEVV::DataStructures::AIF::AIFMesh &sm,
    const typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor k,
    // const typename boost::property_map<typename
    // FEVV::DataStructures::AIF::AIFMesh,
    // boost::vertex_index_t>::const_type>::value_type& v
    const ValueT &v)
{
  typedef typename boost::property_map< FEVV::DataStructures::AIF::AIFMesh,
                                        boost::vertex_index_t >::type Map;
  Map pm = get(prop, sm);
  // auto pm = boost::get(prop, sm);
  put(pm, k, v);
}

/*****************************************************************************/

//! \brief  Specialization of put(pmap, mesh, key, val) for AIF.
inline void
put(const boost::vertex_point_t &prop,
    const FEVV::DataStructures::AIF::AIFMesh &sm,
    const typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor k,
    const typename boost::graph_traits<
        FEVV::DataStructures::AIF::AIFMesh >::vertex_property_type &vp)
{
  assert(false);
  // auto pos_pm = boost::get(boost::vertex_point, sm);
  // boost::put(pos_pm, k, vp);
}

/*****************************************************************************/

// boost MutablePropertyGraph Concept
inline typename boost::graph_traits<
    FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor
add_vertex(typename boost::graph_traits<
               FEVV::DataStructures::AIF::AIFMesh >::vertex_property_type vp,
           FEVV::DataStructures::AIF::AIFMesh &sm)
{
  // typedef FEVV::DataStructures::AIF::AIFTopologyHelpers AIFHelpers; //
  // AIFHelpers type is defined in Graph_traits_aif.h
  typename boost::graph_traits<
      FEVV::DataStructures::AIF::AIFMesh >::vertex_descriptor v =
      AIFHelpers::add_vertex(sm);
  auto pm = boost::get(boost::vertex_point, sm);
  put(pm, v, vp);
  return v;
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV

