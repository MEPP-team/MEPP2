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
#include "FEVV/Wrappings/Wrappings_pcl_point_cloud.h"
#include "FEVV/Wrappings/Geometry_traits_pcl_point_cloud.h"

#include <boost/property_map/property_map.hpp> // for boost property maps


namespace FEVV {

// _PMap_traits specialization for FEVV::PCLPointCloud
// use boost::vector_property_map except for NormalMap
// and ColorMap ; NormalMap and ColorMap are mapping
// over PCL datastructure


// define default property map type (associative property map,
// vector property map...) for each cell type (vertex, face...) ;

// vertex desciptor already is an index, so we can use it
// as boost::vector_property_map key
inline
boost::identity_property_map
get(const boost::vertex_index_t&, const FEVV::PCLPointCloud&)
{
  return boost::identity_property_map();
}

// default vertex property map
template< typename ValueT >
struct Vertex_pmap_traits< FEVV::PCLPointCloud, ValueT >
{
  typedef typename boost::property_map< FEVV::PCLPointCloud,
                                        boost::vertex_index_t >::const_type
      vertex_index_map_type;
  typedef typename boost::vector_property_map< ValueT, vertex_index_map_type >
      pmap_type;

  static pmap_type create(const FEVV::PCLPointCloud &pc)
  {
    auto index_map = get(boost::vertex_index, pc);
    pmap_type pmap(index_map);
    return pmap;
  }
};


// define standard property map types for (property,cell) pair,
// for example vertex-normal

// specialize the property maps traits for vertex-normal
template<>
struct _PMap_traits< FEVV::PCLPointCloud, FEVV::vertex_normal_t >
{
  typedef FEVV::PCLPointCloudNormalMap   pmap_type;

  static pmap_type create(const FEVV::PCLPointCloud &pc)
  {
    return pmap_type(const_cast< FEVV::PCLPointCloud& >(pc));
  }
};

#if 0 //TODO-elo-rm
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
#endif

// specialize the standard property map for vertex-color
template<>
struct _PMap_traits< FEVV::PCLPointCloud, FEVV::vertex_color_t >
{
  typedef FEVV::PCLPointCloudColorMap   pmap_type;

  static pmap_type create(const FEVV::PCLPointCloud &pc)
  {
    return pmap_type(const_cast< FEVV::PCLPointCloud& >(pc));
  }
};

} // namespace FEVV
