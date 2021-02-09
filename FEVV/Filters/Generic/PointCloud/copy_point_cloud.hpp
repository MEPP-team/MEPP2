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
#include "FEVV/Tools/Math/color_conversion.hpp"

#include <iostream>
#include <type_traits> // for std::remove_reference


namespace FEVV {
namespace Filters {


/**
 * Helper type for vertex to vertex map.
 */
template< typename PointCloudS,
          typename PointCloudT >
using PCVertexToVertexMap = std::map<
   typename boost::graph_traits< PointCloudS >::vertex_descriptor,
   typename boost::graph_traits< PointCloudT >::vertex_descriptor >;


/**
 * Implement the named parameters idiom for copy_point_cloud()
 * (for consistency with copy_graph())
 *
 * \param  v2v  vertex to vertex map to be populated
 */
template< typename PointCloudS,
          typename PointCloudT >
struct CopyPCParameters
{
  typedef  PCVertexToVertexMap< PointCloudS, PointCloudT >  V2VMap;

  CopyPCParameters& v2v(V2VMap *_v2v)
  {
    m_v2v = _v2v;
    return *this;
  }

  V2VMap *m_v2v = nullptr;
};


/**
 * \brief  Copy a source point cloud into a target point cloud.
 *         Copy standard properties too.
 *
 * \param  pc_s      the point cloud to copy from
 * \param  pmaps_s   the source property maps bag
 * \param  pc_t      the point cloud to copy to
 * \param  pmaps_t   the target property maps bag
 * \param  v2v       vertex to vertex map to populate (nullptr allowed)
 * \param  gt_s      the geometry traits of source point cloud
 * \param  gt_t      the geometry traits of target point cloud
 *
 * \sa     the simplified variant that use the default geometry traits.
 */
template< typename PointCloudS,
          typename PointCloudT,
          typename Parameters,
          typename GeometryTraitsS,
          typename GeometryTraitsT >
void
copy_point_cloud(const PointCloudS     &pc_s,
                 const PMapsContainer  &pmaps_s,
                 PointCloudT           &pc_t,
                 PMapsContainer        &pmaps_t,
                 const Parameters      &params,
                 const GeometryTraitsS &gt_s,
                 const GeometryTraitsT &gt_t)
{
  // target types
  typedef  typename GeometryTraitsT::Scalar  TScalar;
  typedef  typename GeometryTraitsT::Point   TPoint;
  typedef  typename GeometryTraitsT::Vector  TVector;
  typedef  decltype(make_property_map(FEVV::vertex_color, pc_t))  TColorMap;
  typedef  typename TColorMap::value_type    TColor;
  // note: we need the Color type because it is not a Vector with all
  //       datastructures ; the easiest way to define the Color type is
  //       "the type that is stored in the Color map", hence the declaration
  //       above

  // detect which standard properties are used by source
  bool use_vertex_normal = has_map(pmaps_s, vertex_normal);
  bool use_vertex_color  = has_map(pmaps_s, vertex_color);

  // retrieve point maps
  auto s_point_pm = get(boost::vertex_point, pc_s);
  auto t_point_pm = get(boost::vertex_point, pc_t);

  // create needed property maps in target
  if(use_vertex_normal)
  {
    auto vn_pm = make_property_map(FEVV::vertex_normal, pc_t);
    put_property_map(FEVV::vertex_normal, pc_t, pmaps_t, vn_pm);
  }
  if(use_vertex_color)
  {
    auto vc_pm = make_property_map(FEVV::vertex_color, pc_t);
    put_property_map(FEVV::vertex_color, pc_t, pmaps_t, vc_pm);
  }

  // reset v2v map
  if(params.m_v2v)
    params.m_v2v->clear();

  // copy geometry and properties
  auto iterator_pair = vertices(pc_s);
  auto s_vi = iterator_pair.first;
  auto s_vi_end = iterator_pair.second;
  for(; s_vi != s_vi_end; ++s_vi)
  {
    // create new vertex in target
    auto v = add_vertex(pc_t);

    // copy geometry
    auto s_point = get(s_point_pm, *s_vi);
    auto x = gt_s.get_x(s_point);
    auto y = gt_s.get_y(s_point);
    auto z = gt_s.get_z(s_point);
    put(t_point_pm, v, TPoint(static_cast< TScalar >(x),
                              static_cast< TScalar >(y),
                              static_cast< TScalar >(z)));

    // copy properties
    if(use_vertex_normal)
    {
      // vertex normal
      auto source_pm = get_property_map(FEVV::vertex_normal, pc_s, pmaps_s);
      auto target_pm = get_property_map(FEVV::vertex_normal, pc_t, pmaps_t);
      auto normal    = get(source_pm, *s_vi);
      put(target_pm, v, TVector(static_cast< TScalar >(normal[0]),
                                static_cast< TScalar >(normal[1]),
                                static_cast< TScalar >(normal[2])));
    }

    if(use_vertex_color)
    {
      // vertex color
      auto source_pm = get_property_map(FEVV::vertex_color, pc_s, pmaps_s);
      auto target_pm = get_property_map(FEVV::vertex_color, pc_t, pmaps_t);
      auto color     = get(source_pm, *s_vi);

      TColor t_color;
      typedef  typename std::remove_reference< decltype(t_color[0]) >::type
                                                              TColorValueType;

      put(target_pm,
          v,
          TColor(FEVV::Math::convert_color_value< TColorValueType >(color[0]),
                 FEVV::Math::convert_color_value< TColorValueType >(color[1]),
                 FEVV::Math::convert_color_value< TColorValueType >(color[2])));
    }

    // populate v2v map
    if(params.m_v2v)
      (*params.m_v2v)[*s_vi] = v;
  }

  // display some information

  std::cout << "copy_point_cloud(): input point_cloud has "
            << size_of_vertices(pc_s) << " vertices, "
            << std::endl;
  std::cout << "copy_point_cloud(): input point_cloud has property maps: [";
  for(auto &name: list_property_maps(pmaps_s))
    std::cout << " " << name;
  std::cout << " ]" << std::endl;

  std::cout << "copy_point_cloud(): output point_cloud has "
            << size_of_vertices(pc_t) << " vertices, "
            << std::endl;
  std::cout << "copy_point_cloud(): output point_cloud has property maps: [";
  for(auto &name: list_property_maps(pmaps_t))
    std::cout << " " << name;
  std::cout << " ]" << std::endl;
}


/**
 * \brief  Copy a source point cloud into a target point cloud.
 *         Copy standard properties too.
 *
 * \param  pc_s      the point cloud to copy from
 * \param  pmaps_s   the source property maps bag
 * \param  pc_t      the point cloud to copy to
 * \param  pmaps_t   the target property maps bag
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename PointCloudS,
          typename PointCloudT,
          typename Parameters,
          typename GeometryTraitsS = FEVV::Geometry_traits< PointCloudS >,
          typename GeometryTraitsT = FEVV::Geometry_traits< PointCloudT > >
void
copy_point_cloud(const PointCloudS     &pc_s,
                 const PMapsContainer  &pmap_s,
                 PointCloudT           &pc_t,
                 PMapsContainer        &pmap_t,
                 const Parameters      &params =
                    CopyPCParameters< PointCloudS, PointCloudT >())
{
  GeometryTraitsS gt_s(pc_s);
  GeometryTraitsT gt_t(pc_t);
  copy_point_cloud(pc_s, pmap_s, pc_t, pmap_t, params, gt_s, gt_t);
}


} // namespace Filters
} // namespace FEVV
