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
#include <cassert>
#include <cmath>


namespace FEVV {

// Needs Graph_traits, Geometry_traits and Graph_properties
template<typename PointCloudT>
void
check_point_cloud_kNN_search_concept(void)
{
  // types
  typedef  boost::graph_traits< PointCloudT >         GraphTraits;
  typedef  typename GraphTraits::vertices_size_type   vertices_size_type;
  typedef  typename GraphTraits::vertex_descriptor    vertex_descriptor;
  typedef  typename GraphTraits::vertex_iterator      vertex_iterator;
  typedef  FEVV::Geometry_traits< PointCloudT >       Geometry;
  typedef  typename Geometry::Point                   Point;

  // create point cloud 
  PointCloudT cloud;

  // create points
  vertex_descriptor v0 = add_vertex(cloud); 
  vertex_descriptor v1 = add_vertex(cloud); 
  vertex_descriptor v2 = add_vertex(cloud); 
  vertex_descriptor v3 = add_vertex(cloud); 
  vertex_descriptor v4 = add_vertex(cloud); 
  vertex_descriptor v5 = add_vertex(cloud); 

  // set points geometry
  auto pm = get(boost::vertex_point, cloud);
  put(pm, v0, Point(0, 0, 0));
  put(pm, v1, Point(2, 2, 0));
  put(pm, v2, Point(1, 2, 0));
  put(pm, v3, Point(1, 1, 0));
  put(pm, v4, Point(0, 3, 0));
  put(pm, v5, Point(3, 0, 0));

  // do kNN-search
  unsigned int k = 3;
  auto kd_tree_ptr = create_kd_tree(cloud);
  auto result = kNN_search(*kd_tree_ptr, k, Point(0.6, 0.6, 0.), cloud);

  // check search result
  auto nn_ids = result.first;
  assert(nn_ids.size() == k);
  assert(nn_ids[0] == static_cast< vertex_descriptor >(3));
  assert(nn_ids[1] == static_cast< vertex_descriptor >(0));
  assert(nn_ids[2] == static_cast< vertex_descriptor >(2));
  auto nn_dists = result.second;
  assert(nn_dists.size() == k);
  assert(std::abs(nn_dists[0] - 0.565685) <= 0.0001);
  assert(std::abs(nn_dists[1] - 0.848528) <= 0.0001);
  assert(std::abs(nn_dists[2] - 1.45602) <= 0.0001);
}

} // namespace FEVV
