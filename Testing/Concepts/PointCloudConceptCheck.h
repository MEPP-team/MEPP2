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


namespace FEVV {


template<typename PointCloudT>
void
check_point_cloud_concept(PointCloudT &pc)
{
  // types
  typedef  boost::graph_traits< PointCloudT >         GraphTraits;
  typedef  typename GraphTraits::vertices_size_type   vertices_size_type;
  typedef  typename GraphTraits::vertex_descriptor    vertex_descriptor;
  typedef  typename GraphTraits::vertex_iterator      vertex_iterator;

  // null_vertex()
  vertex_descriptor vd_null = GraphTraits::null_vertex();

  // free functions

  // vertices(pc)
  std::pair< vertex_iterator, vertex_iterator > iterator_pair;
  iterator_pair = vertices(pc);

  // num_vertices(pc)
  vertices_size_type n = num_vertices(pc);

  // add_vertex(pc)
  vertex_descriptor vd = add_vertex(pc);

  // remove_vertex(vd, pc)
  remove_vertex(vd, pc);
}


} // namespace FEVV
