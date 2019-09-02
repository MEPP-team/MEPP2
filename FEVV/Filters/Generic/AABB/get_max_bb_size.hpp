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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"

#include "FEVV/Filters/Generic/AABB/compute_mesh_bounding_box.hpp"


namespace FEVV {
namespace Filters {

/**
 * \brief   Gets the maximum size of the object bounding box.
 *
 * \param   g    input mesh
 * \param   pm   point map
 * \param   gt   geometry traits
 *
 * \return  the maximum size of the bounding box.
*/
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
double
get_max_bb_size(const HalfedgeGraph &g,
             PointMap &pm,
             GeometryTraits &gt)
{
  typename boost::property_traits< PointMap >::value_type  min_aabb,
                                                           max_aabb;
  Filters::compute_mesh_bounding_box< HalfedgeGraph,
                                    PointMap,
                                    GeometryTraits >(g,
                                                     pm,
                                                     min_aabb,
                                                     max_aabb,
                                                     gt);

  double xmax = gt.get_x(max_aabb) - gt.get_x(min_aabb);
  double ymax = gt.get_y(max_aabb) - gt.get_y(min_aabb);
  double zmax = gt.get_z(max_aabb) - gt.get_z(min_aabb);
  double max = xmax > ymax ? xmax : ymax;

  return max > zmax ? max : zmax;
}

/**
 * \brief   Gets the maximum size of the object bounding box.
 *          Use the default geometry traits of the mesh.
 *
 * \param   g    input mesh
 * \param   pm   point map
 *
 * \return  the maximum size of the bounding box.
*/
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
double
get_max_bb_size(const HalfedgeGraph &g,
             PointMap &pm)
{
  GeometryTraits gt(g);
  return get_max_bb_size< HalfedgeGraph,
                       PointMap,
                       GeometryTraits >(g, pm, gt);
}

} //namespace Filters
} //namespace FEVV
