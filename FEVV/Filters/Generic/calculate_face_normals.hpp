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
#include "FEVV/Operators/Generic/calculate_face_normal.hpp"

namespace FEVV {
namespace Filters {
/**
 * \brief  Calculate "some" normal of all the faces of the considered
 *         mesh and populate the argument provided Face Normal Map
 *         with the resulting normals.
 *
 * \param[in]  g    The considered mesh
 * \param[in]  pm   The Point Map holding the vertices geometry
 * \param[out] fnm  The Property Map that will be populated with the resulting
 *                  normal vectors.
 * \param[in]  gt   The Geometry Traits used to perfom the geometrical
 *                  calculations
 *
 * \note: this function name is prefixed with 'calculate_' to mimic
 *        CGAL usage of the term 'compute'.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename FaceNormalMap,
          typename GeometryTraits >
void
calculate_face_normals(const HalfedgeGraph &g,
                       const PointMap &pm,
                       FaceNormalMap fnm,
                       const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;

  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::face_iterator face_iterator;

  face_iterator fb, fe;
  for(boost::tie(fb, fe) = faces(g); fb != fe; ++fb)
  {
    Vector vec = FEVV::Operators::
        calculate_face_normal< HalfedgeGraph, PointMap, GeometryTraits >(
            *fb, g, pm, gt);
    put(fnm, *fb, vec);
  }
}


/**
 * \brief  Compute all faces normals
 *
 * \param  g     mesh
 * \param  pm    point map
 * \param  fnm   face normal map
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename FaceNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
calculate_face_normals(const HalfedgeGraph &g,
                       const PointMap &pm,
                       FaceNormalMap fnm)
{
  GeometryTraits gt(g);
  return calculate_face_normals< HalfedgeGraph,
                                 PointMap,
                                 FaceNormalMap,
                                 GeometryTraits >(g, pm, fnm, gt);
}

} // namespace Filters
} // namespace FEVV
