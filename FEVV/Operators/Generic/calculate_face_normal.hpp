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

namespace FEVV {
namespace Operators {

/**
 * \brief  Calculate "some" normal to the considered face.
 *
 * \param[in] f     The considered face
 * \param[in] g     The mesh to which the face belongs
 * \param[in] pm    The Point Map holding the vertices geometry
 * \param[in] gt    The Geometry Traits associated to the mesh
 *
 * \note: this function name is prefixed with 'calculate_' to mimic
 *        CGAL usage of the term 'compute'.
 */
template< typename HalfedgeGraph, typename PointMap, typename GeometryTraits >
typename GeometryTraits::Vector
calculate_face_normal(
    typename boost::graph_traits< HalfedgeGraph >::face_descriptor f,
    const HalfedgeGraph &g,
    const PointMap &pm,
    const GeometryTraits &gt)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GeometryTraits::Point Point;
  typedef typename GeometryTraits::Vector Vector;

  halfedge_descriptor edg = halfedge(f, g);
  halfedge_descriptor edgb = edg;

  Point p0 = get(pm, target(edg, g));
  edg = next(edg, g);
  Point p1 = get(pm, target(edg, g));
  edg = next(edg, g);
  Point p2 = get(pm, target(edg, g));
  edg = next(edg, g);

  if(edg == edgb)
  {
    // triangle
    try
    {
      Vector normal = gt.unit_normal(p1, p2, p0);
      return normal;
    }
    catch(const std::exception &e)
    {
      //  When 2 edges of the face are collinear an exception is triggered, in
      //  debug mode only. Then the line:
      //     std::cout << "e.what() = '" << e.what() << "'" << std::endl;
      //  produces the following output:
      //     e.what() = 'CGAL ERROR: precondition violation!
      //     Expr: ! K().collinear_3_object()(p,q,r)
      //     File: /home/eric/CGAL-4.11/include/CGAL/Kernel/function_objects.h
      //     Line: 2490'

      if(std::string(e.what()).find("collinear_3_object") == std::string::npos)
      {
        // This is NOT the true 'collinear' exception. Re-throw the exception
        // unchanged.
        throw;
      }

      // This is the true 'collinear' exception. We set normal vector to be
      // (NaN, NaN, NaN) and hope the caller will find out that the normal
      // is snafu.
      Vector normal(std::numeric_limits< double >::quiet_NaN(),
                    std::numeric_limits< double >::quiet_NaN(),
                    std::numeric_limits< double >::quiet_NaN());

      return normal;
    }
  }
  else
  { // not a triangle but a quad, a polygon

    // The correct algorithme would be to compute the optimal plane fitting
    // all the face vertices. However this would be too computationally
    // expensive for such a common operation and hence we fold back to
    // computing the mean of the normals of the triangles build out of
    // triplets of vertices "skifully" taken among the vertices of the face.
    Vector n(gt.normal(p1, p2, p0));

    p0 = p2; // do not move P0 in the following loop
    do
    {
      edg = next(edg, g);
      p1 = get(pm, source(edg, g));
      p2 = get(pm, target(edg, g));

      n = gt.add(n, gt.normal(p1, p2, p0));
    } while(edg != edgb);

    return gt.normalize(n);
  }
}

} // namespace Operators
} // namespace FEVV
