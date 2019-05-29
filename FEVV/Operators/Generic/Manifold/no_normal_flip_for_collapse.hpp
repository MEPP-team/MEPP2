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

#include <CGAL/boost/graph/helpers.h> // for is_border()

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Tools/Math/MatrixOperations.hpp"

namespace FEVV {
namespace Operators {
/*!
 * \brief Function used to detected a local normal flip
 *        due to a vertex position change.
 *
 * \tparam  HalfedgeGraph a Mesh type that provides a Model of the
 *          HalfedgeGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \tparam  PointMap A modifiable point map to manage vertex positions.  
 * \tparam GeometryTraits The geometric kernel when available. This is defaulted
 *         to FEVV::Geometry_traits<HalfedgeGraph>.
 * \param[in] g The HalfedgeGraph instance.
 * \param[in] pm The point map which associates a vertex descriptor with a vertex
 *          position.
 * \param[in] e A candidate edge to collapse.
 * \param[in] v The vertex incident to edge e whose position is modified. 
 * \param[in] new_pos_v The new position for vertex v.
 * \param[in] gt The geometry trait object. 
 * \return  true when no normal flip will occur by changing the position of vertex
            v to new_pos_v.
 */	
template< typename HalfedgeGraph, 
          typename PointMap, 
          typename GeometryTraits = typename FEVV::Geometry_traits<HalfedgeGraph> >
  bool no_normal_flip_for_collapse(
      const HalfedgeGraph& g,
      const PointMap& pm,
      typename boost::graph_traits<HalfedgeGraph>::edge_descriptor e,
      typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor v,
      const typename GeometryTraits::Point& new_pos_v,
      const GeometryTraits &gt)
  {
    typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h1 = halfedge(e, g);
    typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h2 = opposite(h1, g);
    auto iterator_range = CGAL::halfedges_around_target(v, g);
    for (auto h : iterator_range)
    {
      if (CGAL::is_border(h, g))
        continue;

      // cases of the incident faces to the collapsed edge
      typename boost::graph_traits<HalfedgeGraph>::face_descriptor f = face(h, g); // we know that f!=null_face()
      if (f == face(h1, g))
        continue;
      if (f == face(h2, g))
        continue;

      // other faces
      typename GeometryTraits::Point p1 = get(pm, source(h, g));
      typename GeometryTraits::Point p0 = get(pm, source(prev(h, g), g));

      typename GeometryTraits::Vector face_normal, face_normal_after_collapse;
      // get the normal of the current face
      // check if the current face is not null
      bool b = FEVV::Math::Vector::are_collinear<GeometryTraits>(gt.sub(p0, p1), gt.sub(get(pm, v), p1));
      if (!b)
        face_normal = gt.normal(p0, p1, get(pm, v));
      else
      {
        face_normal = typename GeometryTraits::Vector(std::numeric_limits<typename GeometryTraits::Scalar>::quiet_NaN(),
          std::numeric_limits<typename GeometryTraits::Scalar>::quiet_NaN(),
          std::numeric_limits<typename GeometryTraits::Scalar>::quiet_NaN());
      }
      // simulation of the new face after collapse and check if its normal has the same sign
      face_normal_after_collapse = gt.normal(p0, p1, new_pos_v);

      const typename  GeometryTraits::Scalar sign = gt.dot_product(face_normal, face_normal_after_collapse);
      if (sign < 0.0)
        return false;
    }
    return true;
  }

} // namespace Operators
} // namespace FEVV
