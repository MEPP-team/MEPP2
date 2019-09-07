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
#include "FEVV/Operators/Generic/calculate_face_normal.hpp"

#ifndef NDEBUG
#include <cassert>
#endif

namespace FEVV {
namespace Operators {

/**
 * \brief Computes the unit normal at considered vertex v as the average of
 *        the normals of incident faces.
 * \note The normals of the adjacent faces are computed on the fly and droped
 *       on function exit.
 * \sa   When already having a precomputed Property Map of the normals
 *       consider using the other version of the function.
 * \param[in] v     The considered vertex
 * \param[in] g     The mesh to which the vertex belongs
 * \param[in] pm    The Point Map holding the vertices geometry
 * \param[in] gt    The Geometry Traits associated to the mesh
 * \return   The vertex v' normal.
 */
template< typename HalfedgeGraph, typename PointMap, typename GeometryTraits >
typename GeometryTraits::Vector
calculate_vertex_normal(
    typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor v,
    const HalfedgeGraph &g,
    const PointMap &pm,
    const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor
      halfedge_descriptor;

  if(halfedge(v, g) == boost::graph_traits< HalfedgeGraph >::null_halfedge())
  {
    return gt.NULL_VECTOR; // instead throws an exception?
  }

  Vector normal = gt.NULL_VECTOR;
#ifndef NDEBUG
  int cpt = 0;
  int nbBorders = 0;
#endif
try{
  halfedge_descriptor he = halfedge(v, g);
  halfedge_descriptor end = he;
  do
  {
    if(!CGAL::is_border(he, g))
    {
      Vector n = FEVV::Operators::
          calculate_face_normal< HalfedgeGraph, PointMap, GeometryTraits >(
              face(he, g), g, pm, gt);
      normal = gt.add(normal, n);
#ifndef NDEBUG
      ++cpt;
#endif
    }
#ifndef NDEBUG
    else
      ++nbBorders;
#endif
    he = opposite(next(he, g), g);
  } while(he != end);
#ifndef NDEBUG
  assert(cpt == (static_cast< int >(degree(v, g)) - nbBorders));
#endif
}
catch(...)
{
	return gt.NULL_VECTOR;
}
  return gt.normalize(normal);
}

/**
 * \brief Computes the unit normal at considered vertex v as the average of
 *        the precomputed normals (provided as argument) of incident faces.
 * \sa   When the Property Map of the normals is not at hand consider using
 *       the other version of the function.
 * \param[in] v     The considered vertex
 * \param[in] g     The mesh to which the vertex belongs
 * \param[in] pm    The Point Map holding the vertices geometry
 * \param[in] fnm   The Face Normal Map providing the normals to incident faces
 * \param[in] gt    The Geometry Traits associated to the mesh
 * \return   The vertex v' normal.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename FaceNormalMap,
          typename GeometryTraits >
typename GeometryTraits::Vector
calculate_vertex_normal(
    typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor v,
    const HalfedgeGraph &g,
    const PointMap &pm,
    const FaceNormalMap &fnm, /// provided face normals
    const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Vector Vector;
  typedef typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor
      halfedge_descriptor;

  if(halfedge(v, g) ==
     boost::graph_traits< HalfedgeGraph >::null_halfedge()) // if ( degree(v, g)
                                                            // <= 1 )
    return gt.NULL_VECTOR; // instead throws an exception?

  Vector normal = gt.NULL_VECTOR;
#ifndef NDEBUG
  int cpt = 0;
  int nbBorders = 0;
#endif
try{
  halfedge_descriptor he = halfedge(v, g);
  halfedge_descriptor end = he;
  do
  {
    if(!CGAL::is_border(he, g))
    {
      normal = gt.add(normal, get(fnm, face(he, g)));
#ifndef NDEBUG
      ++cpt;
#endif
    }
#ifndef NDEBUG
    else
      ++nbBorders;
#endif
    he = opposite(next(he, g), g);
  } while(he != end);
#ifndef NDEBUG
  assert(cpt == (static_cast< int >(degree(v, g)) - nbBorders));
#endif
}
catch(...)
{
	return gt.NULL_VECTOR;
}
  return gt.normalize(normal);
}

} // namespace Operators
} // namespace FEVV
