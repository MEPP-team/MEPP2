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

#include "boolops_polyhedra.hpp"

namespace FEVV {
namespace Filters {

//--------------------- UNION -------------------------

/**
 * \brief  Computes the union of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 * \param  gA      1st input mesh
 * \param  gB      2nd input mesh
 * \param  g_out   output mesh
 * \param  gt      the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 * 
 * \note   For compliance with the general project policy, the point maps
 *         of the input and output meshes should be passed as parameters.
 *         But the specific implementation of this filter can not use these
 *         property maps directly. In order to avoid passing unused 
 *         parameters, which later triggers compilation warning, it is 
 *         decided to derogate to the general project policy.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_union(HalfedgeGraph &gA,
              HalfedgeGraph &gB,
              HalfedgeGraph &g_out,
              const GeometryTraits &gt)
{
  BoolPolyhedra< HalfedgeGraph >(gA, gB, g_out, UNION);
}

/**
 * \brief  Computes the union of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 *         Use the default geometry traits of the mesh.
 *
 * \param  gA      1st input mesh
 * \param  gB      2nd input mesh
 * \param  g_out   output mesh
 *
 * \sa     the variant that use the geometry traits provided by the user.
 * 
 * \note   For compliance with the general project policy, the point maps
 *         of the input and output meshes should be passed as parameters.
 *         But the specific implementation of this filter can not use these
 *         property maps directly. In order to avoid passing unused 
 *         parameters, which later triggers compilation warning, it is 
 *         decided to derogate to the general project policy.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_union(HalfedgeGraph &gA,
              HalfedgeGraph &gB,
              HalfedgeGraph &g_out)
{
  GeometryTraits gt(gA);
  boolean_union< HalfedgeGraph, GeometryTraits >(gA, gB, g_out, gt);
}


//--------------------- INTERSECTION -------------------------

/**
 * \brief  Computes the intersection of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 * \param  gA      1st input mesh
 * \param  gB      2nd input mesh
 * \param  g_out   output mesh
 * \param  gt      the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 * 
 * \note   For compliance with the general project policy, the point maps
 *         of the input and output meshes should be passed as parameters.
 *         But the specific implementation of this filter can not use these
 *         property maps directly. In order to avoid passing unused 
 *         parameters, which later triggers compilation warning, it is 
 *         decided to derogate to the general project policy.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_inter(HalfedgeGraph &gA,
              HalfedgeGraph &gB,
              HalfedgeGraph &g_out,
              const GeometryTraits &gt)
{
  BoolPolyhedra< HalfedgeGraph >(gA, gB, g_out, INTER);
}

/**
 * \brief  Computes the intersection of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 *         Use the default geometry traits of the mesh.
 *
 * \param  gA      1st input mesh
 * \param  gB      2nd input mesh
 * \param  g_out   output mesh
 *
 * \sa     the variant that use the geometry traits provided by the user.
 * 
 * \note   For compliance with the general project policy, the point maps
 *         of the input and output meshes should be passed as parameters.
 *         But the specific implementation of this filter can not use these
 *         property maps directly. In order to avoid passing unused 
 *         parameters, which later triggers compilation warning, it is 
 *         decided to derogate to the general project policy.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_inter(HalfedgeGraph &gA,
              HalfedgeGraph &gB,
              HalfedgeGraph &g_out)
{
  GeometryTraits gt(gA);
  boolean_inter< HalfedgeGraph, GeometryTraits >(gA, gB, g_out, gt);
}


//--------------------- SUBTRACTION -------------------------

/**
 * \brief  Computes the subtraction of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 * \param  gA      1st input mesh
 * \param  gB      2nd input mesh
 * \param  g_out   output mesh
 * \param  gt      the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 * 
 * \note   For compliance with the general project policy, the point maps
 *         of the input and output meshes should be passed as parameters.
 *         But the specific implementation of this filter can not use these
 *         property maps directly. In order to avoid passing unused 
 *         parameters, which later triggers compilation warning, it is 
 *         decided to derogate to the general project policy.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_minus(HalfedgeGraph &gA,
              HalfedgeGraph &gB,
              HalfedgeGraph &g_out,
              const GeometryTraits &gt)
{
  BoolPolyhedra< HalfedgeGraph >(gA, gB, g_out, MINUS);
}

/**
 * \brief  Computes the subtraction of two polyhedra.
 *
 *         Ref: "Exact and Efficient Booleans for Polyhedra", C. Leconte,
 *              H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018,
 *              2010
 *
 *         Use the default geometry traits of the mesh.
 *
 * \param  gA      1st input mesh
 * \param  gB      2nd input mesh
 * \param  g_out   output mesh
 *
 * \sa     the variant that use the geometry traits provided by the user.
 * 
 * \note   For compliance with the general project policy, the point maps
 *         of the input and output meshes should be passed as parameters.
 *         But the specific implementation of this filter can not use these
 *         property maps directly. In order to avoid passing unused 
 *         parameters, which later triggers compilation warning, it is 
 *         decided to derogate to the general project policy.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
boolean_minus(HalfedgeGraph &gA,
              HalfedgeGraph &gB,
              HalfedgeGraph &g_out)
{
  GeometryTraits gt(gA);
  boolean_minus< HalfedgeGraph, GeometryTraits >(gA, gB, g_out, gt);
}


} // namespace Filters
} // namespace FEVV

