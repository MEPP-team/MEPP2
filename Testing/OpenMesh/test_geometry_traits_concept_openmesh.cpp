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
/**
 * \todo
 * CLEANME
 * The following test is a quest to find the minimum geometrical predicates
 * required in order to write the SRC/Filters without CGAL.
 * For example FEVV/Filters/normal.h currently uses CGAL::unit_normal()
 * which code (see the copy at the bottom of this file) assumes that the
 * "Mesh::Point" coordinates can be accessed through members x, y and z
 * i.e. one can write "Mesh::Point P; and then "auto a = p.x;".
 * But OpenMesh coordinate naming scheme uses the vector notation and
 * coordinate access is made through p[Ø), p[1] and p[2].
 * This is precisely the type of difficulty that Boost.Geometry tries to
 * avoid. For example the
 * [point concept of
 * Boost.Geometry"](http://www.boost.org/doc/libs/1_60_0/libs/geometry/doc/html/geometry/reference/concepts/concept_point.html)
 * should probably be a source of inspiration. In particular for its
 * abstraction (through traits) or coordinate access.
 */

#include "FEVV/Wrappings/Geometry_traits_openmesh.h" // Always in first position for MVS 4996 warnings

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "Testing/Concepts/GeometryConceptCheck.h"


int
main(int narg, char **argv)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> Mesh;
  BOOST_CONCEPT_ASSERT((FEVV::Point_3Concept< Mesh >));

  typedef FEVV::Geometry_traits< Mesh > Geometry;
  Mesh mesh;
  Geometry g(mesh);
  FEVV::GeometryTraits::verify_assumption(g);
  FEVV::GeometryTraits::verify_operator_results(g);

  return 0;
}
