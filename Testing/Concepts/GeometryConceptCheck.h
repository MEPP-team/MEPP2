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

#include <boost/concept_check.hpp>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/Geometry_traits_operators.h"

namespace FEVV {

/**
 * \brief An elementary concept check \ref Geometry_traits specializations
 * [Boost Concept Check
 * Library](http://www.boost.org/doc/libs/1_60_0/libs/concept_check/concept_check.htm)
 * (BCCL)
 */
template< typename Mesh >
struct GeometryConcept
{
  typedef Geometry_traits< Mesh > Geometry;
  typedef Geometry G; // Syntactic shortcut for operators
  typedef typename Geometry::Point Point;
  typedef typename Geometry::Vector Vector;
  typedef typename Geometry::Scalar Scalar;

  BOOST_CONCEPT_USAGE(GeometryConcept)
  {
    Mesh mesh;
    Geometry gt(mesh);
    Vector v1;
    Scalar s;

    //////////////////////// Point manipulations

    // Point constructor
    Point p(1.1, 2.2, 3.3);

    // Point copy construction
    Point q(p);

    // Point assignement operator
    q = p;

    // Point at origin
    p = G::ORIGIN;
    p = gt.ORIGIN;

    // Point coordinate component access
    Scalar x = gt.get_x(p);
    Scalar y = gt.get_y(p);
    Scalar z = gt.get_z(p);

    //////////////////////// Vector manipulations

    // Vector constructor
    Vector v(4.4, 5.5, 6.6);

    // Vector copy construction
    Vector v2(v1);

    // Vector assignement operator:
    v2 = v1;

    // Null vector
    v = G::NULL_VECTOR;
    v = gt.NULL_VECTOR;

    // Vector/Scalar coordinate component access
    x = v1[0];
    y = v1[1];
    z = v1[2];

    // One can NOT overwrite a single vector coordinate component e.g.
    // v[0] = 1; will fail (refer below for the details).
    // In order to alter one or many coordinate component of a Vector one must
    // create a new Vector and reassign it:
    v1 = Vector(v1[0], v1[1], 999.0);

    // Misc vectors operations
    v1 = gt.normalize(v);
    s  = gt.length2(v);
    s  = gt.length(v);
    v = gt.add(v1, v2);
    gt.dot_product(v1, v2);
    gt.cross_product(v1, v2);

    // Scalar/Vector:
    v2 = gt.scalar_mult(v1, s);
    v2 = gt.scalar_mult(s, v1);

    //////////////////////// Vector/Point and Point/Vector operators

    s  = gt.length(p, q);
    v1 = gt.normal(p, p, p);
    v1 = gt.unit_normal(p, p, p);
    p = gt.add_p(p, v1);
    p = gt.sub_p(p, v1);
    v1 = gt.sub(p, p);

#if 0
////////////////////////////////////////////////////////////////
// The following is on debugging purposes, remains as traces
// of various trials, is deprecated or not recommendable.
// Do NOT use it in your generic code.

// SurfaceMesh fail: operator[] returns a const value
#ifndef FEVV_USE_CGAL
    v1[0] = 5;
#endif

    // Not recommendable because the notation is too heavy
    v1 = GeometryTraits::sub< G >(p, p);
    v1 = GeometryTraits::add_v< G >(v1, v2);
    p = GeometryTraits::add_p< G >(p, v1);

/////////////////////////////////////////////////////////////
// Only working for some data structures and thus not generic
#ifndef FEVV_USE_CGAL
// Point r = p + q;
#endif
#endif
  }
};

} // namespace FEVV

