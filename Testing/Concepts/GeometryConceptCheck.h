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
struct Point_3Concept
//  : boost::Assignable<X>
{
  typedef Geometry_traits< Mesh > Geometry;
  typedef Geometry G; // Syntactic shortcut for operators
  typedef typename Geometry::Point Point;
  typedef typename Geometry::Vector Vector;
  typedef typename Geometry::Scalar Scalar;

  BOOST_CONCEPT_USAGE(Point_3Concept)
  {
    Mesh mesh;
    Geometry g(mesh);


    //////////////////////// Point manipulations

    // Point copy construction
    Point q(p);

    // Point copy assignement operator
    q = p;

    // Point coordinate component access
    Scalar x = g.get_x(p);
    Scalar y = g.get_y(p);
    Scalar z = g.get_z(p);

    //////////////////////// Vector manipulations
    // Vector copy construction
    Vector v2(v1);

    // Vector copy assignement:
    Vector v1 = Vector(1, 2, 3);

    // Vector/Scalar coordinate component access
    x = v1[0];
    y = v1[1];
    z = v1[2];

    // One can NOT overwrite a single vector coordinate component e.g.
    // v[0] = 1; will fail (refer below for the details).
    // In order to alter one or many coordinate component of a Vector one must
    // create a new Vector and reassign it:
    v1 = Vector(v1[0], v1[1], 999.0);

    // Vector/Vector addition is generally provided natively by the
    // implementations
    Vector v3 = v1 + v2;
    Vector v4 = v1 - v2;

    g.dot_product(v1, v2);
    g.cross_product(v1, v2);

    // Scalar/Vector:
    s = g.length(v1);
    v2 = g.scalar_mult(v1, s);

    //////////////////////// Vector/Point and Point/Vector operators
    v1 = g.normal(p, p, p);
    v1 = g.unit_normal(p, p, p);

    v3 = g.add(v1, v2);
    p = g.add_p(p, v1);
    p = g.sub_p(p, v1);
    v1 = g.sub(p, p);


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

    // Deprecated, not worth/nice using
    x = g.template get< 0 >(p); // deprecated
    y = g.template get< 1 >(p); // deprecated
    z = g.template get< 2 >(p); // deprecated

/////////////////////////////////////////////////////////////
// Only working for some data structures and thus not generic
#ifndef FEVV_USE_CGAL
// Point r = p + q;
#endif
  }

private:
  Point p;
  Vector v1;
  Scalar s;
};

} // namespace FEVV

