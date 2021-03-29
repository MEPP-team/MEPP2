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

#undef NDEBUG // enable assert in release build
#include <cassert>
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
    Point p(1.1f, 2.2f, 3.3f);

    // Point copy construction
    Point q(p);

    // Point assignement operator
    q = p;

    // Point at origin
    p = G::ORIGIN;
    p = gt.ORIGIN;

    // Point coordinate component access
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
    Scalar x = gt.get_x(p);
    Scalar y = gt.get_y(p);
    Scalar z = gt.get_z(p);
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

    //////////////////////// Vector manipulations

    // Vector constructor
    Vector v(4.4f, 5.5f, 6.6f);

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
    v1 = Vector(v1[0], v1[1], 999.0f);

    // Misc vectors operations
    v1 = gt.normalize(v);
    s  = gt.length2(v);
    s  = gt.length(v);
    v = gt.add_v(v1, v2);
    gt.dot_product(v1, v2);
    gt.cross_product(v1, v2);

    // Scalar/Vector:
    v2 = gt.scalar_mult(v1, s);

    //////////////////////// Vector/Point and Point/Vector operators

    s  = gt.length(p, q);
    v1 = gt.normal(p, p, p);
    v1 = gt.unit_normal(p, p, p);
    p = gt.add_pv(p, v1);
    p = gt.sub_pv(p, v1);
    v1 = gt.sub_p(p, p);

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


/**
 * Check operators results.
 */
template< typename GeometryTraits >
void check_operators_results(void)
{
  typedef typename GeometryTraits::Point  Point;
  typedef typename GeometryTraits::Vector Vector;
  typedef typename GeometryTraits::Scalar Scalar;

  // helper lambda function to test if two FP numbers are almost equal
  auto almost_equal_s = [](Scalar s1, Scalar s2, Scalar epsilon) -> bool
  {
    return std::abs(s1 - s2) <= epsilon;
  };

  // helper lambda function to test if two FP points are almost equal
  auto almost_equal_p = [almost_equal_s](Point p1, Point p2, Scalar epsilon) -> bool
  {
    return almost_equal_s(GeometryTraits::get_x(p1),
                          GeometryTraits::get_x(p2),
                          epsilon) &&
           almost_equal_s(GeometryTraits::get_y(p1),
                          GeometryTraits::get_y(p2),
                          epsilon) &&
           almost_equal_s(GeometryTraits::get_z(p1),
                          GeometryTraits::get_z(p2),
                          epsilon);
  };

  // helper lambda function to test if two FP vectors are almost equal
  auto almost_equal_v = [almost_equal_s](Vector v1, Vector v2, Scalar epsilon) -> bool
  {
    return almost_equal_s(v1[0], v2[0], epsilon) &&
           almost_equal_s(v1[1], v2[1], epsilon) &&
           almost_equal_s(v1[2], v2[2], epsilon);
  };

  Scalar eps = 1e-5f;
  Vector u(8.3f, -3.2f, -3.8f);
  Vector v(1.2f, 2.3f, 3.4f);
  Point p(3.0f, 12.5f, -2.1f);
  Point q(-1.4f, 4.6f, -1.9f);
  Point r(-2.6f, 4.7f, 5.8f);

  {
    Scalar l2 = GeometryTraits::length2(v);
    assert(almost_equal_s(l2, 18.29f, eps));

    Scalar l = GeometryTraits::length(v);
    assert(almost_equal_s(l, 4.276681f, eps));

    Vector n = GeometryTraits::normalize(v);
    assert(almost_equal_v(n, Vector(0.280591f, 0.537800f, 0.795009f), eps));
  }

  {
    Vector w = GeometryTraits::sub_p(p, q);
    assert(almost_equal_v(w, Vector(4.4f, 7.9f, -0.2f), eps));

    Scalar l = GeometryTraits::length(p, q);
    assert(almost_equal_s(l, 9.044888f, eps));
  }

  {
    // these two function are not implemented with static method
    // with OpenMesh, so we must create a mesh object and a geometry traits
    // object
    typename GeometryTraits::Mesh m;
    GeometryTraits gt(m);

    Vector n = gt.normal(p, q, r);
    assert(almost_equal_v(n, Vector(-60.85f, 33.64f, -9.92f), eps) ||
           almost_equal_v(n, Vector(-0.866393f, 0.478972f, -0.141243f), eps) );
           // OpenMesh returns a normalized normal
    
    n = gt.unit_normal(p, q, r);
    assert(almost_equal_v(n, Vector(-0.866393f, 0.478972f, -0.141243f), eps));
  }

  {
    Vector w = GeometryTraits::add_v(u, v);
    assert(almost_equal_v(w, Vector(9.5f, -0.9f, -0.4f), eps));

    Scalar d = GeometryTraits::dot_product(u, v);
    assert(almost_equal_s(d, -10.32f, eps));

    w = GeometryTraits::cross_product(u, v);
    assert(almost_equal_v(w, Vector(-2.14f, -32.78f, 22.93f), eps));
  }

  {
    Point a = GeometryTraits::add_pv(p, v);
    assert(almost_equal_p(a, Point(4.2f, 14.8f, 1.3f), eps));

    a = GeometryTraits::sub_pv(p, v);
    assert(almost_equal_p(a, Point(1.8f, 10.2f, -5.5f), eps));

    Vector w = GeometryTraits::scalar_mult(v, -2.345f);
    assert(almost_equal_v(w, Vector(-2.814f, -5.3935f, -7.973f), eps));
  }
}

} // namespace FEVV

