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

#include "FEVV/Filters/Generic/normalize_vector_map.h"

#include <iostream>
#include <cmath>

#undef NEAR
#define NEAR(val1, val2, epsilon)   (fabs((val1)-(val2)) <= (epsilon))

#undef NDEBUG // enable assert in release mode
              // must be the last include to avoid side effect
#include <cassert>


template< typename MeshT >
int
test_normalize_vector_map(void)
{
  typedef FEVV::Geometry_traits< MeshT >   GeometryTraits;
  typedef typename GeometryTraits::Vector  Vector;

  MeshT m;
  GeometryTraits gt(m);

  float eps = 1e-5f;

  // test vertex vector map normalization

  {
    auto v1 = add_vertex(m);
    auto v2 = add_vertex(m);

    auto vn_pm = make_property_map(FEVV::vertex_normal, m);
    put(vn_pm, v1, Vector(2.0f, 3.0f, 4.0f));
    put(vn_pm, v2, Vector(3.0f, 4.0f, 2.0f));

    FEVV::Filters::normalize_vector_map_vertices(m, vn_pm);

    auto vec1 = get(vn_pm, v1);
    std::cout << "v vec1 = " << vec1[0] << " " << vec1[1] << " "  << vec1[2]
              << std::endl;
    assert(NEAR(vec1[0], 0.371391f, eps));
    assert(NEAR(vec1[1], 0.557086f, eps));
    assert(NEAR(vec1[2], 0.742781f, eps));

    auto vec2 = get(vn_pm, v2);
    std::cout << "v vec2 = " << vec2[0] << " " << vec2[1] << " "  << vec2[2]
              << std::endl;
    assert(NEAR(vec2[0], 0.557086f, eps));
    assert(NEAR(vec2[1], 0.742781f, eps));
    assert(NEAR(vec2[2], 0.371391f, eps));
  }

  // test face vector map normalization

  {
    auto f1 = add_face(m);
    auto f2 = add_face(m);

    auto fn_pm = make_property_map(FEVV::face_normal, m);
    put(fn_pm, f1, Vector(4.0f, 2.0f, 3.0f));
    put(fn_pm, f2, Vector(2.0f, 3.0f, 4.0f));

    FEVV::Filters::normalize_vector_map_faces(m, fn_pm);

    auto vec1 = get(fn_pm, f1);
    std::cout << "f vec1 = " << vec1[0] << " " << vec1[1] << " "  << vec1[2]
              << std::endl;
    assert(NEAR(vec1[0], 0.742781f, eps));
    assert(NEAR(vec1[1], 0.371391f, eps));
    assert(NEAR(vec1[2], 0.557086f, eps));

    auto vec2 = get(fn_pm, f2);
    std::cout << "f vec2 = " << vec2[0] << " " << vec2[1] << " "  << vec2[2]
              << std::endl;
    assert(NEAR(vec2[0], 0.371391f, eps));
    assert(NEAR(vec2[1], 0.557086f, eps));
    assert(NEAR(vec2[2], 0.742781f, eps));
  }

  // test halfedge vector map normalization

  {
    auto e1 = add_edge(m);
    auto h1 = halfedge(e1, m);
    auto h2 = opposite(h1, m);

    auto hn_pm = make_property_map(FEVV::halfedge_normal, m);
    put(hn_pm, h1, Vector(3.0f, 4.0f, 2.0f));
    put(hn_pm, h2, Vector(4.0f, 2.0f, 3.0f));

    FEVV::Filters::normalize_vector_map_halfedges(m, hn_pm);

    auto vec1 = get(hn_pm, h1);
    std::cout << "h vec1 = " << vec1[0] << " " << vec1[1] << " "  << vec1[2]
              << std::endl;
    assert(NEAR(vec1[0], 0.557086f, eps));
    assert(NEAR(vec1[1], 0.742781f, eps));
    assert(NEAR(vec1[2], 0.371391f, eps));

    auto vec2 = get(hn_pm, h2);
    std::cout << "h vec2 = " << vec2[0] << " " << vec2[1] << " "  << vec2[2]
              << std::endl;
    assert(NEAR(vec2[0], 0.742781f, eps));
    assert(NEAR(vec2[1], 0.371391f, eps));
    assert(NEAR(vec2[2], 0.557086f, eps));
  }

  return 0;
}
