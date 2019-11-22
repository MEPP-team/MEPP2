// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

#include "FEVV/Wrappings/properties.h"
#include "FEVV/Wrappings/Geometry_traits.h"

#include <sstream> // std::ostringstream
// must be the last include to avoid side effect of disabling NDEBUG
#undef NDEBUG // enable assert in release mode
#include <cassert>

//----------------------------------------------------------

template< typename PointCloud,
          typename Point,
          typename Normal,
          typename Color >
void
test_property_maps_point_cloud(PointCloud& pc)
{
  FEVV::Geometry_traits< PointCloud > gt(pc);

  // create property maps bag
  FEVV::PMapsContainer pmaps_bag;

  // add some points
  for(int i = 0; i < 4; i++)
    add_vertex(pc);

  //--------------------------------------------------------------------------
  // test Point type
  //--------------------------------------------------------------------------
  //
  {
    Point p1; // ctor with no parameter
    Point p2(1., 2., 3.); // ctor with 3 parameters
    auto tmp = p1[0]; // operator[] as RValue
    std::ostringstream ss;
    ss << p2; // operator<<
  }

  //--------------------------------------------------------------------------
  // initialize PointMap
  //--------------------------------------------------------------------------
  //
  {
    auto pm = get(boost::vertex_point, pc);

    auto iterator_pair = vertices(pc);
    int cnt = 1;
    auto vi = iterator_pair.first;
    auto vi_end = iterator_pair.second;
    for(; vi != vi_end; ++vi)
    {
      put(pm, *vi, Point(cnt, cnt+1, cnt+2));
      cnt += 3;
    }

    assert(get(pm, 0) == Point(1, 2, 3));
    assert(get(pm, 1) == Point(4, 5, 6));
    assert(get(pm, 2) == Point(7, 8, 9));
    assert(get(pm, 3) == Point(10, 11, 12));
  }

  //--------------------------------------------------------------------------
  // read-write PointMap
  //--------------------------------------------------------------------------
  //
  {
    auto pm = get(boost::vertex_point, pc);

    auto iterator_pair = vertices(pc);
    auto vi = iterator_pair.first;
    auto vi_end = iterator_pair.second;
    for(; vi != vi_end; ++vi)
    {
      auto p = get(pm, *vi);
      Point p2(gt.get_x(p)*2, gt.get_y(p)*2, gt.get_z(p)*2);
      put(pm, *vi, p2);
    }

    assert(get(pm, 0) == Point(2, 4, 6));
    assert(get(pm, 1) == Point(8, 10, 12));
    assert(get(pm, 2) == Point(14, 16, 18));
    assert(get(pm, 3) == Point(20, 22, 24));
  }

  //--------------------------------------------------------------------------
  // test Normal type
  //--------------------------------------------------------------------------
  //
  {
    Normal n1; // ctor with no parameter
    Normal n2(1., 2., 3.); // ctor with 3 parameters
    auto tmp = n1[0]; // operator[] as RValue
    std::ostringstream ss;
    ss << n2; // operator<<
  }

  //--------------------------------------------------------------------------
  // create NormalMap
  //--------------------------------------------------------------------------
  //
  {
    auto v_nm = make_property_map(FEVV::vertex_normal, pc);
    put_property_map(FEVV::vertex_normal, pc, pmaps_bag, v_nm);
  }

  //--------------------------------------------------------------------------
  // initialize NormalMap
  //--------------------------------------------------------------------------
  //
  {
    auto v_nm = get_property_map(FEVV::vertex_normal, pc, pmaps_bag);

    auto iterator_pair = vertices(pc);
    int cnt = 10;
    auto vi = iterator_pair.first;
    auto vi_end = iterator_pair.second;
    for(; vi != vi_end; ++vi)
    {
      put(v_nm, *vi, Normal(cnt, cnt+1, cnt+2));
      cnt += 3;
    }

    assert(get(v_nm, 0) == Normal(10, 11, 12));
    assert(get(v_nm, 1) == Normal(13, 14, 15));
    assert(get(v_nm, 2) == Normal(16, 17, 18));
    assert(get(v_nm, 3) == Normal(19, 20, 21));
  }

  //--------------------------------------------------------------------------
  // read-write NormalMap
  //--------------------------------------------------------------------------
  //
  {
    auto v_nm = get_property_map(FEVV::vertex_normal, pc, pmaps_bag);

    auto iterator_pair = vertices(pc);
    auto vi = iterator_pair.first;
    auto vi_end = iterator_pair.second;
    for(; vi != vi_end; ++vi)
    {
      auto n = get(v_nm, *vi);
      Normal n2(n[0]*10, n[1]*10, n[2]*10);
      put(v_nm, *vi, n2);
    }

    assert(get(v_nm, 0) == Normal(100, 110, 120));
    assert(get(v_nm, 1) == Normal(130, 140, 150));
    assert(get(v_nm, 2) == Normal(160, 170, 180));
    assert(get(v_nm, 3) == Normal(190, 200, 210));
  }

  //--------------------------------------------------------------------------
  // test Color type
  //--------------------------------------------------------------------------
  //
  {
    Color c1; // ctor with no parameter
    Color c2(1., 2., 3.); // ctor with 3 parameters
    auto tmp = c1[0]; // operator[] as RValue
    std::ostringstream ss;
    ss << c2; // operator<<
  }

  //--------------------------------------------------------------------------
  // create ColorMap
  //--------------------------------------------------------------------------
  //
  {
    auto v_cm = make_property_map(FEVV::vertex_color, pc);
    put_property_map(FEVV::vertex_color, pc, pmaps_bag, v_cm);
  }

  //--------------------------------------------------------------------------
  // initialize ColorMap
  //--------------------------------------------------------------------------
  //
  {
    auto v_cm = get_property_map(FEVV::vertex_color, pc, pmaps_bag);

    auto iterator_pair = vertices(pc);
    int cnt = 210;
    auto vi = iterator_pair.first;
    auto vi_end = iterator_pair.second;
    for(; vi != vi_end; ++vi)
    {
      put(v_cm, *vi, Color(cnt, cnt+1, cnt+2));
      cnt += 3;
    }

    assert(get(v_cm, 0) == Color(210, 211, 212));
    assert(get(v_cm, 1) == Color(213, 214, 215));
    assert(get(v_cm, 2) == Color(216, 217, 218));
    assert(get(v_cm, 3) == Color(219, 220, 221));
  }

  //--------------------------------------------------------------------------
  // read-write ColorMap
  //--------------------------------------------------------------------------
  //
  {
    auto v_cm = get_property_map(FEVV::vertex_color, pc, pmaps_bag);

    auto iterator_pair = vertices(pc);
    auto vi = iterator_pair.first;
    auto vi_end = iterator_pair.second;
    for(; vi != vi_end; ++vi)
    {
      auto c = get(v_cm, *vi);
      Color c2(c[0]/10, c[1]/10, c[2]/10);
      put(v_cm, *vi, c2);
    }

    assert(get(v_cm, 0) == Color(21.0, 21.1, 21.2));
    assert(get(v_cm, 1) == Color(21.3, 21.4, 21.5));
    assert(get(v_cm, 2) == Color(21.6, 21.7, 21.8));
    assert(get(v_cm, 3) == Color(21.9, 22.0, 22.1));
  }

  //--------------------------------------------------------------------------
  // create non-standard property map
  //--------------------------------------------------------------------------
  //
  auto special_prop_map =
      FEVV::make_vertex_property_map< PointCloud, int >(pc);

  //--------------------------------------------------------------------------
  // initialize non-standard property map
  //--------------------------------------------------------------------------
  //
  {
    auto iterator_pair = vertices(pc);
    int cnt = 333;
    auto vi = iterator_pair.first;
    auto vi_end = iterator_pair.second;
    for(; vi != vi_end; ++vi)
    {
      put(special_prop_map, *vi, cnt);
      cnt++;
    }

    assert(get(special_prop_map, 0) == 333);
    assert(get(special_prop_map, 1) == 334);
    assert(get(special_prop_map, 2) == 335);
    assert(get(special_prop_map, 3) == 336);
  }

  //--------------------------------------------------------------------------
  // read-write non-standard property map
  //--------------------------------------------------------------------------
  //
  {
    auto iterator_pair = vertices(pc);
    auto vi = iterator_pair.first;
    auto vi_end = iterator_pair.second;
    for(; vi != vi_end; ++vi)
    {
      auto v = get(special_prop_map, *vi);
      put(special_prop_map, *vi, v*2);
    }

    assert(get(special_prop_map, 0) == 666);
    assert(get(special_prop_map, 1) == 668);
    assert(get(special_prop_map, 2) == 670);
    assert(get(special_prop_map, 3) == 672);
  }

  //--------------------------------------------------------------------------
  // all tests passed
  //--------------------------------------------------------------------------
  //
  std::cout << "Test passed.\n";
}
