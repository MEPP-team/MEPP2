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
#pragma once

#if defined _MSC_VER
#pragma warning(disable : 4244)
  // workaround issue https://github.com/CGAL/cgal/issues/4367
  // remove the pragma when issue is fixed
#endif
#include <CGAL/version_macros.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/IO/Color.h>


namespace FEVV {


using CGALPointSetKernel = CGAL::Cartesian< float >;
using CGALPointSetPoint  = CGALPointSetKernel::Point_3;
using CGALPointSetVector = CGALPointSetKernel::Vector_3;

// note:
// - CGAL::Color in CGAL 4.14 misses operator[] that is used in MEPP2
//   so we define our own color type derived from CGAL::Color
// - starting from CGAL-5.0-beta1, CGAL::Color has everything needed
//   by MEPP2 and we must then use it directly
//TODO-elo-CGAL-5.0-beta1  using CGALPointSetColor = CGAL::Color;
class CGALPointSetColor : public CGAL::Color
{
public:
  // must redefine default constructor because
  // constructor(r,g,b) is re-defined
  CGALPointSetColor(void) : CGAL::Color()
  {}

  // must redefine constructor(r,g,b) to use
  // base class constructor(r,g,b)
  CGALPointSetColor(unsigned char r, unsigned char g, unsigned char b)
      : CGAL::Color(r, g, b)
  {}

  // color must behave as a vector
  unsigned char operator[](int i)	const
  {
    if(i == 0)
      return r();
    else if(i == 1)
      return g();
    else if(i == 2)
      return b();
    else
    {
      throw std::runtime_error(
          "FEVV::CGALPointSetColor: invalid index in operator[]");
      return 0; // dummy return to avoid a warning
    }
  }
};

using CGALPointSet = CGAL::Point_set_3< CGALPointSetPoint >;


} // namespace FEVV

