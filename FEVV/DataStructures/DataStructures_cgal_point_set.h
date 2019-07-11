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

#include <CGAL/version_macros.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_set_3.h>


namespace FEVV {


using CGALPointSetKernel = CGAL::Cartesian< double >;
using CGALPointSetPoint  = CGALPointSetKernel::Point_3;
using CGALPointSetvector = CGALPointSetKernel::Vector_3;

using CGALPointSet = CGAL::Point_set_3< CGALPointSetPoint >;


} // namespace FEVV

