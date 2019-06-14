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
#include <vector>

namespace FEVV {

using CGALKernel = CGAL::Cartesian< double >;
using CGALPoint = CGALKernel::Point_3;

using CGALPointSet = std::vector< CGALPoint >;
//Note: could be a std::list too

//TODO-elo: generalization
//template< typename PointT >
//using CGALVectorPointSet = std::vector< PointT >;
//using CGALPointSet = CGALVectorPointSet< CGALPoint >;

//TODO-elo: later point+normal
//  using CGALVector = CGALKernel::Vector_3;
//  using CGALPointVectorPair = std::pair< CGALPoint, CGALVector >;
//  using CGALPointSet = std::vector< CGALPointVectorPair >;

} // namespace FEVV

