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

// IMPORTANT !!!
// NOTE : if we include <CGAL/version_macros.h> we have macro redefinition of
// CGAL_VERSION_STR
#ifndef CGAL_VERSION_MAJOR
#define CGAL_VERSION_MAJOR (CGAL_VERSION_NR / 10000000 % 100)
#endif

#ifndef CGAL_VERSION_MINOR
#define CGAL_VERSION_MINOR (CGAL_VERSION_NR / 100000 % 100)
#endif

#ifndef CGAL_VERSION_PATCH
#define CGAL_VERSION_PATCH (CGAL_VERSION_NR / 10000 % 10)
#endif
// IMPORTANT !!!

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

namespace FEVV {

using CGALKernel = CGAL::Cartesian< double >;
using MeshPolyhedron =
    CGAL::Polyhedron_3< CGALKernel, CGAL::Polyhedron_items_with_id_3 >;

} // namespace FEVV

