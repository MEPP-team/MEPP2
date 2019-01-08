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
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

namespace FEVV {

using CGALKernel = CGAL::Cartesian< double >;
using CGALPoint = CGALKernel::Point_3;

using MeshSurface = CGAL::Surface_mesh< CGALPoint >;

} // namespace FEVV

