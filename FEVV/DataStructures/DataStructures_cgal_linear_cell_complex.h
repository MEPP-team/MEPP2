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

#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

// FIXME: the following inclusion can be removed once the generic
//        reader becomes the only one used (as opposed to using the
//        LCC native ones).
#include <CGAL/Linear_cell_complex_incremental_builder.h>

namespace FEVV {

using CGALKernel = CGAL::Cartesian< double >;
using CGALLCCTraits = CGAL::Linear_cell_complex_traits< 3, CGALKernel >;

struct CGALItem
{
  template< class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Tag_true Darts_with_id;
    using Vertex_attribute = CGAL::Cell_attribute_with_point_and_id< Refs >;
    using Face_attribute = CGAL::Cell_attribute_with_id< Refs >;
    using Attributes =
        CGAL::cpp11::tuple< Vertex_attribute, void, Face_attribute >;
  };
};

using MeshLCC = CGAL::
    Linear_cell_complex_for_combinatorial_map< 2, 3, CGALLCCTraits, CGALItem >;

} // namespace FEVV

