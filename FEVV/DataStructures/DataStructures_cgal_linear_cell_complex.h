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

#include <CGAL/version_macros.h>

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

