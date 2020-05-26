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
#include "FEVV/DataStructures/DataStructures_cgal_point_set.h"

#include "FEVV/Wrappings/Graph_traits_extension_cgal_point_set.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_point_set.h"
#include "FEVV/Wrappings/properties_cgal_point_set.h"

#include "FEVV/Filters/CGAL/Point_set/cgal_point_set_reader.hpp"
#include "FEVV/Filters/CGAL/Point_set/cgal_point_set_writer.hpp"
#include "FEVV/Filters/CGAL/Point_set/copy_graph_cgal_point_set.hpp"

#include "Testing/Generic/test_copy_graph.inl"


int
main(int argc, const char **argv)
{
  return test_copy_graph< FEVV::CGALPointSet >(argc, argv);
}
