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
  // for FEVV::CGALPointSet
#include "FEVV/Wrappings/Graph_traits_cgal_point_set.h"
  // for boost::graph_traits< FEVV::CGALPointSet >
#include "FEVV/Wrappings/Geometry_traits_cgal_point_set.h"
  // for FEVV::RetrieveKernel< FEVV::CGALPointSet >
#include "FEVV/Wrappings/Graph_properties_cgal_point_set.h"
  // for get(FEVV::CGALPointSetPointMap&, ...)
#include "FEVV/Wrappings/properties_cgal_point_set.h"
  // for FEVV::PMap_traits< FEVV::CGALPointSet > and CGAPPointSet-property-maps

#include "FEVV/Filters/CGAL/Point_set/cgal_point_set_reader.hpp"
  // for FEVV::Filters::read_mesh< FEVV::CGALPointSet >
#include "FEVV/Filters/CGAL/Point_set/cgal_point_set_writer.hpp"
  // for FEVV::Filters::write_mesh< FEVV::CGALPointSet >

#include "Examples/Generic/example_property_maps.hpp"


// main
int main(int argc, char *argv[])
{
  return example_property_maps< FEVV::CGALPointSet >(argc, argv);
}
