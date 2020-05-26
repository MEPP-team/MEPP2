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
#include "FEVV/DataStructures/DataStructures_pcl_point_cloud.h"

#include "FEVV/Wrappings/Graph_traits_extension_pcl_point_cloud.h"
#include "FEVV/Wrappings/Geometry_traits_pcl_point_cloud.h"
#include "FEVV/Wrappings/properties_pcl_point_cloud.h"

#include "FEVV/Filters/PCL/pcl_point_cloud_reader.hpp"
#include "FEVV/Filters/PCL/pcl_point_cloud_writer.hpp"
#include "FEVV/Filters/PCL/copy_graph_pcl.hpp"

#include "Testing/Generic/test_copy_graph.inl"


int
main(int argc, const char **argv)
{
  return test_copy_graph< FEVV::PCLPointCloud >(argc, argv);
}
