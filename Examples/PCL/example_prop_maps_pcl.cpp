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
  // for FEVV::PCLPointCloud
#include "FEVV/Wrappings/Graph_traits_pcl_point_cloud.h"
  // for boost::graph_traits< FEVV::PCLPointCloud >
#include "FEVV/Wrappings/Geometry_traits_pcl_point_cloud.h"
  // for FEVV::RetrieveKernel< FEVV::PCLPointCloud >
#include "FEVV/Wrappings/Graph_properties_pcl_point_cloud.h"
  // for get(FEVV::PCLPointCloudPointMap&, ...)
#include "FEVV/Wrappings/properties_pcl_point_cloud.h"
  // for FEVV::PMap_traits< FEVV::PCLPointCloud >

#include "FEVV/Filters/PCL/pcl_point_cloud_reader.hpp"
  // for FEVV::Filters::read_mesh< FEVV::PCLPointCloud >
#include "FEVV/Filters/PCL/pcl_point_cloud_writer.hpp"
  // for FEVV::Filters::write_mesh< FEVV::PCLPointCloud >

#include "Examples/Generic/example_property_maps.hpp"


// main
int main(int argc, char *argv[])
{
  return example_property_maps< FEVV::PCLPointCloud >(argc, argv);
}
