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
#include <iostream>

#if defined _MSC_VER
#pragma warning(disable : 4996) // MT
#endif

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

int
main(int argc, char **argv)
{
  pcl::PointCloud< pcl::PointXYZ > cloud;

  // Fill in the cloud data
  cloud.width = 5;
  cloud.height = 1;
  cloud.is_dense = false;
  cloud.points.resize(cloud.width * cloud.height);

  for(size_t i = 0; i < cloud.points.size(); ++i)
  {
    cloud.points[i].x = 1024 * rand() / (RAND_MAX + 1.0f);
    cloud.points[i].y = 1024 * rand() / (RAND_MAX + 1.0f);
    cloud.points[i].z = 1024 * rand() / (RAND_MAX + 1.0f);
  }

  pcl::io::savePCDFileASCII("test_pcd.pcd", cloud);
  std::cerr << "Saved " << cloud.points.size()
            << " data points to test_pcd.pcd." << std::endl;

  for(size_t i = 0; i < cloud.points.size(); ++i)
    std::cerr << "    " << cloud.points[i].x << " " << cloud.points[i].y << " "
              << cloud.points[i].z << std::endl;

  return (0);
}
