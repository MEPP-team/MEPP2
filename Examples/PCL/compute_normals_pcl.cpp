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

#if 1

// *** MEPP2 code calling PCL native filters.            ***
// *** Point and normal stored in a FEVV::PCLPointCloud. ***

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

#include <pcl/search/impl/search.hpp>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree_flann.h>


// main
int main(int argc, char *argv[])
{
  // parse arguments
  if(argc != 3)
  {
    std::cout << "Usage:  " << argv[0]
              << "  input_mesh_filename  output_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0] << "  ../Testing/Data/tetra.xyz  tetra.out.ply"
              << std::endl;
    return EXIT_FAILURE;
  }
  std::string input_file = argv[1];
  std::string output_file = argv[2];

  //----------------------------------

  // create point cloud
  typedef  FEVV::PCLPointCloud  PointCloudT;
  PointCloudT pc;

  // load point cloud
  std::cout << "load point cloud" << std::endl;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file, pc, pmaps_bag);

  //----------------------------------

  // retrieve Point Map
  auto pm = get(boost::vertex_point, pc);

  // create vertex-normal property map
  using VertexNormalMap = typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                                      PointCloudT >::pmap_type;
  VertexNormalMap v_nm;
  std::cout << "create vertex-normal map" << std::endl;
  v_nm = make_property_map(FEVV::vertex_normal, pc);

  // store property map in property maps bag
  put_property_map(FEVV::vertex_normal, pc, pmaps_bag, v_nm);

  // compute normals
  std::cout << "compute normals" << std::endl;

  //----------------------------------
  // PCL specific section - begin
  //----------------------------------

  typedef  FEVV::PCLEnrichedPoint  PointType;
  typedef  FEVV::PCLEnrichedPoint  NormalType;
    // use the cloud point to store the normal

  pcl::NormalEstimation< PointType, NormalType > normal_estimator;
  normal_estimator.setInputCloud(pc.makeShared());

  pcl::search::KdTree< PointType >::Ptr kdtree(
      new pcl::search::KdTree< PointType >);
  normal_estimator.setSearchMethod(kdtree);

  normal_estimator.setRadiusSearch(0.03);
  normal_estimator.compute(pc);

  //----------------------------------
  // PCL specific section - end
  //----------------------------------

  // save point cloud with normals
  std::cout << "save point cloud with normals" << std::endl;
  FEVV::Filters::write_mesh(output_file, pc, pmaps_bag);
  
  return 0;
}

#else

// *** PCL Native code. ***

/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2014-, Open Perception, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */

// partially inspired by
// - pcl-1.8.1-src/doc/tutorials/content/sources/
//       global_hypothesis_verification/global_hypothesis_verification.cpp
// - pcl-1.8.1-src/doc/tutorials/content/sources/
//       correspondence_grouping/correspondence_grouping.cpp

#include <iostream>

#include <pcl/search/impl/search.hpp>

#include <pcl/io/pcd_io.h>
#include <pcl/io/ascii_io.h>
#include <pcl/point_cloud.h>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <string>


// help
void
showHelp(char * program_name)
{
  std::cout << std::endl;
  std::cout << "Usage: " << program_name << " cloud_filename.xyz" << std::endl;
}


// print cloud
template< typename CloudType >
void
print_cloud(const CloudType &cloud)
{
  std::cout << "*** PointCloud = \n" << cloud << std::endl;

  std::cout << "*** points = " << std::endl;
  for(int i = 0; i < cloud.points.size(); i++)
  {
    std::cout << "cloud.points[" << i << "] = " << cloud.points[i]
              << std::endl;
  }

#if 0
  // how to access point and normal

  std::cout << "*** points.data = " << std::endl;
  for(int i = 0; i < cloud.points.size(); i++)
  {
    for(int j = 0; j < 4; j++)
    {
      std::cout << "cloud.points[" << i << "].data[" << j << "] = " << cloud.points[i].data[j]
        << std::endl;
    }
  }

  std::cout << "*** points.normal = " << std::endl;
  for(int i = 0; i < cloud.points.size(); i++)
  {
    for(int j = 0; j < 4; j++)
    {
      std::cout << "cloud.points[" << i << "].normal[" << j << "] = " << cloud.points[i].normal[j]
        << std::endl;
    }
  }
#endif
}

// main function
int
main (int argc, char** argv)
{
  // show help
  if(argc != 2)
  {
    showHelp(argv[0]);
    return 0;
  }

  // fetch point cloud filename in arguments
  std::string input_filename(argv[1]);

  // create datastructure to store points
  typedef  pcl::PointXYZ  PointType;
  pcl::PointCloud< PointType >::Ptr cloud(new pcl::PointCloud< PointType >());

  // load file
  pcl::ASCIIReader reader;
  reader.setInputFields< PointType >();
  reader.setExtension(".xyz");

  if(reader.read(input_filename, *cloud) < 0)
  {
    std::cout << "Error loading point cloud from " << input_filename << std::endl
              << std::endl;
    return -1;
  }
  else
  {
    std::cout << "Point cloud successfully loaded from " << input_filename << std::endl;
  }

  // print cloud after loading
  std::cout << std::endl << "cloud after loading :" << std::endl;
  print_cloud(*cloud);

  // compute normals
  //typedef  pcl::Normal  NormalType;
  typedef  pcl::PointNormal  NormalType;
  pcl::NormalEstimation< PointType, NormalType > normal_estimation;
  normal_estimation.setInputCloud(cloud);

  pcl::search::KdTree< PointType >::Ptr kdtree(
      new pcl::search::KdTree< PointType >);
  normal_estimation.setSearchMethod(kdtree);

  pcl::PointCloud< NormalType > normals;
  normal_estimation.setRadiusSearch(0.03);
  normal_estimation.compute(normals);

#if 1
  // Copy the xyz info from cloud_xyz and add it to cloud_normals as the xyz field in PointNormals estimation is zero
  for(size_t i = 0; i < normals.points.size(); ++i)
  {
    normals.points[i].x = cloud->points[i].x;
    normals.points[i].y = cloud->points[i].y;
    normals.points[i].z = cloud->points[i].z;
  }
#endif

  // print cloud with normal
  std::cout << std::endl << "cloud with normals :" << std::endl;
  print_cloud(normals);


  // write Point cloud to file
  std::string output_filename("test.pcd");
  if(pcl::io::savePCDFile(output_filename, normals) < 0)
  {
    std::cout << "Error saving point cloud to " << output_filename << std::endl;
    return -1;
  }
  else
  {
    std::cout << "Point cloud successfully saved to " << output_filename << std::endl;
  }

  return 0;
}

#endif