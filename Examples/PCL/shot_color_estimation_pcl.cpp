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

// *** MEPP2 code calling PCL native filters. ***
// *** Use FEVV::PCLPointCloud datastructure. ***

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
#include <pcl/features/shot_omp.h>


// main
int main(int argc, char *argv[])
{
  // parse arguments
  if(argc != 2)
  {
    std::cout << "Usage:  " << argv[0] << "  input_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0] << "  bun0.pcd" << std::endl;
    return EXIT_FAILURE;
  }
  std::string input_file = argv[1];

  //----------------------------------

  // create point cloud
  typedef  FEVV::PCLPointCloud  PointCloudT;
  PointCloudT pc;

  // load point cloud
  std::cout << "load point cloud" << std::endl;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file, pc, pmaps_bag);

  //----------------------------------

  // create vertex-normal property map
  using VertexNormalMap = typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                                      PointCloudT >::pmap_type;
  VertexNormalMap v_nm;
  std::cout << "create vertex-normal map" << std::endl;
  v_nm = make_property_map(FEVV::vertex_normal, pc);

  // store property map in property maps bag
  put_property_map(FEVV::vertex_normal, pc, pmaps_bag, v_nm);

  // compute shot color estimation
  std::cout << "compute shot color estimation" << std::endl;

  //----------------------------------
  // PCL specific section - begin
  //----------------------------------

  typedef  FEVV::PCLEnrichedPoint  PointType;
  typedef  FEVV::PCLEnrichedPoint  NormalType;
    // use the cloud point to store the normal

  // init
  std::vector< int > indices;
  indices.resize(pc.points.size());
  for(size_t i = 0; i < indices.size(); ++i)
    indices[i] = static_cast< int >(i);
  pcl::search::KdTree< PointType >::Ptr tree;
  tree.reset(new pcl::search::KdTree< PointType >(false));
  tree->setInputCloud(pc.makeShared());

  // estimate normals first
  double mr = 0.002;
  pcl::NormalEstimation< PointType, NormalType > n;
  n.setInputCloud(pc.makeShared());
  boost::shared_ptr< std::vector< int > > indicesptr(new std::vector< int >(indices));
  n.setIndices(indicesptr);
  n.setSearchMethod(tree);
  n.setRadiusSearch(20 * mr);
  n.compute(pc);

  // init kdtree
  pcl::search::KdTree< PointType >::Ptr rgbaTree;
  rgbaTree.reset(new pcl::search::KdTree< PointType >(false));

  // SHOT1344
  pcl::SHOTColorEstimationOMP< PointType, NormalType, pcl::SHOT1344 > shot1344(
      true, true);
  shot1344.setInputNormals(pc.makeShared());
  shot1344.setRadiusSearch(20 * mr);
  pcl::PointCloud< pcl::SHOT1344 >::Ptr shots1344(new pcl::PointCloud< pcl::SHOT1344 >());
  shot1344.setInputCloud(pc.makeShared());
  shot1344.setIndices(indicesptr);
  shot1344.setSearchMethod(rgbaTree);
  shot1344.compute(*shots1344);

  //----------------------------------
  // PCL specific section - end
  //----------------------------------

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

// inspired by pcl-1.8.1-src/test/features/test_shot_estimation.cpp

#include <pcl/point_cloud.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/io/pcd_io.h>
#include <pcl/features/shot.h>
#include <pcl/features/shot_omp.h>
#include "pcl/features/shot_lrf.h"
#include <pcl/features/3dsc.h>
#include <pcl/features/usc.h>

#include <iostream>
#include <cmath> // for fabs()
// must be the last include to avoid side effect of disabling NDEBUG
#undef NDEBUG // enable assert in release mode
#include <cassert>


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
}


#define EXPECT_NEAR(v1, v2, epsilon)   assert(fabs((v1)-(v2)) <= (epsilon))


// main function
int
main(void)
{
  // fetch point cloud filename in arguments
  std::string input_filename("bun0.pcd");

  // create datastructure to store points
  pcl::PointCloud<pcl::PointXYZ> cloud;

  // load file
  if(pcl::io::loadPCDFile(input_filename, cloud) < 0)
  {
    std::cout << "Failed to load file '" << input_filename << "'." << std::endl;
    std::cout << "This file is provided in PCL sources." << std::endl;
    return -1; // failure
  }

  // print cloud after loading
  //std::cout << std::endl << "cloud after loading :" << std::endl;
  //print_cloud(cloud);

  // init
  std::vector< int > indices;
  typedef pcl::search::KdTree<pcl::PointXYZ>::Ptr KdTreePtr;
  KdTreePtr tree;
  indices.resize(cloud.points.size());
  for(size_t i = 0; i < indices.size(); ++i)
    indices[i] = static_cast< int >(i);
  tree.reset(new pcl::search::KdTree< pcl::PointXYZ >(false));
  tree->setInputCloud(cloud.makeShared());

  // estimate normals first
  double mr = 0.002;
  pcl::NormalEstimation< pcl::PointXYZ, pcl::Normal > n;
  pcl::PointCloud< pcl::Normal >::Ptr normals(new pcl::PointCloud< pcl::Normal >());
  n.setInputCloud(cloud.makeShared());
  boost::shared_ptr< std::vector< int > > indicesptr(new std::vector< int >(indices));
  n.setIndices(indicesptr);
  n.setSearchMethod(tree);
  n.setRadiusSearch(20 * mr);
  n.compute(*normals);
  std::cout << "normals->points.size() = " << normals->points.size() << std::endl;

  // Create fake point cloud with colors
  pcl::PointCloud< pcl::PointXYZRGBA > cloudWithColors;
  for(int i = 0; i < static_cast< int >(cloud.points.size()); ++i)
  {
    pcl::PointXYZRGBA p;
    p.x = cloud.points[i].x;
    p.y = cloud.points[i].y;
    p.z = cloud.points[i].z;

    //p.rgba = ((i % 255) << 16) + (((255 - i) % 255) << 8) + ((i * 37) % 255);
    p.r = (i % 255);
    p.g = ((255 - i) % 255);
    p.b = ((i * 37) % 255);
    cloudWithColors.push_back(p);
  }

  // init kdtree
  pcl::search::KdTree< pcl::PointXYZRGBA >::Ptr rgbaTree;
  rgbaTree.reset(new pcl::search::KdTree< pcl::PointXYZRGBA >(false));

  // SHOT1344
  pcl::SHOTColorEstimationOMP< pcl::PointXYZRGBA, pcl::Normal, pcl::SHOT1344 > shot1344(true,
                                                                    true);
  shot1344.setInputNormals(normals);
  shot1344.setRadiusSearch(20 * mr);
  pcl::PointCloud< pcl::SHOT1344 >::Ptr shots1344(new pcl::PointCloud< pcl::SHOT1344 >());
  shot1344.setInputCloud(cloudWithColors.makeShared());
  shot1344.setIndices(indicesptr);
  shot1344.setSearchMethod(rgbaTree);

  // estimate
  shot1344.compute(*shots1344);

  // check result
  EXPECT_NEAR(shots1344->points[103].descriptor[10], 0.0020453099, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[11], 0.0021887729, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[21], 0.057930067, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[42], 0.011778189, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[53], 0.0065085669, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[54], 0.012025614, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[55], 0.0044803056, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[64], 0.0644530654, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[65], 0.0465045683, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[86], 0.011518310, 1e-5);

  EXPECT_NEAR(shots1344->points[103].descriptor[357], 0.0020453099, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[360], 0.0027993850, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[386], 0.0451327376, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[387], 0.0544394031, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[389], 0.0047547864, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[453], 0.0051176427, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[481], 0.0053625242, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[482], 0.012025614, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[511], 0.0057367259, 1e-5);
  EXPECT_NEAR(shots1344->points[103].descriptor[512], 0.048375979, 1e-5);

  std::cout << "All checks passed." << std::endl;

  return 0;
}

#endif