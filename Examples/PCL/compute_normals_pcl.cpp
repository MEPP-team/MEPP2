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
