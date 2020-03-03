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

// -------------------------------------------------------------------------------------
//        Iterative weighted PCA for robust and edge-aware normal vector estimation
//--------------------------------------------------------------------------------------
// Julia Sanchez, Florence Denis, David Coeurjolly, Florent dupont, Laurent Trassoudaine, Paul Checchin
// Liris (Lyon), Institut Pascal (Clermont Ferrand)
// Région Auvergne Rhône Alpes ARC6
// Private use for reviewers only
// --------------------------------------------------------------------------------------


#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <utility>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <flann/flann.hpp>

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

using namespace std;

class RawCloud{
public:
    RawCloud(){}
    int Read(std::string filepath);
    void buildTree();
    Eigen::Matrix<float, Eigen::Dynamic, 3> * getPC(){return pointcloud_;}
    flann::Index<flann::L2<float>>* getTree(){return tree_;}
    float getResolution ();
    void rescale();

    Eigen::Matrix<float, Eigen::Dynamic, 3> *pointcloud_;

private:
    int ReadCSV(std::string filepath);
    FEVV::PMapsContainer pmaps_bag;
    flann::Index<flann::L2<float>>* tree_;
    float getNearestNeighborDistance(Eigen::Vector3f pt);
    void SearchFLANNTree(flann::Index<flann::L2<float>>* index,
                                Eigen::Vector3f& input,
                                std::vector<int>& indices,
                                std::vector<float>& dists,
                                int nn);

};
