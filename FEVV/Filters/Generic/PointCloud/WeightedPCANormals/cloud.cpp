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


#include "cloud.h"
#include <pcl/search/impl/search.hpp>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree_flann.h>

// Reads file and puts it into pointcloud_
int RawCloud::ReadCSV(std::string filepath)
{
    int n = 0;
    vector<Eigen::Vector3f> res;
    std::cout<<"reading file : "<<filepath<<std::endl<<std::endl;
    std::ifstream fin(filepath);
    if (fin.is_open())
    {
        string test;
        std::getline ( fin, test, '\n' );

        int count = 0;
        for(int i =0; i<test.size(); ++i)
        {
            if(test[i]==',')
                ++count;
        }

        if(count==2)
        {
            while (std::getline ( fin, test, ',' ))
            {
                Eigen::Vector3f pt;
                pt(0) = stof(test);
                std::getline ( fin, test, ',' );
                pt(1) = stof(test);
                std::getline ( fin, test, '\n' );
                pt(2) = stof(test);

                res.push_back(pt);
                ++n;
            }
        }
        else if(count==5)
        {
            while (std::getline ( fin, test, ',' ))
            {
                Eigen::Vector3f pt;
                pt(0) = stof(test);
                std::getline ( fin, test, ',' );
                pt(1) = stof(test);
                std::getline ( fin, test, ',' );
                pt(2) = stof(test);
                std::getline ( fin, test, '\n' );

                res.push_back(pt);
                ++n;
            }
        }
        fin.close();
        pointcloud_ = new Eigen::Matrix<float, Eigen::Dynamic, 3>(res.size(), 3);
        for(int i = 0; i<res.size(); ++i)
            pointcloud_->row(i) = res[i];
    }
    else
    {
       std::cout<<"did not find file"<<std::endl<<std::endl;
    }
    std::cout<<"number of points  : "<<n<<std::endl<<std::endl;
    return n;
}

//builds tree to extract neighbors
void RawCloud::buildTree()
{
    int dim = pointcloud_->cols();
    int pc_size = pointcloud_->rows();

    std::vector<float> pc(pc_size * dim);
    flann::Matrix<float> flann_mat(&pc[0], pc_size, dim);

    int n = 0;

    for (int i =0; i < pc_size; ++i)
    {
        for (int j =0; j < dim; ++j)
        {
            pc[n] = (*pointcloud_)(i,j);
            ++n;
        }
    }

    tree_ = new flann::Index<flann::L2<float>>(flann_mat, flann::KDTreeSingleIndexParams(15));
    tree_->buildIndex();
}

//computes resolution as mean distance between nearest neighbors
float RawCloud::getResolution ()
{
  float res = 0.0;

  for (int i = 0; i<pointcloud_->rows(); ++i)
      res += getNearestNeighborDistance(pointcloud_->row(i));

  res /= pointcloud_->rows();
  return res;
}

void RawCloud::rescale ()
{
    float resolution = getResolution();
    std::cout<<"current resolution : "<<resolution<<std::endl<<std::endl;
    *pointcloud_ *= (float)(0.001)/resolution;
}

// Gets first neighbor distance
float RawCloud::getNearestNeighborDistance(Eigen::Vector3f pt)
{
    std::vector<float> dis;
    std::vector<int> neigh;

    SearchFLANNTree(tree_, pt, neigh, dis, 3);
    if(dis[1] != 0)
        return sqrt(dis[1]);
    else
        return sqrt(dis[2]);
}


//Search function in the tree to get the nearest neighbors
void RawCloud::SearchFLANNTree(flann::Index<flann::L2<float>>* index,
                            Eigen::Vector3f& input,
                            std::vector<int>& indices,
                            std::vector<float>& dists,
                            int nn)
{
    int dim = input.size();

    std::vector<float> query;
    query.resize(dim);
    for (int i = 0; i < dim; i++)
        query[i] = input(i);
    flann::Matrix<float> query_mat(&query[0], 1, dim);

    indices.resize(nn);
    dists.resize(nn);
    flann::Matrix<int> indices_mat(&indices[0], 1, nn);
    flann::Matrix<float> dists_mat(&dists[0], 1, nn);

    index->knnSearch(query_mat, indices_mat, dists_mat, nn, flann::SearchParams(128));
}
