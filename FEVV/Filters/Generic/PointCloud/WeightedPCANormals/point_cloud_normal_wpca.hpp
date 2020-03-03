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


#include <stdio.h>
#include <iostream>
#include <fstream>
//#include <omp.h> // Open mp disabled to avoid dependency 

#include "core.h"
#include "cloud.h"

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

namespace FEVV {
namespace Filters {

template< typename PointCloud,
          typename PointMap,
          typename NormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< PointCloud >
          >
void compute_weighted_pca(const PointCloud& pointCloud,
                          const PointMap&  	pointMap,
                          NormalMap&  normalMap,
                          int n_neight,
                          float noise,
                          float curvature,
                          const GeometryTraits &gt){

    /// Convert generic cloud into raw cloud
    /// Raw cloud is the structure used by the original normal computation algorithm
    RawCloud c;
    int cloud_size = num_vertices(pointCloud);
    auto iterator_pair = vertices(pointCloud);
    auto vi = iterator_pair.first;
    auto vi_end = iterator_pair.second;

    c.pointcloud_ = new Eigen::Matrix<float, Eigen::Dynamic, 3>(cloud_size+1, 3);

    vi = iterator_pair.first;
    int i = 0;
    for(; vi != vi_end; ++vi){
        auto p = get(pointMap, *vi);
        c.pointcloud_->row(i)[0] = gt.get_x(p);
        c.pointcloud_->row(i)[1] = gt.get_y(p);
        c.pointcloud_->row(i)[2] = gt.get_z(p);
        ++i;
    }

    /// Original algorithm
    c.buildTree();
    float res = c.getResolution();
    std::cout<<"cloud resolution : "<<res<<std::endl<<std::endl;

    Eigen::Matrix<float, Eigen::Dynamic, 3>* pc = c.getPC();
    flann::Index<flann::L2<float>>* tree = c.getTree();

    int n = 0;
    int inv = 0;
    int nan = 0;
    noise = std::max(noise, noise_min);

    std::cout<<"estimated noise="<<noise<<std::endl<<std::endl;
    std::cout<<"minimum curvature radius tolerated="<<curvature<<std::endl<<std::endl;
    std::cout<<"division factor="<<div_fact<<std::endl<<std::endl;

    std::vector<Eigen::Vector3f> normals(cloud_size);
    std::vector<Eigen::Vector3f> points(cloud_size);
    std::vector<int> onEdge(cloud_size);

    std::cout<<"------------------------------Computing normals------------------------------"<<std::endl<<std::endl;

    std::cout<<"Avancement : "<<std::endl;
    auto t_tot1 = std::chrono::high_resolution_clock::now();

// #pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads()) shared(onEdge, pc, tree, noise, n_neigh, curvature, normals, points) // comment if one thread used
    for (int i = 0; i < cloud_size; ++i)
    {
        CApp app(pc, tree, i, noise);
        Eigen::Vector3f point_ref = app.getPoint();
        app.setParams(div_fact, curvature);
        app.selectNeighborsKnn(n_neight);

        if(!(i%1000))
            std::cout<<((float)i/(float)cloud_size) * 100<<"%"<<std::endl;
        app.init1();

        if( app.isOnEdge() ) // when the mean projection error is greater than the noise-------------------------------------------------------------------------------------------------------
        {
            onEdge[i] = 1;

            //Compute first solution n_1------------------------------------------------------------------------------------------------------
            bool first = true;
            app.Optimize(first);
            app.OptimizePos(first, thresh_weight);

            //Compute second solution n_2------------------------------------------------------------------------------------------------------
            app.reinitPoint();
            app.init2();
            first = false;
            if(app.SuspectedOnEdge_)
            {
                app.Optimize(first);
                app.OptimizePos(first, thresh_weight);
            }

            //save result ------------------------------------------------------------------------------------------------------

            if( !app.isNan())
            {
                normals[i] = app.finalNormal_;
                points[i] = app.finalPos_;
            }
            else
            {
                normals[i]={0,0,0};
                points[i] = {0,0,0};
                std::cout<<"nan normal or point :"<<i<<std::endl;
                ++nan;
            }
        }
        else
        {
            onEdge[i] = 0;

            if( !app.isNan())
            {
                normals[i] = app.finalNormal_;
                points[i] = app.finalPos_;
            }
            else
            {
                std::cout<<"nan normal or point :"<<i<<std::endl;
                ++nan;
            }
        }

        points[i] = point_ref;

        if(normals[i].dot(points[i])>0)
            normals[i]=-normals[i];
    }

    auto t_tot2 = std::chrono::high_resolution_clock::now();
    typedef typename boost::property_traits< NormalMap >::value_type Normal;

    vi = iterator_pair.first;
    for(int i = 0; i < normals.size(); ++i){
        Normal normal(normals[i][0], normals[i][1], normals[i][2]);
        put(normalMap, *vi, normal);
        ++vi;
    }

    std::cout<<"total time to get normals :" <<std::chrono::duration_cast<std::chrono::milliseconds>(t_tot2-t_tot1).count()<<" milliseconds"<<std::endl<<std::endl;
}

template< typename PointCloud,
          typename PointMap,
          typename NormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< PointCloud >
          >
void compute_weighted_pca(const PointCloud& pointCloud,
                          const PointMap&  	pointMap,
                          NormalMap&  normalMap,
                          int n_neight,
                          float noise,
                          float curvature)
{
    GeometryTraits gt(pointCloud);
    compute_weighted_pca(pointCloud, pointMap, normalMap, n_neight, noise, curvature, gt);
}

}
}
