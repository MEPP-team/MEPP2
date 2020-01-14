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
#pragma once

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"
#include "FEVV/Filters/Generic/color_mesh.h"

#include <cmath>  // for std::pow, std::sqrt, std::abs
#include <Eigen/Dense>


namespace FEVV {
namespace Filters {

/**
 * \brief  Compute mean and stdev of curvature.
 * 
 * \param  pc       the point cloud
 * \param  v_curvm  the vertex curvature map
 * \param  mean     output value of mean
 * \param  stdev    output value of stdev
 */
template< typename PointCloud,
          typename VertexCurvatureMap >
void
compute_mean_stdev_curvature(const PointCloud &pc,
                             const VertexCurvatureMap &v_curvm,
                             double &mean,
                             double &stdev)
{
  // init
  double sum_x = 0.0;
  double sum_x2 = 0.0;
  size_t vertices_nbr = 0;

  // compute sum_x and sum_x2 in a single pass
  auto iterator_pair = vertices(pc);
  auto vi = iterator_pair.first;
  auto vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    auto curv = get(v_curvm, *vi);
    sum_x  += curv;
    sum_x2 += curv*curv;
    vertices_nbr++;
  }

  // compute mean
  mean = sum_x / vertices_nbr;

  // compute stdev
  stdev = std::sqrt((sum_x2 / vertices_nbr) - mean*mean);
}


/**
 * \brief  Convert a curvature map to a color map.
 *
 * \param  pc       the point cloud
 * \param  v_curvm  the vertex curvature map
 * \param  v_cm     a vertex color map to store the curvature as a color
 */
template< typename PointCloud,
          typename VertexCurvatureMap,
          typename VertexColorMap >
void
curvature_to_color(const PointCloud &pc,
                   const VertexCurvatureMap &v_curvm,
                   VertexColorMap &v_cm)
{
  // compute absolute value of curvature
  auto v_curvm_abs = FEVV::make_vertex_property_map< PointCloud, double >(pc);
  auto iterator_pair = vertices(pc);
  auto vi = iterator_pair.first;
  auto vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    auto curv = get(v_curvm, *vi);
    put(v_curvm_abs, *vi, std::abs(curv));
  }

  // compute mean and stdev of absolute curvature
  double curv_mean, curv_stdev;
  compute_mean_stdev_curvature(pc, v_curvm_abs, curv_mean, curv_stdev);
  std::cout << "curv_mean = " << curv_mean << std::endl; //TODO-elo-rm
  std::cout << "curv_stdev = " << curv_stdev << std::endl; //TODO-elo-rm

  // compute a restricted curvature range to eliminate outliers
  double curv_range_min = 0;
  double curv_range_max = curv_mean + 2*curv_stdev;
  std::cout << "curv_range_min= " << curv_range_min << std::endl; //TODO-elo-rm
  std::cout << "curv_range_max= " << curv_range_max << std::endl; //TODO-elo-rm

  // populate color map
  FEVV::Filters::color_vertices_from_map(
      pc, v_curvm_abs, v_cm, curv_range_min, curv_range_max);
  
  // DBG //TODO-elo-rm
  vi = iterator_pair.first;
  for(; vi != vi_end; ++vi)
  {
    auto curv = get(v_curvm, *vi);
    auto curv_abs = get(v_curvm_abs, *vi);
    auto color = get(v_cm, *vi);
    std::cout << "curv= " << curv << "  curv_abs=" << curv_abs << "  color=" << color << std::endl; //TODO-elo-rm
  }
}


/**
 * 
 * \brief  Compute the curvature at origin point using a set of neighbors
 *         given by a list of indices.
 *
 * \param origin : Point to be projected.
 * \param refpoints : Contains all points from ref points cloud.
 * \param indices : Index of points in refpoints cloud used to compute the projection.
 * \param proj : Reference containing the point resulting from the projection.
 * \param H : Reference containing the mean curvature of the projected point.
 *
 * \return Returns both the projection and the mean curvature (Referenced variables).
 */
template< typename Point,
          typename PointMap,
          typename NNIdsVector,
          typename GeometryTraits >
double
curvature_at_point(const Point &origin,
                   const PointMap &pm,
                   const NNIdsVector &neighbors_ids,
                   const GeometryTraits &gt)
{
	Eigen::Matrix3d M;
	M.setZero();
	Eigen::Vector3d mu;
	mu.setZero();
	size_t nneighbors = neighbors_ids.size();
  std::cout << "nneighbors = " << nneighbors << std::endl; //TODO-elo-rm

	for (size_t i = 0; i < nneighbors; ++i)
	{
		Point p = get(pm, neighbors_ids[i]);
		Eigen::Vector3d neighbor(gt.get_x(p), gt.get_y(p), gt.get_z(p));
		mu = mu + neighbor;
		M = M + neighbor*neighbor.transpose();
	}

	mu = mu / ((double)nneighbors);
	M = 1. / ((double)nneighbors)*M - mu*mu.transpose();

	//get local frame
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(M);

	Eigen::Vector3d t1 = eig.eigenvectors().col(2);
	Eigen::Vector3d t2 = eig.eigenvectors().col(1);
	Eigen::Vector3d  n = eig.eigenvectors().col(0);

	Eigen::MatrixXd A(nneighbors, 6);
	Eigen::VectorXd B(nneighbors);

	//build linear system
	for (size_t i = 0; i < nneighbors; ++i)
	{
		Point p = get(pm, neighbors_ids[i]);
		double xglob = gt.get_x(p) - gt.get_x(origin);
		double yglob = gt.get_y(p) - gt.get_y(origin);
		double zglob = gt.get_z(p) - gt.get_z(origin);
		Eigen::Vector3d v(xglob, yglob, zglob);
		double x = v.transpose()*t1;
		double y = v.transpose()*t2;
		double z = v.transpose()*n;

		A(i, 0) = x*x;
		A(i, 1) = y*y;
		A(i, 2) = x*y;
		A(i, 3) = x;
		A(i, 4) = y;
		A(i, 5) = 1;

		B(i) = z;
	}

	Eigen::VectorXd coeffs = A.colPivHouseholderQr().solve(B);

	//TODO-elo-rm  //corresponding point:
	//TODO-elo-rm  Eigen::Vector3d delta = coeffs(5)*n;
	//TODO-elo-rm  proj = origin + delta;

	//corresponding curvature
	double fxx = 2 * coeffs(0);
	double fyy = 2 * coeffs(1);
	double fxy = coeffs(2);
	double fx = coeffs(3);
	double fy = coeffs(4);

	double H = 0.5*((1 + fx*fx)*fyy + (1 + fy*fy)*fxx - 2 * fxy*fx*fy) / pow(1 + fx*fx + fy*fy, 1.5);

  return H;
}


/**
 * \brief Compute the curvature for each point of the point cloud using
 *        the k nearest neighbors.
 * 
 * \param  pc  the point cloud
 * \param  pm  the point map of the point cloud
 * \param  k   number of nearest neighbors
 * \param  v_curvm  a vertex property map to store the curvature value
 *                  (double) of each vertex
 * \param  v_cm     a vertex color map to store the curvature as a color
 * \param  gt  the geometry traits
 */
template< typename PointCloud,
          typename PointMap,
          typename VertexCurvatureMap,
          typename VertexColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< PointCloud > >
void
point_cloud_curvature(const PointCloud &pc,
                      const PointMap &pm,
                      unsigned int k, 
                      VertexCurvatureMap &v_curvm,
                      VertexColorMap &v_cm,
                      const GeometryTraits &gt)
{
  // create k-d tree
  auto kd_tree_ptr = create_kd_tree(pc);

  // compute curvature at each vertex
  //TODO-elo-rm  double cmin = std::numeric_limits<float>::max();
  //TODO-elo-rm  double cmax = -cmin;
  auto iterator_pair = vertices(pc);
  auto vi = iterator_pair.first;
  auto vi_end = iterator_pair.second;
  for(; vi != vi_end; ++vi)
  {
    // retrieve current point
    auto point = get(pm, *vi);

    // do kNN-search around current point
    auto result = kNN_search(*kd_tree_ptr, k, point, pc);
    auto neighbors_ids = result.first;

    // compute curvature
    double curv = curvature_at_point(point, pm, neighbors_ids, gt);
    std::cout << "curv = " << curv << std::endl; //TODO-elo-rm

    // store curvature into property map
    put(v_curvm, *vi, curv);

    //TODO-elo-rm  // keep min and max values of curvature
    //TODO-elo-rm  cmin = std::min(cmin, curv);
    //TODO-elo-rm  cmax = std::max(cmax, curv);
  }

  // convert curvature to color for display
  curvature_to_color(pc, v_curvm, v_cm);
}


/**
 * \brief Compute the curvature for each point of the point cloud.
 *        Helper function to simplify the call to the filter.
 * 
 * \param  pc  the point cloud
 * \param  pm  the point map of the point cloud
 * \param  k   number of nearest neighbors
 * \param  v_curvm  a vertex property map to store the curvature value
 *                  (double) of each vertex
 * \param  v_cm     a vertex color map to store the curvature 
 *                  as a color
 */
template< typename PointCloud,
          typename PointMap,
          typename VertexCurvatureMap,
          typename VertexColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< PointCloud > >
void
point_cloud_curvature(const PointCloud &pc,
                      const PointMap &pm,
                      unsigned int k, 
                      VertexCurvatureMap &v_curvm,
                      VertexColorMap &v_cm)
{
  GeometryTraits gt(pc);
  point_cloud_curvature(pc, pm, k, v_curvm, v_cm, gt);
}

} // namespace Filters
} // namespace FEVV
