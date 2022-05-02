// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include "FEVV/Wrappings/Geometry_traits.h"

#include <list>
#include <vector>
#include <utility>

#if defined(CGAL_LINKED_WITH_TBB)
#define TAG CGAL::Parallel_tag
#else
#define TAG CGAL::Sequential_tag
#endif


namespace FEVV {
	namespace Filters {


		template <typename HalfedgeGraph,
			typename PointMap,
			typename face_iterator = typename boost::graph_traits<HalfedgeGraph>::face_iterator,
			typename Vector = typename FEVV::Geometry_traits<HalfedgeGraph>::Vector,
			typename Point = typename FEVV::Geometry_traits<HalfedgeGraph>::Point>

    class GeometricMetrics
		{
		public:
      //typedef typename FEVV::Geometry_traits<HalfedgeGraph>::Kernel K;
			typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
			typedef CGAL::AABB_face_graph_triangle_primitive<HalfedgeGraph> AABB_Primitive;
			typedef CGAL::AABB_traits<CGALKernel, AABB_Primitive>           AABB_Traits;
			typedef CGAL::AABB_tree<AABB_Traits>                            AABB_Tree;

			GeometricMetrics(const HalfedgeGraph& LoD_init, const PointMap& pm_init) : _LoD_init(LoD_init), _pm_init(pm_init)
			{
              initialize_AABB_tree_for_init_LoD();

              subsample_LoD_init();
			}
			
			double compute_hausdorff(const HalfedgeGraph& LoD, int cas)
			{
        if (cas == 0) {
          return CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance<TAG>(LoD,
                                                                                    _LoD_init,
                 CGAL::Polygon_mesh_processing::parameters::number_of_points_per_area_unit(_number_of_points_per_area_unit_LoD_init));
        }
        else
        { // Computes the approximate Hausdorff distance from _LoD_init to LoD by returning the distance of the farthest point from LoD 
          // amongst a sampling of _LoD_init generated with the function sample_triangle_mesh() with _LoD_init and np1 as parameter.
          return CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance<TAG>(_LoD_init,
                                                                                    LoD,
                 CGAL::Polygon_mesh_processing::parameters::number_of_points_per_area_unit(_number_of_points_per_area_unit_LoD_init));
        }
			}
    protected:
			void initialize_AABB_tree_for_init_LoD()
			{
				std::pair<face_iterator, face_iterator> f_range_it= faces(_LoD_init);
				face_iterator orig_begin = f_range_it.first;
				face_iterator orig_end = f_range_it.second;

        _LoD_init_tree.clear();

				_LoD_init_tree.insert(orig_begin, orig_end, _LoD_init);
				_LoD_init_tree.accelerate_distance_queries();
			}

			void subsample_LoD_init()
			{
              _vec_distorsion.clear();
              std::vector<double> empty_vector;
              std::swap(_vec_distorsion, empty_vector);
              _vec_distorsion.push_back(0.0);

              _samples_LoD_init.clear();
              std::list<Point> empty_list;
              std::swap(_samples_LoD_init, empty_list);

              _area = CGAL::Polygon_mesh_processing::area(_LoD_init);
              auto nb_v = FEVV::size_of_vertices(_LoD_init);
              _grid_spacing = std::sqrt(_area / FEVV::size_of_faces(_LoD_init) * 4 / std::sqrt(3));// we assume that they are only triangles

              // We assume the input mesh has enough points on its surface for a distance
              // computation, therefore only vertex positions are used as samples
              CGAL::Polygon_mesh_processing::sample_triangle_mesh(_LoD_init, 
                                                            std::back_inserter(_samples_LoD_init), 
                                                            // the following set of parameters allow an RMSE computation
                                                            CGAL::Polygon_mesh_processing::parameters::use_grid_sampling(true)
                                                            .grid_spacing(_grid_spacing)
                                                             // the following set of parameters does not allow an RMSE computation:
                                                             // CGAL::Polygon_mesh_processing::parameters::use_random_uniform_sampling(true)
                                                             // .do_sample_faces(false)
                                                             // .do_sample_edges(false)
                                                             // .do_sample_vertices(true)
                                                            );
              // the value of this attribute depends on
              // the selected parameters for sample_triangle_mesh
              
              _number_of_points_per_area_unit_LoD_init = nb_v / std::max(0.001, _area) ;
              //std::cout << "area = " << CGAL::Polygon_mesh_processing::area(_LoD_init) << " ; _number_of_points_per_area_unit_LoD_init = " << _number_of_points_per_area_unit_LoD_init << std::endl;
              
			}
    public:
			double compute_L2(const HalfedgeGraph& mesh_to_sample, 
                        int cas,
                        bool compute_RMSE_instead_of_max // note that mean is different from root mean square
                        )
			{
				FEVV::Geometry_traits<HalfedgeGraph> gt(_LoD_init);

				using Point_and_primitive_id = typename AABB_Tree::Point_and_primitive_id;

					double L2_error = 0.0;
					if (cas == 0) // use AABB_tree of init mesh (the finest one) and sample the mesh_to_sample mesh 
					{ // distance D(mesh_to_sample, _LoD_init) = min/mean/max over all samples p of min ||p - p_LoD_init||
						// 1) subsampling of the current LoD
						std::list<Point> samples;
                        CGAL::Polygon_mesh_processing::sample_triangle_mesh(mesh_to_sample,
                          std::back_inserter(samples),

                          // the following set of parameters allow an RMSE computation
                          CGAL::Polygon_mesh_processing::parameters::use_grid_sampling(true)
                          .grid_spacing(_grid_spacing)
                          // the following set of parameters does not allow an RMSE computation:
                          //CGAL::Polygon_mesh_processing::parameters::use_random_uniform_sampling(true)
                          //.do_sample_faces(true)
                          //.do_sample_edges(true)
                          //.do_sample_vertices(true)
                          //.number_of_points_per_area_unit(_number_of_points_per_area_unit_LoD_init)
                                                  );
                        std::cout << "L2 cas 0: " << samples.size() << " samples" << std::endl;
						// 2) compute distance between points and init surface (_LoD_init)
						for (auto it = samples.begin(); it != samples.end(); ++it)
						{
							// computes closest point and primitive id
							Point_and_primitive_id pp = _LoD_init_tree.closest_point_and_primitive(*it);
							// get euclidean distance from point *it to surface mesh _LoD_init
							Vector d = gt.sub_p(pp.first, *it);
              if(compute_RMSE_instead_of_max)
                L2_error += gt.length2(d); 
              else
              { // max is better than either mean or min
                auto dist = gt.length(d);
                if (dist > L2_error)
                  L2_error = dist;
              }
						}
            if (compute_RMSE_instead_of_max)
              L2_error = std::sqrt( L2_error / (double)samples.size() ); // RMSE
					}
					else // opposite case
					{ // distance D(_LoD_init, mesh_to_sample) = min/mean/max over all samples p of min ||p - p_mesh_to_sample||
						std::pair<face_iterator, face_iterator> f_range_it = faces(mesh_to_sample);
						face_iterator orig_begin = f_range_it.first;
						face_iterator orig_end = f_range_it.second;

						AABB_Tree tree(orig_begin, orig_end, mesh_to_sample);
						tree.accelerate_distance_queries();
            std::cout << "L2 cas 1: " << _samples_LoD_init.size() << " init samples" << std::endl;
						// Searching for the closest point and facet for each vertex and compute L2
						for (auto it = _samples_LoD_init.begin(); it != _samples_LoD_init.end(); ++it)
						{
							// computes closest point and primitive id
							Point_and_primitive_id pp = tree.closest_point_and_primitive(*it);
                            // get euclidean distance from point *it to surface mesh mesh_to_sample
							Vector d = gt.sub_p(pp.first, *it);
              if (compute_RMSE_instead_of_max)
                L2_error += gt.length2(d); 
              else
              { // max is better than either mean or min
                auto dist = gt.length(d);
                if (dist > L2_error)
                  L2_error = dist;
              }
						}
            if (compute_RMSE_instead_of_max)
              L2_error = std::sqrt( L2_error / (double)_samples_LoD_init.size() ); // RMSE
					}
          //std::cout << "computer L2 error = " << L2_error << " ; CGAL Hausdorff value = " << compute_hausdorff(mesh_to_sample, cas) << std::endl;
          //std::cout << "computer L2 error = " << L2_error << std::endl;
					return L2_error;
			}

			double compute_symmetric_L2(const HalfedgeGraph& LoD, bool compute_RMSE_instead_of_max)
			{
				double d_LoD_to_mesh_init = this->compute_L2(LoD, 0, compute_RMSE_instead_of_max);
				double d_mesh_init_to_LoD = this->compute_L2(LoD, 1, compute_RMSE_instead_of_max);

                // Hausdorf or RMSE distances: the symmetrical distance takes the max
               double d_max = (d_LoD_to_mesh_init < d_mesh_init_to_LoD) ? d_mesh_init_to_LoD : d_LoD_to_mesh_init;

              // save distorsion values
              _vec_distorsion.push_back(d_max);

              return d_max;
			}

			const std::vector<double>& get_vec_distorsion() const { return _vec_distorsion; }

		private:
			HalfedgeGraph _LoD_init; // copy of init mesh: cannot be const (error with AABB otherwise)
			PointMap _pm_init; // copy of init pm
			AABB_Tree _LoD_init_tree;
			std::list<Point> _samples_LoD_init;

      double _number_of_points_per_area_unit_LoD_init;
      double _area, _grid_spacing;

			std::vector<double> _vec_distorsion; // distortion values for all batches

		};
	}
}
