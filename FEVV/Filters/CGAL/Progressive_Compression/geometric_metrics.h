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

			GeometricMetrics(const HalfedgeGraph& LoD_init,
				const PointMap& pm_init) : _LoD_init(LoD_init), _pm_init(pm_init)
			{
				_vec_distorsion.push_back(0.0);
			}
			
			// Computes the approximate Hausdorff distance from tm1 to tm2 by returning the distance of the farthest point from tm2 
			// amongst a sampling of tm1 generated with the function sample_triangle_mesh() with tm1 and np1 as parameter.
			double compute_hausdorff(const HalfedgeGraph& LoD) 
			{
				return CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance<TAG>(_LoD_init, 
				                                                                          LoD, 
																						  CGAL::Polygon_mesh_processing::parameters::number_of_points_per_area_unit(4000));
			}

			void initialize_AABB_tree_for_init_LoD()
			{
				std::pair<face_iterator, face_iterator> f_range_it= faces(_LoD_init);
				face_iterator orig_begin = f_range_it.first;
				face_iterator orig_end = f_range_it.second;

				_LoD_init_tree.insert(orig_begin, orig_end, _LoD_init);
				_LoD_init_tree.accelerate_distance_queries();
			}

			void subsample_LoD_init()
			{
				CGAL::Polygon_mesh_processing::sample_triangle_mesh(_LoD_init, 
                                                            std::back_inserter(_samples_LoD_init), 
                                                            CGAL::Polygon_mesh_processing::parameters::number_of_points_per_area_unit(4000)
                                                            //CGAL::Polygon_mesh_processing::parameters::number_of_points_on_faces(10)
                                                            );
			}

			double compute_L2(HalfedgeGraph& mesh_AABB_tree, 
                        HalfedgeGraph& mesh_to_sample, 
                        int cas,
                        bool compute_mean_instead_of_max // note that mean is different from root mean square
                        )
			{
				FEVV::Geometry_traits<HalfedgeGraph> gt(_LoD_init);

				using Point_and_primitive_id = typename AABB_Tree::Point_and_primitive_id;

					double L2_error = 0.0;
					if (cas == 0) // use AABB_tree of init mesh (the finest one) and sample the mesh_to_sample mesh 
					{
						// 1) subsampling of the current LoD
						std::list<Point> samples;
						CGAL::Polygon_mesh_processing::sample_triangle_mesh(mesh_to_sample, 
                                                                std::back_inserter(samples), 
                                                                CGAL::Polygon_mesh_processing::parameters::number_of_points_per_area_unit(4000)
                                                                //CGAL::Polygon_mesh_processing::parameters::number_of_points_on_faces(100000)
                                                               );

						// 2) compute distance between points and init surface (mesh_AABB_tree)
						for (auto it = samples.begin(); it != samples.end(); ++it)
						{
							// computes closest point and primitive id
							Point_and_primitive_id pp = _LoD_init_tree.closest_point_and_primitive(*it);
							// get euclidean distance from point *it to surface mesh mesh_AABB_tree
							Vector d = gt.sub_p(pp.first, *it);
              auto dist = gt.length(d);
              if(compute_mean_instead_of_max)
                L2_error += gt.length(d); // to get the true mean L2, a surface integral must be computed
              else
              { // max is better than either mean or min
                if (dist > L2_error)
                  L2_error = dist;
              }
						}
            if (compute_mean_instead_of_max)
						  L2_error /= (double)samples.size(); // mean L2 approximation if the sampling is regular
                                                  // to get the true mean L2, a surface integral must be computed
					}
					else // opposite case
					{
						std::pair<face_iterator, face_iterator> f_range_it = faces(mesh_to_sample);
						face_iterator orig_begin = f_range_it.first;
						face_iterator orig_end = f_range_it.second;

						AABB_Tree tree(orig_begin, orig_end, mesh_to_sample);
						tree.accelerate_distance_queries();

						// Searching for the closest point and facet for each vertex and compute L2
						for (auto it = _samples_LoD_init.begin(); it != _samples_LoD_init.end(); ++it)
						{
							// computes closest point and primitive id
							Point_and_primitive_id pp = tree.closest_point_and_primitive(*it);
              // get euclidean distance from point *it to surface mesh mesh_to_sample
							Vector d = gt.sub_p(pp.first, *it);
              auto dist = gt.length(d);
              if (compute_mean_instead_of_max)
                L2_error += gt.length(d); // to get the true mean L2, a surface integral must be computed
              else
              { // max is better than either mean or min
                if (dist > L2_error)
                  L2_error = dist;
              }
						}
            if (compute_mean_instead_of_max)
						  L2_error /= (double)_samples_LoD_init.size(); // mean L2 approximation if the sampling is regular
                                                            // to get the true L2, a surface integral must be computed 
					}
					return L2_error;
			}

			double compute_symmetric_L2(HalfedgeGraph& LoD, bool firstCall, bool compute_mean_instead_of_max = false)
			{
				if (firstCall)
				{
					this->initialize_AABB_tree_for_init_LoD();
					this->subsample_LoD_init();
				}

				double d_LoD_to_mesh_init = this->compute_L2(_LoD_init, LoD, 0, compute_mean_instead_of_max);
				double d_mesh_init_to_LoD = this->compute_L2(LoD, _LoD_init, 1, compute_mean_instead_of_max);

        if (compute_mean_instead_of_max)
        {
          // save distorsion values
          _vec_distorsion.push_back((d_LoD_to_mesh_init + d_mesh_init_to_LoD) / 2.0);

          return (d_LoD_to_mesh_init + d_mesh_init_to_LoD) / 2.0;
        }
        else
        {
          // Hausdorf distance
          double d_max = (d_LoD_to_mesh_init < d_mesh_init_to_LoD) ? d_mesh_init_to_LoD : d_LoD_to_mesh_init;

          // save distorsion values
          _vec_distorsion.push_back(d_max);

          return d_max;
        }
			}

			const std::vector<double>& get_vec_distorsion() const { return _vec_distorsion; }

		private:
			HalfedgeGraph _LoD_init; // copy of init mesh
			PointMap _pm_init; // copy of init pm
			AABB_Tree _LoD_init_tree;
			std::list<Point> _samples_LoD_init;

			std::vector<double> _vec_distorsion; // distortion values for all batches

		};
	}
}
