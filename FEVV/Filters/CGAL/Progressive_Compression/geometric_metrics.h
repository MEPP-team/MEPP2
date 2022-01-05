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
			
			//computes the approximate Hausdorff distance from tm1 to tm2 by returning the distance of the farthest point from tm2 
			//amongst a sampling of tm1 generated with the function sample_triangle_mesh() with tm1 and np1 as parameter.
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
				CGAL::Polygon_mesh_processing::sample_triangle_mesh(_LoD_init, std::back_inserter(_samples_LoD_init), CGAL::Polygon_mesh_processing::parameters::number_of_points_on_faces(10));
			}

			double compute_L2(HalfedgeGraph& mesh_AABB_tree, HalfedgeGraph& mesh_to_sample, int cas)
			{
				FEVV::Geometry_traits<HalfedgeGraph> gt(_LoD_init);

					//using Object_and_primitive_id = typename AABB_Tree::Object_and_primitive_id;
					using Point_and_primitive_id = typename AABB_Tree::Point_and_primitive_id;


					double L2_error = 0.0;

					if (cas == 0) // correspond au cas où on fait l'AABB-tree sur le mesh de haut niveau et o� l'on sample le mesh de bas niveau
					{
						// subsample the current LoD
						std::list<Point> samples;
						CGAL::Polygon_mesh_processing::sample_triangle_mesh(mesh_to_sample, std::back_inserter(samples), CGAL::Polygon_mesh_processing::parameters::number_of_points_on_faces(100000));

						//Searching for the closest point and facet for each vertex and compute L2
						double L2_error = 0.0;
						for (auto it = samples.begin(); it != samples.end(); ++it)
						{
							// computes closest point and primitive id
							Point_and_primitive_id pp = _LoD_init_tree.closest_point_and_primitive(*it);
							//get distance from mesh1 to mesh2
							Vector d = gt.sub_p(pp.first, *it);
							L2_error += gt.length(d);

						}

						L2_error /= (double)samples.size();
					}
					else // correspond au cas opposé
					{

						std::pair<face_iterator, face_iterator> f_range_it = faces(mesh_AABB_tree);
						face_iterator orig_begin = f_range_it.first;
						face_iterator orig_end = f_range_it.second;

						AABB_Tree tree(orig_begin, orig_end, mesh_AABB_tree);
						tree.accelerate_distance_queries();

						//Searching for the closest point and facet for each vertex and compute L2
						for (auto it = _samples_LoD_init.begin(); it != _samples_LoD_init.end(); ++it)
						{
							// computes closest point and primitive id
							Point_and_primitive_id pp = tree.closest_point_and_primitive(*it);
							//get distance from mesh1 to mesh2
							Vector d = gt.sub_p(pp.first, *it);
							L2_error += gt.length(d);

						}

						L2_error /= (double)_samples_LoD_init.size();
					}

					return L2_error;
			
			}

			double compute_symetric_L2(HalfedgeGraph& LoD, bool firstCall)
			{
				if (firstCall)
				{
					this->initialize_AABB_tree_for_init_LoD();
					this->subsample_LoD_init();

				}

				double d_LoD_to_mesh_init = this->compute_L2(_LoD_init, LoD, 0);
				double d_mesh_init_to_LoD = this->compute_L2(LoD, _LoD_init, 1);

				// save distorsion values
				_vec_distorsion.push_back((d_LoD_to_mesh_init + d_mesh_init_to_LoD) / 2.0);

				return (d_LoD_to_mesh_init + d_mesh_init_to_LoD) / 2.0;
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
