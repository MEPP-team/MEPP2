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

#include "FEVV/Wrappings/Geometry_traits.h"

#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdint>
#include <cmath>
#include <set>
#include <list>
#include <vector>
#include <utility>
#include <cstdio>
#include <cstdlib>

namespace FEVV {
namespace Filters {
	
/**
 *  \brief Uniform_quantization is a class dedicated to the XYZ uniform
 *         quantization of vertex coordinates stored in the mesh point map.
 */
template <typename HalfedgeGraph,
          typename PointMap,
          typename Vector = typename FEVV::Geometry_traits<HalfedgeGraph>::Vector,
          typename Point = typename FEVV::Geometry_traits<HalfedgeGraph>::Point,
          typename Geometry = FEVV::Geometry_traits<HalfedgeGraph>,
          typename vertex_descriptor = typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor,
          typename vertex_iterator = typename boost::graph_traits<HalfedgeGraph>::vertex_iterator>
			
class Uniform_quantization 
{

			public:
				Uniform_quantization(const HalfedgeGraph& g, PointMap& pm, int nb_bits) : _g(g), _pm(pm), _nb_bits(nb_bits) 
				{
					this->find_min_and_max();
					this->set_bounding_box();
					this->set_quantization_step();
				}

			protected:
				void find_min_and_max()
				{
                    Geometry gt(_g);
                    vertex_iterator vi = vertices(_g).first;
                    if(vertices(_g).first == vertices(_g).second)
                      return;
                    _p_min = get(_pm, *vi);
                    _p_max = _p_min;
                    ++vi;
                    double x_min = gt.get_x(_p_min);
                    double y_min = gt.get_y(_p_min);
                    double z_min = gt.get_z(_p_min);
                    double x_max = gt.get_x(_p_max);
                    double y_max = gt.get_y(_p_max);
                    double z_max = gt.get_z(_p_max);
                    double x_current;
                    double y_current;
                    double z_current;
					for( ; vi != vertices(_g).second; ++vi)
					{
						Point p_current = get(_pm, *vi);

						// get coord
						x_current = gt.get_x(p_current);
						y_current = gt.get_y(p_current);
						z_current = gt.get_z(p_current);

						// get min
						if(x_current < x_min)
							x_min = x_current;
						if(y_current < y_min)
							y_min = y_current;
						if(z_current < z_min)
							z_min = z_current;

						// get max
						if(x_current > x_max)
							x_max = x_current;
						if(y_current > y_max)
							y_max = y_current;
						if(z_current > z_max)
							z_max = z_current;
					}
                    _p_min = Point(x_min, y_min, z_min);
                    _p_max = Point(x_max, y_max, z_max);
				}

				
				struct boundingBox
				{
					Point starting_point;
					Vector vl;
					Vector vh;
					Vector vp;

				};

				void set_bounding_box()
				{
                    Geometry gt(_g);

					_bb.starting_point = _p_min;
					_bb.vl = Vector(0, gt.get_y(_p_max)-gt.get_y(_p_min), 0);
					_bb.vh = Vector(0, 0, gt.get_z(_p_max)-gt.get_z(_p_min));
					_bb.vp = Vector(gt.get_x(_p_max)-gt.get_x(_p_min), 0, 0);

				}

				void set_quantization_step()
				{
                    Geometry gt(_g);
					double l = gt.length(_bb.vl);
					double h = gt.length(_bb.vh);
                    double p = gt.length(_bb.vp);

					_max_length.first = l;
					_max_length.second = 'y';
					if(_max_length.first < p)
					{
						_max_length.first = p;
						_max_length.second = 'x';
					}
					if(_max_length.first < h)
					{
						_max_length.first = h;
						_max_length.second = 'z';

					}
				
					_quantization_step = _max_length.first / pow(2.0, _nb_bits);
#if(DEBUG)					
					std::cout << "Uniform_quantization: max length : " << _max_length.first << std::endl;
					std::cout << "Uniform_quantization: quantization step : " << _quantization_step << std::endl;
#endif					
				}
			public:
				double get_quantization_step() const { return _quantization_step; }

                /// Quantizes all vertex positions stored in the point map.
				void point_quantization(bool recompute_quanti_param = false)
				{
                    Geometry gt(_g);
					///////////////////////////////////////////////////////////
					if(recompute_quanti_param) 
					{
						this->find_min_and_max();
						this->set_bounding_box();
						this->set_quantization_step();
					}
					///////////////////////////////////////////////////////////
					double p_min_x = gt.get_x(_p_min);
					double p_min_y = gt.get_y(_p_min);
					double p_min_z = gt.get_z(_p_min);

					vertex_iterator vi = vertices(_g).first;
					for( ; vi != vertices(_g).second; ++vi)
					{
						Point point = get(_pm, *vi);

						double px = gt.get_x(point);
						double py = gt.get_y(point);
						double pz = gt.get_z(point);

						uint32_t pq_x = static_cast<uint32_t>(round((px - p_min_x) / _quantization_step));
						if(pq_x >= pow(2.0, _nb_bits)) // may occur when px = p_max_x 
							pq_x = static_cast<uint32_t>(pow(2.0, _nb_bits) - 1);
						uint32_t pq_y = static_cast<uint32_t>(round((py - p_min_y) / _quantization_step));
						if(pq_y >= pow(2.0, _nb_bits)) // may occur when py = p_max_y 
							pq_y = static_cast<uint32_t>(pow(2.0, _nb_bits) - 1);
						uint32_t pq_z = static_cast<uint32_t>(round((pz - p_min_z) / _quantization_step));
						if(pq_z >= pow(2.0, _nb_bits)) // may occur when pz = p_max_z 
                            pq_z = static_cast<uint32_t>(pow(2.0, _nb_bits) - 1);

						Point new_position = Point(pq_x, pq_y, pq_z);
						put(_pm, *vi, new_position);
					}
				}
                /// Quantizes a Point p according to quantization parameters.
                /// The mesh point map is not modified.
                /// Returns the quantized XYZ point. 
				Point quantize(const Point& p)
				{
					Geometry gt(_g);

					double p_min_x = gt.get_x(_p_min);
					double p_min_y = gt.get_y(_p_min);
					double p_min_z = gt.get_z(_p_min);

					double px = gt.get_x(p);
					double py = gt.get_y(p);
					double pz = gt.get_z(p);

					uint32_t pq_x = static_cast<uint32_t>(round((px - p_min_x) / _quantization_step));
					if(pq_x >= pow(2.0, _nb_bits)) // may occur when px = p_max_x 
						pq_x = static_cast<uint32_t>(pow(2.0, _nb_bits) - 1);
					uint32_t pq_y = static_cast<uint32_t>(round((py - p_min_y) / _quantization_step));
					if(pq_y >= pow(2.0, _nb_bits)) // may occur when py = p_max_y 
						pq_y = static_cast<uint32_t>(pow(2.0, _nb_bits) - 1);
					uint32_t pq_z = static_cast<uint32_t>(round((pz - p_min_z) / _quantization_step));
					if(pq_z >= pow(2.0, _nb_bits)) // may occur when pz = p_max_z 
						pq_z = static_cast<uint32_t>(pow(2.0, _nb_bits) - 1);

					Point pq = Point(pq_x, pq_y, pq_z);
					return pq;
				}
																						   


				std::vector<double> get_bb_dimension() const
				{
                    Geometry gt(_g);
					double l = gt.length(_bb.vl);
					double p = gt.length(_bb.vp);
					double h = gt.length(_bb.vh);

					std::vector<double> bb_dimension{ p, l, h };
					return bb_dimension;
				}

				std::vector<double> get_init_coord() const
				{
                    Geometry gt(_g);
					double x = gt.get_x(_p_min);
					double y = gt.get_y(_p_min);
					double z = gt.get_z(_p_min);

					std::vector<double> init_coord{ x, y, z };
					return init_coord;
				}

				double get_diagonal() const
                {
                  Geometry gt(_g);
                  Vector diagonal = _bb.vl + _bb.vh + _bb.vp;
                  return gt.length(diagonal);
                }
				
                int get_nb_bits_quantization() const { return _nb_bits; }
			private:
				const HalfedgeGraph& _g; /// Topology remains the same
				PointMap& _pm;           /// Point map changes
				const int _nb_bits;
				double _quantization_step;
				std::pair<double,char> _max_length;
				Point _p_min;
				Point _p_max;
				boundingBox _bb;
};
} // namespace Filters
} // namespace FEVV