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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <list>
#include <vector>
#include <utility>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

namespace FEVV {
namespace Filters {

/**
 *  \brief Uniform_dequantization is a class dedicated to the XYZ uniform
 *         dequantization of vertex coordinates stored in the mesh point map.
 */
template<typename HalfedgeGraph,
         typename PointMap,
         typename Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector,
         typename Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point,
         typename vertex_descriptor =
         typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
         typename vertex_iterator =
         typename boost::graph_traits< HalfedgeGraph >::vertex_iterator >

class Uniform_dequantization
{

		public:
			Uniform_dequantization(const HalfedgeGraph &g,
				PointMap &pm,
				int nb_bits_quantization,
				const std::vector< double > &bb_dimension,
				const std::vector< double > &init_point)
				: _g(g), _pm(pm), _nb_bits_quantization(nb_bits_quantization),
				_bb_dimension(bb_dimension), _init_point(init_point)
			{
				this->set_max_length();
				this->set_quantization_step();
			}
		protected:
			void set_max_length()
			{
				_max_length = _bb_dimension[0];
				if(_max_length < _bb_dimension[1])
					_max_length = _bb_dimension[1];
				if(_max_length < _bb_dimension[2])
					_max_length = _bb_dimension[2];
			}

			void set_quantization_step()
			{
				_quantization_step = _max_length / (pow(2.0, _nb_bits_quantization) - 1);
				std::cout << "quantization step : " << _quantization_step << std::endl;
			}
		public:
            /// Dequantizes a vertex position according to dequantization 
            /// parameters.
            /// The mesh point map is not modified.
            /// Returns the dequantized XYZ vertex position. 		
			Point dequantize(vertex_descriptor v)
			{
				const Point& pq = get(_pm, v);

				return dequantize(pq);
			}
            /// Dequantizes a Point p according to dequantization parameters.
            /// The mesh point map is not modified.
            /// Returns the dequantized XYZ point. 
			Point dequantize(const Point& pq)
			{
				double p_min_x = _init_point[0];
				double p_min_y = _init_point[1];
				double p_min_z = _init_point[2];

				FEVV::Geometry_traits< HalfedgeGraph > gt(_g);
				uint32_t pq_x = static_cast<uint32_t>(gt.get_x(pq));
				uint32_t pq_y = static_cast<uint32_t>(gt.get_y(pq));
				uint32_t pq_z = static_cast<uint32_t>(gt.get_z(pq));

				// Reconstruction
				double p_x = (double)pq_x * _quantization_step + p_min_x;
				double p_y = (double)pq_y * _quantization_step + p_min_y;
				double p_z = (double)pq_z * _quantization_step + p_min_z;

				Point new_position(p_x, p_y, p_z);
				return new_position;
			}
			/// Dequantizes all vertex positions stored in the point map.
			void point_dequantization()
			{

				FEVV::Geometry_traits< HalfedgeGraph > gt(_g);
				double p_min_x = _init_point[0];
				double p_min_y = _init_point[1];
				double p_min_z = _init_point[2];

				vertex_iterator vi = vertices(_g).first;

				for( ; vi != vertices(_g).second; ++vi)
				{
					const Point& pq = get(_pm, *vi);

					uint32_t pq_x = static_cast<uint32_t>(gt.get_x(pq));
					uint32_t pq_y = static_cast<uint32_t>(gt.get_y(pq));
					uint32_t pq_z = static_cast<uint32_t>(gt.get_z(pq));

					// Left reconstruction: why do no use instead the centered reconstruction?
					double p_x = (double)pq_x * _quantization_step + p_min_x;
					double p_y = (double)pq_y * _quantization_step + p_min_y;
					double p_z = (double)pq_z * _quantization_step + p_min_z;

					Point new_position = Point(p_x, p_y, p_z);

					put(_pm, *vi, new_position);
				}
			}

            int get_nb_bits_quantization() const { return _nb_bits_quantization; }
		private:
			const HalfedgeGraph &_g; /// Topology remains the same
			PointMap &_pm;           /// Point map changes
			const int _nb_bits_quantization;
			const std::vector< double > _bb_dimension;
			const std::vector< double > _init_point;
			double _max_length;
			double _quantization_step;
};

} // namespace Filters
} // namespace FEVV
