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

#include "Predictor.h"

#include <vector>

namespace FEVV {
namespace Filters {

template<typename HalfedgeGraph,
		typename PointMap>
class Raw_positions : public Predictor<HalfedgeGraph,
                                      PointMap>
{
public:
    using vertex_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor;
    using halfedge_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor;
    using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;
    using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
    using Geometry = typename FEVV::Geometry_traits< HalfedgeGraph >;
	
	typedef Predictor<HalfedgeGraph,PointMap> 
	   Super_class;
	
  Raw_positions(HalfedgeGraph &g,
               Kept_position<HalfedgeGraph,PointMap> * kp,
               PointMap &pm): Super_class(g,kp,pm)
			   {
				   Super_class::_type = PREDICTION_TYPE::POSITION; 
				   Super_class::_nbResiduals= 2;
               }

  ~Raw_positions() {}


  std::vector< Vector > compute_residuals(vertex_descriptor v1,
                                         vertex_descriptor v2, 
                                         Point /*kept_position*/)
  {
    const Point& p1 =  get(Super_class::_pm, v1);
    const Point& p2 =  get(Super_class::_pm, v2);

    Vector vec1 = Super_class::_gt.sub_p(p1,Super_class::_gt.ORIGIN);
    Vector vec2 = Super_class::_gt.sub_p(p2, Super_class::_gt.ORIGIN);
    std::vector<Vector> residuals;
    residuals.reserve(2);
    residuals.push_back(std::move(vec1));
    residuals.push_back(std::move(vec2));

	  return residuals;
  }

  std::vector< Vector > compute_residuals(Collapse_info<HalfedgeGraph, PointMap> &mem) override
  {
    Vector vec1 = Super_class::_gt.sub_p(mem.get_pos_v1(), Super_class::_gt.ORIGIN);
    Vector vec2 = Super_class::_gt.sub_p(mem.get_pos_v2(), Super_class::_gt.ORIGIN);
    std::vector< Vector > residuals;
    residuals.reserve(2);
    residuals.push_back(std::move(vec1));
    residuals.push_back(std::move(vec2));
    if(mem.get_reverse())
    {
      Vector save = std::move(residuals[0]);
      residuals[0] = std::move(residuals[1]);
      residuals[1] = std::move(save);
	  }
    mem.record_error_prediction(residuals);
    return residuals;
  }

  std::pair< Point, Point > place_points(const std::vector<Vector> &residuals, 
                                        vertex_descriptor vkept, 
                                        halfedge_descriptor h1, 
                                        halfedge_descriptor h2, 
                                        bool is_reverse) override
  {
    Point p1 = Super_class::_gt.add_pv(Super_class::_gt.ORIGIN, residuals[0]); 
    Point p2 = Super_class::_gt.add_pv(Super_class::_gt.ORIGIN, residuals[1]);
    return std::make_pair(p1,p2);
  }

  std::vector<Vector> reverse() { 
    Vector save = std::move(Super_class::_residuals[0]); 
    Super_class::_residuals[0] = std::move(Super_class::_residuals[1]);
    Super_class::_residuals[1] = std::move(save);
    return Super_class::_residuals;
  }

  void set_rev(bool /*b*/) override {}

  std::string get_as_string() const override { return "position"; }
};
} // namespace Filters
} // namespace FEVV
