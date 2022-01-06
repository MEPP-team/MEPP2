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
class RawPositions : public Predictor<HalfedgeGraph,
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
	   SuperClass;
	
  RawPositions(HalfedgeGraph &g,
               KeptPosition<HalfedgeGraph,PointMap> * kp,
               PointMap &pm): SuperClass(g,kp,pm){
				   SuperClass::_type = PREDICTION_TYPE::POSITION; 
				   SuperClass::_nbResiduals= 2;}

  ~RawPositions() {}


  std::vector< Vector > ComputeResiduals(vertex_descriptor v1,
                                         vertex_descriptor v2, 
                                         Point /*kept_position*/) override
  {

	const Point& p1 =  get(SuperClass::_pm, v1);
	const Point& p2 =  get(SuperClass::_pm, v2);

    Vector vec1 = SuperClass::_gt.sub_p(p1,SuperClass::_gt.ORIGIN);
	Vector vec2 = SuperClass::_gt.sub_p(p2, SuperClass::_gt.ORIGIN);
	std::vector<Vector> residuals;
	residuals.reserve(2);
	residuals.push_back(std::move(vec1));
    residuals.push_back(std::move(vec2));
	return residuals;
  }

  std::vector< Vector > ComputeResiduals(CollapseInfo<HalfedgeGraph, PointMap> &mem) override
  {
    Vector vec1 = SuperClass::_gt.sub_p(mem._pos_v1, SuperClass::_gt.ORIGIN);
    Vector vec2 = SuperClass::_gt.sub_p(mem._pos_v2, SuperClass::_gt.ORIGIN);
    std::vector< Vector > residuals;
    residuals.reserve(2);
    residuals.push_back(std::move(vec1));
    residuals.push_back(std::move(vec2));
    if(mem._reverse)
    {
      Vector save = std::move(residuals[0]);
      residuals[0] = std::move(residuals[1]);
      residuals[1] = std::move(save);
	}
    return residuals;
  }

  std::pair< Point, Point > PlacePoints(const std::vector<Vector> &residuals, 
                                        vertex_descriptor vkept, 
                                        halfedge_descriptor h1, 
                                        halfedge_descriptor h2, 
                                        bool is_reverse) override
  {
    Point p1 = SuperClass::_gt.add_pv(SuperClass::_gt.ORIGIN, residuals[0]); 
    Point p2 = SuperClass::_gt.add_pv(SuperClass::_gt.ORIGIN, residuals[1]);
    return std::make_pair(p1,p2);
  }

  std::vector<Vector> Reverse() override { 
    Vector save = std::move(SuperClass::_residuals[0]); 
    SuperClass::_residuals[0] = std::move(SuperClass::_residuals[1]);
    SuperClass::_residuals[1] = std::move(save);
    return SuperClass::_residuals;
  }

};
}
}
