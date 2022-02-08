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
#include <map>

namespace FEVV {
namespace Filters {

template<
    typename HalfedgeGraph,
    typename PointMap >
class DeltaPredictor : public Predictor< HalfedgeGraph,
                                         PointMap >
{
public:
    using vertex_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor;
    using halfedge_descriptor =
            typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor;
    using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;
    using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
    using Geometry = typename FEVV::Geometry_traits< HalfedgeGraph >;
	
    typedef Predictor< HalfedgeGraph, PointMap > 
            SuperClass;
  DeltaPredictor(HalfedgeGraph &_g,
                 KeptPosition< HalfedgeGraph, PointMap > *kp,
                 PointMap &pm)
      : SuperClass(_g, kp, pm)
  {
    SuperClass::_type = PREDICTION_TYPE::DELTA;
    _round_midpoint =
        std::make_tuple< bool, bool, bool, bool >(false, false, false, false);
  }


  std::vector< Vector >  ComputeResiduals(CollapseInfo< HalfedgeGraph, PointMap > &mem) override
  { // vt <- v1, vs <- v2

    Point A = mem.get_pos_v1();
    Point B = mem.get_pos_v2();

    if(mem.get_reverse())
    {
      Point save = A;
      A = B;
      B = save;
    }
    Apred = A;
    Bpred = B;

    Vector delta = SuperClass::_gt.sub_p(B, A);

    Vector vec1 = // vt->vkept
       SuperClass::_gt.sub_p(mem.get_pos_v1(), mem.get_pos_vkept()); 
    Vector vec2 = // vs->vkept
       SuperClass::_gt.sub_p(mem.get_pos_v2(), mem.get_pos_vkept()); 
    std::vector< Vector > residuals;
    Vector to_be_pushed, revers;
    //auto copy_round = mem._midpoint_rounds;
    if(SuperClass::_gt.length(vec1) <= SuperClass::_gt.length(vec2))
    {
      to_be_pushed = std::move(vec1);
      revers = std::move(vec2);
    }
    else
    {
      to_be_pushed = std::move(vec2);
      revers = std::move(vec1);
    }
    if(SuperClass::_kp->getType() == VKEPT_POSITION::MIDPOINT)
    {
      residuals = {delta};
    }
    if(SuperClass::_kp->getType() == VKEPT_POSITION::HALFEDGE)
    {
      residuals = {delta};
    }
    mem.record_error_prediction(residuals);
    return residuals;
  }

  std::vector< Vector > Reverse(vertex_descriptor v)
  {
    SuperClass::_residuals[0] = _vkept_to_other_residual[v];
    return SuperClass::_residuals;
  }


  std::pair< Point, Point > PlacePoints(const std::vector< Vector > &residuals,
                                        vertex_descriptor vkept,
                                        halfedge_descriptor /*h1*/,
                                        halfedge_descriptor /*h2*/,
                                        bool /*is_reverse*/) override
  {
    const Point& current_point = get(SuperClass::_pm, vkept);

    Point p1, p2;
    if(residuals.size() == 2)
    {
      p1 = SuperClass::_gt.add_pv(current_point, residuals[0]);
      p2 = SuperClass::_gt.add_pv(current_point, residuals[1]);
    }
    else
    {
      if(SuperClass::_kp->getType() == VKEPT_POSITION::MIDPOINT)
      {
        Vector D = residuals[0];
        Vector DA = SuperClass::_gt.scalar_mult(D, -0.5f);
        auto&& tmp = SuperClass::_gt.add_pv(SuperClass::_gt.ORIGIN, DA);
        DA = Vector(std::floor(SuperClass::_gt.get_x(tmp)),
                    std::floor(SuperClass::_gt.get_y(tmp)),
                    std::floor(SuperClass::_gt.get_z(tmp)) );
        Vector DB = SuperClass::_gt.add_v(D, DA);

        p1 = SuperClass::_gt.add_pv(current_point, DA);

        p2 = SuperClass::_gt.add_pv(current_point, DB);
      }
      else
      {

        if(!_rev)
        {
          p1 = current_point;
          p2 = SuperClass::_gt.add_pv(current_point, residuals[0]);
        }
        else
        {
          p1 = SuperClass::_gt.sub_pv(current_point, residuals[0]);
          p2 = current_point;
        }
      }
    }


    return std::make_pair(p1, p2);
  }

   


  const std::tuple< bool, bool, bool, bool >& getInfoMidPoint() const
  {
    return _round_midpoint;
  }
  void set_bit_info(bool b1, bool b2, bool b3, bool b4) 
  {
    _round_midpoint = std::make_tuple(b1, b2, b3, b4);
  }
  void set_rev(bool b)  { _rev = b; }

  std::string getMethodasString() const override { return "delta"; }

  ~DeltaPredictor() {}

private:
  Point Apred, Bpred;
  Vector DApred, DBpred;
  Vector Dpred;
  std::tuple< bool, bool, bool, bool > _round_midpoint;
  Point _source;
  Point _target;
  std::map< vertex_descriptor, Vector > _vkept_to_other_residual;
  bool _rev;
};


} // namespace Filters

} // namespace FEVV
