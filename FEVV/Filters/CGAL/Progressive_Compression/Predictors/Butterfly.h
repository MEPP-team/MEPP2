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

#include "DeltaPredictor.h"

#include <vector>
#include <list>

namespace FEVV {
namespace Filters {

/** Butterfly prediction scheme, as implemented in the paper by Rossignac
 */
template< typename HalfedgeGraph,
          typename PointMap >
class Butterfly : public Predictor< HalfedgeGraph,
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
  Butterfly(
      HalfedgeGraph &g,
      KeptPosition< HalfedgeGraph, PointMap >
          *kp,
      PointMap &pm)
      : SuperClass(g, kp, pm), _delta_for_borders(g, kp, pm)
  {
    _alpha = 1.15;
    k_a = 0;
    k_b = 0;
    SuperClass::_type = FEVV::Filters::PREDICTION_TYPE::BUTTERFLY;
  }

  bool isOnBorder(vertex_descriptor v) const 
  {
    boost::iterator_range<
        CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
        iterator_range = CGAL::halfedges_around_target(v, SuperClass::_g);
    for(auto h_v : iterator_range)
    {
      if(CGAL::is_border(source(h_v, SuperClass::_g), SuperClass::_g))
      {
        return true;
      }
    }
    return false;
  }


  std::vector< Vector >
  ComputeResiduals(CollapseInfo< HalfedgeGraph, PointMap > &mem) override
  {

    bool border_case = isOnBorder(mem.get_vkept());
    if(!border_case)
    {

      halfedge_descriptor v3_to_vkept =
        std::get< 0 >(halfedge(mem.get_v3(), mem.get_vkept(), SuperClass::_g));

      halfedge_descriptor v4_to_vkept =
        std::get< 0 >(halfedge(mem.get_v4(), mem.get_vkept(), SuperClass::_g));

      Point A = mem.get_pos_v1();
      Point B = mem.get_pos_v2();

      k_a = 0;
      k_b = 0;

      halfedge_descriptor begin = v3_to_vkept, end = v4_to_vkept;
      if(mem.get_reverse())
      {
        begin = v4_to_vkept;
        end = v3_to_vkept;
        Point save = std::move(A);
        A = std::move(B);
        B = std::move(save);
      }
      APred = A; // no std::move here!
      BPred = B; // no std::move here!
      // Get vertices from each stencil (for both A and B)
      std::vector< Point > a_vertices, b_vertices, crown_a, crown_b;
      fill_around_and_crown(a_vertices, crown_a, v3_to_vkept, v4_to_vkept, k_a);
      fill_around_and_crown(b_vertices, crown_b, v4_to_vkept, v3_to_vkept, k_b);

      // Compute SA, SB (barycenters of direct neighbours)  CA, CB (barycenter
      // of crown neighbours)
      Point SA = ComputeS(a_vertices,
                          get(SuperClass::_pm, mem.get_v3()),
                          get(SuperClass::_pm, mem.get_v4()),
                          k_a);
      Point SB = ComputeS(b_vertices,
                          get(SuperClass::_pm, mem.get_v3()),
                          get(SuperClass::_pm, mem.get_v4()),
                          k_b);


      Point CA = ComputeC(crown_a, b_vertices, k_a);
      Point CB = ComputeC(crown_b, a_vertices, k_b);

      CApred = CA;
      CBpred = CB;
      SApred = SA;
      SBpred = SB;

      Vector D = SuperClass::_gt.sub_p(B, A);
      // computing DA' and DB' (vector between kept position and the prediction)
      std::vector< Vector > result =
          ComputeDPrimeAandB(CA, SA, CB, SB, mem.get_pos_vkept());
      Vector DA = SuperClass::_gt.sub_p(A, mem.get_pos_vkept());
      Vector DB = SuperClass::_gt.sub_p(B, mem.get_pos_vkept());
      DApred = std::move(DA);
      DBpred = std::move(DB);
      Vector Dprime = 
        SuperClass::_gt.scalar_mult(SuperClass::_gt.add_v(result[0], result[1]), 0.5f);

      Point Aprime = SuperClass::_gt.add_pv(mem.get_pos_vkept(), 
                       SuperClass::_gt.scalar_mult(Dprime, 0.5f));
      Point Bprime = SuperClass::_gt.add_pv(mem.get_pos_vkept(), 
                       SuperClass::_gt.scalar_mult(Dprime, -0.5f));

      Aprime = Point(std::round(SuperClass::_gt.get_x(Aprime)),
                     std::round(SuperClass::_gt.get_y(Aprime)),
                     std::round(SuperClass::_gt.get_z(Aprime)));

      Bprime = Point(std::round(SuperClass::_gt.get_x(Bprime)),
                     std::round(SuperClass::_gt.get_y(Bprime)),
                     std::round(SuperClass::_gt.get_z(Bprime)));

      if(mem.get_reverse())
      {
        Point save = std::move(Aprime);
        Aprime = std::move(Bprime);
        Bprime = std::move(save);
      }


      Dprime = SuperClass::_gt.sub_p(Bprime, Aprime);
      Vector E = SuperClass::_gt.sub_v(D, Dprime);
      auto tmp = SuperClass::_gt.add_pv(SuperClass::_gt.ORIGIN, E);
      E = Vector(std::round(SuperClass::_gt.get_x(tmp)), 
                 std::round(SuperClass::_gt.get_y(tmp)), 
                 std::round(SuperClass::_gt.get_z(tmp)));
      Dpred = D;
      std::vector< Vector > pushed_residual = {E};
      mem.record_error_prediction(pushed_residual);
      return pushed_residual;
    }
    else
    {
      return _delta_for_borders.ComputeResiduals(mem);
    }
  }


  std::vector< Vector > ComputeDPrimeAandB(const Point &CA,
                                           const Point &SA,
                                           const Point &CB,
                                           const Point &SB,
                                           const Point &kept_position)
  {

    // Convert SA,SB,CA,CB to Vector (easier for computation)
    Vector SA_vec = Vector(SuperClass::_gt.get_x(SA),
                           SuperClass::_gt.get_y(SA),
                           SuperClass::_gt.get_z(SA));
    Vector SB_vec = Vector(SuperClass::_gt.get_x(SB),
                           SuperClass::_gt.get_y(SB),
                           SuperClass::_gt.get_z(SB));


    Vector CA_vec = Vector(SuperClass::_gt.get_x(CA),
                           SuperClass::_gt.get_y(CA),
                           SuperClass::_gt.get_z(CA));
    Vector CB_vec = Vector(SuperClass::_gt.get_x(CB),
                           SuperClass::_gt.get_y(CB),
                           SuperClass::_gt.get_z(CB));

    Vector V = Vector(SuperClass::_gt.get_x(kept_position),
                      SuperClass::_gt.get_y(kept_position),
                      SuperClass::_gt.get_z(kept_position));


    // following the formula from the paper

    Vector alphaSA = SuperClass::_gt.scalar_mult(SA_vec, _alpha);
    Vector alphaSB = SuperClass::_gt.scalar_mult(SB_vec, _alpha);

    Vector opp_alphaCA = SuperClass::_gt.scalar_mult(CA_vec, (1 - _alpha));
    Vector opp_alphaCB = SuperClass::_gt.scalar_mult(CB_vec, (1 - _alpha));

    Vector left_part_A =
        SuperClass::_gt.scalar_mult(SuperClass::_gt.add_v(opp_alphaCA, alphaSA), k_a + 3);
    Vector left_part_B =
        SuperClass::_gt.scalar_mult(SuperClass::_gt.add_v(opp_alphaCB, alphaSB), k_b + 3);

    Vector right_part_A = SuperClass::_gt.scalar_mult(V, _alpha - k_a - 3);
    Vector right_part_B = SuperClass::_gt.scalar_mult(V, _alpha - k_b - 3);

    Vector DA_prime = SuperClass::_gt.add_v(left_part_A, right_part_A);
    Vector DB_prime = SuperClass::_gt.add_v(left_part_B, right_part_B);

    DA_prime = SuperClass::_gt.scalar_mult(DA_prime, 1 / (k_a + 3 + _alpha));
    DB_prime = SuperClass::_gt.scalar_mult(DB_prime, 1 / (k_b + 3 + _alpha));
    auto tmp  = SuperClass::_gt.add_pv(SuperClass::_gt.ORIGIN, DA_prime);
    auto tmp2 = SuperClass::_gt.add_pv(SuperClass::_gt.ORIGIN, DB_prime);

    DA_prime = Vector(std::round(SuperClass::_gt.get_x(tmp) * (-2)),
                      std::round(SuperClass::_gt.get_y(tmp) * (-2)),
                      std::round(SuperClass::_gt.get_z(tmp) * (-2)));
    DB_prime = Vector(std::round(SuperClass::_gt.get_x(tmp2) * (2)),
                      std::round(SuperClass::_gt.get_y(tmp2) * (2)),
                      std::round(SuperClass::_gt.get_z(tmp2) * (2)));

    std::vector< Vector > Ds = {DA_prime, DB_prime};
    return Ds;
  }

  std::pair< Point, Point > PlacePoints(const std::vector< Vector > &residuals,
                                        vertex_descriptor vkept,
                                        halfedge_descriptor h1,
                                        halfedge_descriptor h2,
                                        bool is_reverse) override
  {
    bool is_border_case = isOnBorder(vkept);
    if(!is_border_case)
    {

      const Point& current_point = get(SuperClass::_pm, vkept);

      // Case using 1 residual

      Vector E = residuals[0];

      // fill "a" vertices (vertices around a), count k_a (number of vertices
      // around a) and compute c_a ("crown" around a) vertices
      std::vector< Point > a_vertices, c_a_vertices;
      k_a = 0;
      fill_around_and_crown(a_vertices, c_a_vertices, h1, h2, k_a);
      // fill "b" vertices (vertices around b), count k_b (number of vertices
      // around b) and compute c_b
      // ("crown" around b ) vertices
      std::vector< Point > b_vertices, c_b_vertices;
      k_b = 0;
      fill_around_and_crown(b_vertices, c_b_vertices, h2, h1, k_b);
      // Compute CA, CB (prediction according to crown), SA and SB (prediction
      // according to direct neighbours)
      Point SA = ComputeS(a_vertices,
                          get(SuperClass::_pm, source(h1, SuperClass::_g)),
                          get(SuperClass::_pm, source(h2, SuperClass::_g)),
                          k_a);

      Point SB = ComputeS(b_vertices,
                          get(SuperClass::_pm, source(h1, SuperClass::_g)),
                          get(SuperClass::_pm, source(h2, SuperClass::_g)),
                          k_b);

      Point CA = ComputeC(c_a_vertices, b_vertices, k_a);
      Point CB = ComputeC(c_b_vertices, a_vertices, k_b);

      std::vector< Vector > Ds =
          ComputeDPrimeAandB(CA, SA, CB, SB, current_point);

      //Vector DprimeA = Ds[0];
      //Vector DprimeB = Ds[1];
      // we get correct DA and DB
      Vector Dprime = 
           SuperClass::_gt.scalar_mult(SuperClass::_gt.add_v(Ds[0], Ds[1]), 0.5f);

      Point Aprime = SuperClass::_gt.add_pv(current_point, SuperClass::_gt.scalar_mult(Dprime, 0.5f));
      Point Bprime = SuperClass::_gt.add_pv(current_point, SuperClass::_gt.scalar_mult(Dprime, -0.5f));

      Aprime = Point(std::round(SuperClass::_gt.get_x(Aprime)),
                     std::round(SuperClass::_gt.get_y(Aprime)),
                     std::round(SuperClass::_gt.get_z(Aprime)));
      Bprime = Point(std::round(SuperClass::_gt.get_x(Bprime)),
                     std::round(SuperClass::_gt.get_y(Bprime)),
                     std::round(SuperClass::_gt.get_z(Bprime)));
      Dprime = SuperClass::_gt.sub_p(Bprime, Aprime);

      Vector D = SuperClass::_gt.add_v(Dprime, E);

      if(SuperClass::_kp->getType() == VKEPT_POSITION::MIDPOINT)
      {
        Vector DA = SuperClass::_gt.scalar_mult(D, -0.5f);
        auto tmp = SuperClass::_gt.add_pv(SuperClass::_gt.ORIGIN, DA);
        DA = Vector(std::floor(SuperClass::_gt.get_x(tmp)),
                    std::floor(SuperClass::_gt.get_y(tmp)),
                    std::floor(SuperClass::_gt.get_z(tmp)));
        Vector DB = SuperClass::_gt.add_v(D, DA);

        Point A = SuperClass::_gt.add_pv(current_point, DA);

        Point B = SuperClass::_gt.add_pv(current_point, DB);

        std::pair< Point, Point > resulting_points;

        resulting_points = std::make_pair(A, B);
        return resulting_points;
      }
      if(SuperClass::_kp->getType() == VKEPT_POSITION::HALFEDGE)
      {
        Vector residual;
        Point A;
        Point B;

        if(_rev)
        {
          D = SuperClass::_gt.scalar_mult(D, -1);
        }
        if(_rev)
        {

          A = SuperClass::_gt.add_pv(current_point, D); 
          B = current_point;
        }
        else
        {
          A = current_point;
          B = SuperClass::_gt.add_pv(current_point, D); 
        }


        std::pair< Point, Point > resulting_points;
        resulting_points = std::make_pair(A, B);
        return resulting_points;
      }
      return std::make_pair(Point(0, 0, 0), Point(0, 0, 0));
    }
    else
    {
      return _delta_for_borders.PlacePoints(residuals, 
                                            vkept, 
                                            h1, 
                                            h2, 
                                            is_reverse);    
	}
  }

  const std::tuple< bool, bool, bool, bool >& getInfoMidPoint() const
  {
    return _round_midpoint;
  }
  void set_bit_info(bool b1, bool b2, bool b3, bool b4)
  {
    _round_midpoint = std::make_tuple(b1, b2, b3, b4);
  }
  void set_rev(bool b) {
	  _rev = b; _delta_for_borders.set_rev(b);
  }

  std::string getMethodasString() const override { return "butterfly"; }


  // used for debug: to record the stencil around a vertex, in order to color
  // them later
  void record_stencil(std::vector< vertex_descriptor > &around,
                      std::vector< vertex_descriptor > &crown,
                      halfedge_descriptor begin,
                      halfedge_descriptor end,
                      int &k)
  {
    halfedge_descriptor begin_b =
        opposite(next(begin, SuperClass::_g), SuperClass::_g);
    halfedge_descriptor first_crown_b = next(
        opposite(prev(begin, SuperClass::_g), SuperClass::_g), SuperClass::_g);


    k = 0;
    crown.push_back(target(first_crown_b, SuperClass::_g));
    while(begin_b != end)
    {
      around.push_back(source(begin_b, SuperClass::_g));
      halfedge_descriptor current_crown =
          next(opposite(prev(begin_b, SuperClass::_g), SuperClass::_g),
               SuperClass::_g);
      crown.push_back(target(current_crown, SuperClass::_g));
      begin_b = opposite(next(begin_b, SuperClass::_g), SuperClass::_g);
      k++;
    }
  }


private:
  FEVV::Filters::DeltaPredictor< HalfedgeGraph,
                                 PointMap >
      _delta_for_borders;
  Point _estimated_a;
  Point _estimated_b;

  double _alpha;

  int k_a, k_b;
  bool _rev;

  Point CApred, CBpred, SApred, SBpred, APred, BPred;
  Vector DApred, DBpred;
  Vector Dpred;

  std::tuple< bool, bool, bool, bool > _round_midpoint;

  void fill_around_and_crown(std::vector< Point > &around,
                             std::vector< Point > &crown,
                             halfedge_descriptor begin,
                             halfedge_descriptor end,
                             int &k)
  {
    halfedge_descriptor begin_b =
        opposite(next(begin, SuperClass::_g), SuperClass::_g);
    halfedge_descriptor first_crown_b = next(
        opposite(prev(begin, SuperClass::_g), SuperClass::_g), SuperClass::_g);


    k = 0;
    crown.push_back(
        get(SuperClass::_pm, target(first_crown_b, SuperClass::_g)));
    while(begin_b != end)
    {
      around.push_back(get(SuperClass::_pm, source(begin_b, SuperClass::_g)));
      halfedge_descriptor current_crown =
          next(opposite(prev(begin_b, SuperClass::_g), SuperClass::_g),
               SuperClass::_g);
      crown.push_back(
          get(SuperClass::_pm, target(current_crown, SuperClass::_g)));
      begin_b = opposite(next(begin_b, SuperClass::_g), SuperClass::_g);
      k++;
    }
  }

  Point ComputeS(std::vector< Point > &around, Point v1, Point v2, int &k)
  {
    // Compute S
    Point S(0, 0, 0);
    for(auto it = around.begin(); it != around.end(); it++)
    {
      S = Point(SuperClass::_gt.get_x(S) + SuperClass::_gt.get_x(*it),
                SuperClass::_gt.get_y(S) + SuperClass::_gt.get_y(*it),
                SuperClass::_gt.get_z(S) + SuperClass::_gt.get_z(*it));
    }

    S = Point(SuperClass::_gt.get_x(S) + SuperClass::_gt.get_x(v1) +
                  SuperClass::_gt.get_x(v2),
              SuperClass::_gt.get_y(S) + SuperClass::_gt.get_y(v1) +
                  SuperClass::_gt.get_y(v2),
              SuperClass::_gt.get_z(S) + SuperClass::_gt.get_z(v1) +
                  SuperClass::_gt.get_z(v2));
    S = Point(SuperClass::_gt.get_x(S) / static_cast< double >((k + 3)),
              SuperClass::_gt.get_y(S) / static_cast< double >((k + 3)),
              SuperClass::_gt.get_z(S) / static_cast< double >((k + 3)));
    return S;
  }

  Point ComputeC(std::vector< Point > &crown,
                 std::vector< Point > &around_opposite,
                 int &k)
  {
    Point b_one = Point(0, 0, 0);
    Point b_kb = Point(0, 0, 0);
    if(around_opposite.size() > 0)
    {

      b_one = around_opposite[0];
      b_kb = around_opposite.back();
    }
    Point C = Point(SuperClass::_gt.get_x(b_one) + SuperClass::_gt.get_x(b_kb),
                    SuperClass::_gt.get_y(b_one) + SuperClass::_gt.get_y(b_kb),
                    SuperClass::_gt.get_z(b_one) + SuperClass::_gt.get_z(b_kb));

    for(auto it = crown.begin(); it != crown.end(); it++)
    {
      C = Point(SuperClass::_gt.get_x(C) + SuperClass::_gt.get_x(*it),
                SuperClass::_gt.get_y(C) + SuperClass::_gt.get_y(*it),
                SuperClass::_gt.get_z(C) + SuperClass::_gt.get_z(*it));
    }
    C = Point(SuperClass::_gt.get_x(C) / static_cast< double >((k + 3)),
              SuperClass::_gt.get_y(C) / static_cast< double >((k + 3)),
              SuperClass::_gt.get_z(C) / static_cast< double >((k + 3)));
    return C;
  }


  void fill_crown_around_vertex(std::list< vertex_descriptor > &crown,
                                std::list< vertex_descriptor > &around,
                                vertex_descriptor center,
                                int &k)
  {
    k = 0;
    boost::iterator_range<
        CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
        iterator_range = CGAL::halfedges_around_target(center, SuperClass::_g);
    for(auto h_v : iterator_range)
    {
      around.push_back(source(h_v, SuperClass::_g));

      halfedge_descriptor next_h =
          next(opposite(h_v, SuperClass::_g), SuperClass::_g);

      halfedge_descriptor opp_next_h = opposite(next_h, SuperClass::_g);


      if(!CGAL::is_border(opp_next_h, SuperClass::_g) &&
         !CGAL::is_border(
             next_h, SuperClass::_g)) // second verification might be overkill
      {

        halfedge_descriptor next_to_crown = next(opp_next_h, SuperClass::_g);
        crown.push_back(target(next_to_crown, SuperClass::_g));
      }
      k += 1;
    }
  }
};
} // namespace Filters
} // namespace FEVV
