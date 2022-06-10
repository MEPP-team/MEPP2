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

#include "Delta_predictor.h"

#include <vector>
#include <list>

namespace FEVV {
namespace Filters {

/** Butterfly prediction scheme, as implemented in the paper by Rossignac.
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
      Super_class;
  Butterfly(
      HalfedgeGraph &g,
      Kept_position< HalfedgeGraph, PointMap >
          *kp,
      PointMap &pm)
      : Super_class(g, kp, pm), _delta_for_borders(g, kp, pm)
  {
    _alpha = 1.15;
    k_a = 0;
    k_b = 0;
    Super_class::_type = FEVV::Filters::PREDICTION_TYPE::BUTTERFLY;
  }

  std::vector< Vector >
  compute_residuals(Collapse_info< HalfedgeGraph, PointMap > &mem) override
  {

    bool border_case = is_on_border(mem.get_vkept());
    if(!border_case)
    {

      halfedge_descriptor v3_to_vkept =
        std::get< 0 >(halfedge(mem.get_v3(), mem.get_vkept(), Super_class::_g));

      halfedge_descriptor v4_to_vkept =
        std::get< 0 >(halfedge(mem.get_v4(), mem.get_vkept(), Super_class::_g));

      Point A = mem.get_pos_vt();
      Point B = mem.get_pos_vs();

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
      Point SA = compute_S(a_vertices,
                          get(Super_class::_pm, mem.get_v3()),
                          get(Super_class::_pm, mem.get_v4()),
                          k_a);
      Point SB = compute_S(b_vertices,
                          get(Super_class::_pm, mem.get_v3()),
                          get(Super_class::_pm, mem.get_v4()),
                          k_b);


      Point CA = compute_C(crown_a, b_vertices, k_a);
      Point CB = compute_C(crown_b, a_vertices, k_b);

      CApred = CA;
      CBpred = CB;
      SApred = SA;
      SBpred = SB;

      Vector D = Super_class::_gt.sub_p(B, A);
      // computing DA' and DB' (vector between kept position and the prediction)
      std::vector< Vector > result =
          compute_D_prime_A_and_B(CA, SA, CB, SB, mem.get_pos_vkept());
      Vector DA = Super_class::_gt.sub_p(A, mem.get_pos_vkept());
      Vector DB = Super_class::_gt.sub_p(B, mem.get_pos_vkept());
      DApred = std::move(DA);
      DBpred = std::move(DB);
      Vector Dprime = 
        Super_class::_gt.scalar_mult(Super_class::_gt.add_v(result[0], result[1]), 0.5f);

      Point Aprime = Super_class::_gt.add_pv(mem.get_pos_vkept(), 
                       Super_class::_gt.scalar_mult(Dprime, 0.5f));
      Point Bprime = Super_class::_gt.add_pv(mem.get_pos_vkept(), 
                       Super_class::_gt.scalar_mult(Dprime, -0.5f));

      Aprime = Point(std::round(Super_class::_gt.get_x(Aprime)),
                     std::round(Super_class::_gt.get_y(Aprime)),
                     std::round(Super_class::_gt.get_z(Aprime)));

      Bprime = Point(std::round(Super_class::_gt.get_x(Bprime)),
                     std::round(Super_class::_gt.get_y(Bprime)),
                     std::round(Super_class::_gt.get_z(Bprime)));

      if(mem.get_reverse())
      {
        Point save = std::move(Aprime);
        Aprime = std::move(Bprime);
        Bprime = std::move(save);
      }


      Dprime = Super_class::_gt.sub_p(Bprime, Aprime);
      Vector E = Super_class::_gt.sub_v(D, Dprime);
      auto tmp = Super_class::_gt.add_pv(Super_class::_gt.ORIGIN, E);
      E = Vector(std::round(Super_class::_gt.get_x(tmp)), 
                 std::round(Super_class::_gt.get_y(tmp)), 
                 std::round(Super_class::_gt.get_z(tmp)));
      Dpred = D;
      std::vector< Vector > pushed_residual = {E};
      mem.record_error_prediction(pushed_residual);
      return pushed_residual;
    }
    else
    {
      return _delta_for_borders.compute_residuals(mem);
    }
  }

  std::pair< Point, Point > place_points(const std::vector< Vector > &residuals,
                                        vertex_descriptor vkept,
                                        halfedge_descriptor h1,
                                        halfedge_descriptor h2,
                                        bool is_reverse) override
  {
    bool is_border_case = is_on_border(vkept);
    if(!is_border_case)
    {

      const Point& current_point = get(Super_class::_pm, vkept);

      // Case using 1 residual

      Vector E = residuals[0];

      // Fill "a" vertices (vertices around a), count k_a (number of vertices
      // around a) and compute c_a ("crown" around a) vertices.
      std::vector< Point > a_vertices, c_a_vertices;
      k_a = 0;
      fill_around_and_crown(a_vertices, c_a_vertices, h1, h2, k_a);
      // Fill "b" vertices (vertices around b), count k_b (number of vertices
      // around b) and compute c_b ("crown" around b ) vertices.
      std::vector< Point > b_vertices, c_b_vertices;
      k_b = 0;
      fill_around_and_crown(b_vertices, c_b_vertices, h2, h1, k_b);
      // Compute CA, CB (prediction according to crown), SA and SB (prediction
      // according to direct neighbours)
      Point SA = compute_S(a_vertices,
                          get(Super_class::_pm, source(h1, Super_class::_g)),
                          get(Super_class::_pm, source(h2, Super_class::_g)),
                          k_a);

      Point SB = compute_S(b_vertices,
                          get(Super_class::_pm, source(h1, Super_class::_g)),
                          get(Super_class::_pm, source(h2, Super_class::_g)),
                          k_b);

      Point CA = compute_C(c_a_vertices, b_vertices, k_a);
      Point CB = compute_C(c_b_vertices, a_vertices, k_b);

      std::vector< Vector > Ds =
          compute_D_prime_A_and_B(CA, SA, CB, SB, current_point);

      //Vector DprimeA = Ds[0];
      //Vector DprimeB = Ds[1];
      // we get correct DA and DB
      Vector Dprime = 
           Super_class::_gt.scalar_mult(Super_class::_gt.add_v(Ds[0], Ds[1]), 0.5f);

      Point Aprime = Super_class::_gt.add_pv(current_point, Super_class::_gt.scalar_mult(Dprime, 0.5f));
      Point Bprime = Super_class::_gt.add_pv(current_point, Super_class::_gt.scalar_mult(Dprime, -0.5f));

      Aprime = Point(std::round(Super_class::_gt.get_x(Aprime)),
                     std::round(Super_class::_gt.get_y(Aprime)),
                     std::round(Super_class::_gt.get_z(Aprime)));
      Bprime = Point(std::round(Super_class::_gt.get_x(Bprime)),
                     std::round(Super_class::_gt.get_y(Bprime)),
                     std::round(Super_class::_gt.get_z(Bprime)));
      Dprime = Super_class::_gt.sub_p(Bprime, Aprime);

      Vector D = Super_class::_gt.add_v(Dprime, E);

      if(Super_class::_kp->get_type() == VKEPT_POSITION::MIDPOINT)
      {
        Vector DA = Super_class::_gt.scalar_mult(D, -0.5f);
        auto tmp = Super_class::_gt.add_pv(Super_class::_gt.ORIGIN, DA);
        DA = Vector(std::floor(Super_class::_gt.get_x(tmp)),
                    std::floor(Super_class::_gt.get_y(tmp)),
                    std::floor(Super_class::_gt.get_z(tmp)));
        Vector DB = Super_class::_gt.add_v(D, DA);

        Point A = Super_class::_gt.add_pv(current_point, DA);

        Point B = Super_class::_gt.add_pv(current_point, DB);

        std::pair< Point, Point > resulting_points;

        resulting_points = std::make_pair(A, B);
        return resulting_points;
      }
      if(Super_class::_kp->get_type() == VKEPT_POSITION::HALFEDGE)
      {
        Vector residual;
        Point A;
        Point B;

        if(_rev)
        {
          D = Super_class::_gt.scalar_mult(D, -1);
        }
        if(_rev)
        {

          A = Super_class::_gt.add_pv(current_point, D); 
          B = current_point;
        }
        else
        {
          A = current_point;
          B = Super_class::_gt.add_pv(current_point, D); 
        }


        std::pair< Point, Point > resulting_points;
        resulting_points = std::make_pair(A, B);
        return resulting_points;
      }
      return std::make_pair(Point(0, 0, 0), Point(0, 0, 0));
    }
    else
    {
      return _delta_for_borders.place_points(residuals, 
                                            vkept, 
                                            h1, 
                                            h2, 
                                            is_reverse);    
	}
  }

  const std::tuple< bool, bool, bool, bool >& get_midpoint() const
  {
    return _round_midpoint;
  }
  void set_bit_info(bool b1, bool b2, bool b3, bool b4)
  {
    _round_midpoint = std::make_tuple(b1, b2, b3, b4);
  }
  void set_rev(bool b) override 
  {
	  _rev = b; _delta_for_borders.set_rev(b);
  }

  std::string get_as_string() const override { return "butterfly"; }


  // used for debug: to record the stencil around a vertex, in order to color
  // them later
  void record_stencil(std::vector< vertex_descriptor > &around,
                      std::vector< vertex_descriptor > &crown,
                      halfedge_descriptor begin,
                      halfedge_descriptor end,
                      int &k)
  {
    halfedge_descriptor begin_b =
        opposite(next(begin, Super_class::_g), Super_class::_g);
    halfedge_descriptor first_crown_b = next(
        opposite(prev(begin, Super_class::_g), Super_class::_g), Super_class::_g);


    k = 0;
    crown.push_back(target(first_crown_b, Super_class::_g));
    while(begin_b != end)
    {
      around.push_back(source(begin_b, Super_class::_g));
      halfedge_descriptor current_crown =
          next(opposite(prev(begin_b, Super_class::_g), Super_class::_g),
               Super_class::_g);
      crown.push_back(target(current_crown, Super_class::_g));
      begin_b = opposite(next(begin_b, Super_class::_g), Super_class::_g);
      k++;
    }
  }


private:
  bool is_on_border(vertex_descriptor v) const 
  {
    boost::iterator_range<
        CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
        iterator_range = CGAL::halfedges_around_target(v, Super_class::_g);
    for(auto h_v : iterator_range)
    {
      if(CGAL::is_border(source(h_v, Super_class::_g), Super_class::_g))
      {
        return true;
      }
    }
    return false;
  }

  std::vector< Vector > compute_D_prime_A_and_B(const Point &CA,
                                           const Point &SA,
                                           const Point &CB,
                                           const Point &SB,
                                           const Point &kept_position)
  {

    // Convert SA,SB,CA,CB to Vector (easier for computation)
    Vector SA_vec = Vector(Super_class::_gt.get_x(SA),
                           Super_class::_gt.get_y(SA),
                           Super_class::_gt.get_z(SA));
    Vector SB_vec = Vector(Super_class::_gt.get_x(SB),
                           Super_class::_gt.get_y(SB),
                           Super_class::_gt.get_z(SB));


    Vector CA_vec = Vector(Super_class::_gt.get_x(CA),
                           Super_class::_gt.get_y(CA),
                           Super_class::_gt.get_z(CA));
    Vector CB_vec = Vector(Super_class::_gt.get_x(CB),
                           Super_class::_gt.get_y(CB),
                           Super_class::_gt.get_z(CB));

    Vector V = Vector(Super_class::_gt.get_x(kept_position),
                      Super_class::_gt.get_y(kept_position),
                      Super_class::_gt.get_z(kept_position));


    // following the formula from the paper

    Vector alphaSA = Super_class::_gt.scalar_mult(SA_vec, _alpha);
    Vector alphaSB = Super_class::_gt.scalar_mult(SB_vec, _alpha);

    Vector opp_alphaCA = Super_class::_gt.scalar_mult(CA_vec, (1 - _alpha));
    Vector opp_alphaCB = Super_class::_gt.scalar_mult(CB_vec, (1 - _alpha));

    Vector left_part_A =
        Super_class::_gt.scalar_mult(Super_class::_gt.add_v(opp_alphaCA, alphaSA), k_a + 3);
    Vector left_part_B =
        Super_class::_gt.scalar_mult(Super_class::_gt.add_v(opp_alphaCB, alphaSB), k_b + 3);

    Vector right_part_A = Super_class::_gt.scalar_mult(V, _alpha - k_a - 3);
    Vector right_part_B = Super_class::_gt.scalar_mult(V, _alpha - k_b - 3);

    Vector DA_prime = Super_class::_gt.add_v(left_part_A, right_part_A);
    Vector DB_prime = Super_class::_gt.add_v(left_part_B, right_part_B);

    DA_prime = Super_class::_gt.scalar_mult(DA_prime, 1 / (k_a + 3 + _alpha));
    DB_prime = Super_class::_gt.scalar_mult(DB_prime, 1 / (k_b + 3 + _alpha));
    auto tmp  = Super_class::_gt.add_pv(Super_class::_gt.ORIGIN, DA_prime);
    auto tmp2 = Super_class::_gt.add_pv(Super_class::_gt.ORIGIN, DB_prime);

    DA_prime = Vector(std::round(Super_class::_gt.get_x(tmp) * (-2)),
                      std::round(Super_class::_gt.get_y(tmp) * (-2)),
                      std::round(Super_class::_gt.get_z(tmp) * (-2)));
    DB_prime = Vector(std::round(Super_class::_gt.get_x(tmp2) * (2)),
                      std::round(Super_class::_gt.get_y(tmp2) * (2)),
                      std::round(Super_class::_gt.get_z(tmp2) * (2)));

    std::vector< Vector > Ds = {DA_prime, DB_prime};
    return Ds;
  }

  FEVV::Filters::Delta_predictor< HalfedgeGraph,
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
        opposite(next(begin, Super_class::_g), Super_class::_g);
    halfedge_descriptor first_crown_b = next(
        opposite(prev(begin, Super_class::_g), Super_class::_g), Super_class::_g);


    k = 0;
    crown.push_back(
        get(Super_class::_pm, target(first_crown_b, Super_class::_g)));
    while(begin_b != end)
    {
      around.push_back(get(Super_class::_pm, source(begin_b, Super_class::_g)));
      halfedge_descriptor current_crown =
          next(opposite(prev(begin_b, Super_class::_g), Super_class::_g),
               Super_class::_g);
      crown.push_back(
          get(Super_class::_pm, target(current_crown, Super_class::_g)));
      begin_b = opposite(next(begin_b, Super_class::_g), Super_class::_g);
      k++;
    }
  }

  Point compute_S(std::vector< Point > &around, Point v1, Point v2, int &k)
  {
    // Compute S
    Point S(0, 0, 0);
    for(auto it = around.begin(); it != around.end(); it++)
    {
      S = Point(Super_class::_gt.get_x(S) + Super_class::_gt.get_x(*it),
                Super_class::_gt.get_y(S) + Super_class::_gt.get_y(*it),
                Super_class::_gt.get_z(S) + Super_class::_gt.get_z(*it));
    }

    S = Point(Super_class::_gt.get_x(S) + Super_class::_gt.get_x(v1) +
                  Super_class::_gt.get_x(v2),
              Super_class::_gt.get_y(S) + Super_class::_gt.get_y(v1) +
                  Super_class::_gt.get_y(v2),
              Super_class::_gt.get_z(S) + Super_class::_gt.get_z(v1) +
                  Super_class::_gt.get_z(v2));
    S = Point(Super_class::_gt.get_x(S) / static_cast< double >((k + 3)),
              Super_class::_gt.get_y(S) / static_cast< double >((k + 3)),
              Super_class::_gt.get_z(S) / static_cast< double >((k + 3)));
    return S;
  }

  Point compute_C(std::vector< Point > &crown,
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
    Point C = Point(Super_class::_gt.get_x(b_one) + Super_class::_gt.get_x(b_kb),
                    Super_class::_gt.get_y(b_one) + Super_class::_gt.get_y(b_kb),
                    Super_class::_gt.get_z(b_one) + Super_class::_gt.get_z(b_kb));

    for(auto it = crown.begin(); it != crown.end(); it++)
    {
      C = Point(Super_class::_gt.get_x(C) + Super_class::_gt.get_x(*it),
                Super_class::_gt.get_y(C) + Super_class::_gt.get_y(*it),
                Super_class::_gt.get_z(C) + Super_class::_gt.get_z(*it));
    }
    C = Point(Super_class::_gt.get_x(C) / static_cast< double >((k + 3)),
              Super_class::_gt.get_y(C) / static_cast< double >((k + 3)),
              Super_class::_gt.get_z(C) / static_cast< double >((k + 3)));
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
        iterator_range = CGAL::halfedges_around_target(center, Super_class::_g);
    for(auto h_v : iterator_range)
    {
      around.push_back(source(h_v, Super_class::_g));

      halfedge_descriptor next_h =
          next(opposite(h_v, Super_class::_g), Super_class::_g);

      halfedge_descriptor opp_next_h = opposite(next_h, Super_class::_g);


      if(!CGAL::is_border(opp_next_h, Super_class::_g) &&
         !CGAL::is_border(
             next_h, Super_class::_g)) // second verification might be overkill
      {

        halfedge_descriptor next_to_crown = next(opp_next_h, Super_class::_g);
        crown.push_back(target(next_to_crown, Super_class::_g));
      }
      k += 1;
    }
  }
};
} // namespace Filters
} // namespace FEVV
