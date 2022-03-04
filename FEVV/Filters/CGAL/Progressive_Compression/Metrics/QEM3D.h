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

#include "ErrorMetric.h"

#include <Eigen/Dense>

#include <map>
#include <tuple>

namespace FEVV {
namespace Filters {
template< typename HalfedgeGraph,
          typename PointMap >
class QEM3D : public ErrorMetric< HalfedgeGraph,
                                  PointMap >
{
public:
  using vertex_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor;
  using halfedge_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor;
  using edge_iterator =
      typename boost::graph_traits< HalfedgeGraph >::edge_iterator;
  using edge_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::edge_descriptor;
  using face_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::face_descriptor;
  using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;
  using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
  using Geometry = typename FEVV::Geometry_traits< HalfedgeGraph >;

  typedef ErrorMetric< HalfedgeGraph,
                       PointMap >
      SuperClass;

  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::VectorXd VectorX;
  typedef typename std::tuple< Matrix, VectorX, double > Quadric;


  QEM3D(
      HalfedgeGraph &g,
      PointMap &pm,
      KeptPosition< HalfedgeGraph, PointMap >
          *vkept,
      FEVV::Filters::UniformDequantization< HalfedgeGraph, PointMap >
          &dequantiz)
      : SuperClass(g, pm, vkept, dequantiz)
  {
  }
  ~QEM3D(){}

  double ComputeAreaTriangle(const Point &p1, const Point &p2, const Point &p3) const
  {
    auto a2 = SuperClass::_gt.length2(SuperClass::_gt.sub_p(p2, p1));
    auto b2 = SuperClass::_gt.length2(SuperClass::_gt.sub_p(p3, p2));
    auto c2 = SuperClass::_gt.length2(SuperClass::_gt.sub_p(p1, p3));
    double area2 = (2*(a2*b2 + b2*c2 + a2*c2) - (a2*a2 + b2*b2 + c2*c2))*0.0625f; // expanded Heron's formula
    return sqrt(area2);
  }

  double ComputeAreaTriangleAfterDequantization(face_descriptor f) const
  {
    halfedge_descriptor h1 = halfedge(f, SuperClass::_g);
    halfedge_descriptor h2 = next(h1, SuperClass::_g);
    const Point& p1 = get(SuperClass::_pm, source(h1, SuperClass::_g));
    const Point& p2 = get(SuperClass::_pm, target(h1, SuperClass::_g));
    const Point& p3 = get(SuperClass::_pm, target(h2, SuperClass::_g));
    Point p1_tmp = SuperClass::_dequantiz.dequantize_vertex(p1);
    Point p2_tmp = SuperClass::_dequantiz.dequantize_vertex(p2);
    Point p3_tmp = SuperClass::_dequantiz.dequantize_vertex(p3);

    return ComputeAreaTriangle(p1_tmp, p2_tmp, p3_tmp);
  }

  virtual void ComputeError() override
  {
    if (!SuperClass::_queue.empty())
    {
      typename SuperClass::priority_queue_edges empty2;
      std::swap(SuperClass::_queue, empty2);
    }
    SuperClass::_edges_cost.clear();
    std::swap(SuperClass::_edges_cost, SuperClass::edge2cost_map());
    if (SuperClass::_edges_cost.empty())
    {
      faces_quadrics.clear();
      auto face_pair = faces(SuperClass::_g);
      auto face_it = face_pair.first;
      auto face_it_end = face_pair.second;
      // Compute every face quadric
      for (; face_it != face_it_end; ++face_it)
      {
        Quadric curr_quadric = ComputeFaceQuadric(*face_it);
        faces_quadrics.emplace(*face_it, curr_quadric);
      }

      auto vertices_pair = vertices(SuperClass::_g);
      auto vertices_it = vertices_pair.first;
      auto vertices_it_end = vertices_pair.second;
      // Compute every vertex quadric
      for (; vertices_it != vertices_it_end; ++vertices_it)
      {
        Quadric curr_vertex_quadric = ComputeVertexQuadric(*vertices_it);
        vertices_quadrics.emplace(*vertices_it, curr_vertex_quadric);
      }
      SuperClass::_pm = get(boost::vertex_point, SuperClass::_g);

      auto edge_iterator_pair = edges(SuperClass::_g);
      auto edge_ite = edge_iterator_pair.first;
      int count = 0;
      SuperClass::_threshold = 0;
      for (; edge_ite != edge_iterator_pair.second;
        ++edge_ite)
      {
        Point collapsePos;
        collapsePos = SuperClass::_vkept->ComputePosition(*edge_ite);
        double weight = ComputeCostEdge(*edge_ite, collapsePos);

        SuperClass::_threshold += weight;
        ++count;
        ///////////////////////////////////////////////////////////////////////
        SuperClass::_queue.push(
          std::make_tuple(*edge_ite, weight, collapsePos));
        SuperClass::_edges_cost.emplace(
          *edge_ite, std::make_pair(weight, collapsePos));
      }
      SuperClass::_threshold /= count;
    }
  }


  double ComputeCostEdge(edge_descriptor e, const Point &collapsePos) override
  {
    auto current_halfedge = halfedge(e, SuperClass::_g);

    vertex_descriptor v1 = target(current_halfedge, SuperClass::_g);
    vertex_descriptor v2 = source(current_halfedge, SuperClass::_g);
    Point dequantized_collapsePos =
        SuperClass::_dequantiz.dequantize_vertex(collapsePos);
    Eigen::Vector4d v =
        Eigen::Vector4d(SuperClass::_gt.get_x(dequantized_collapsePos),
                        SuperClass::_gt.get_y(dequantized_collapsePos),
                        SuperClass::_gt.get_z(dequantized_collapsePos),
                        1);
    auto vt = v.transpose();

    Matrix Q(4, 4);
    Quadric Q1 = vertices_quadrics[v1];
    Quadric Q2 = vertices_quadrics[v2];
    Matrix A1 = std::get< 0 >(Q1);
    Matrix A2 = std::get< 0 >(Q2);
    VectorX B1 = std::get< 1 >(Q1);
    VectorX B2 = std::get< 1 >(Q2);
    double d1 = std::get< 2 >(Q1);
    double d2 = std::get< 2 >(Q2);
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        Q(i, j) = A1(i, j) + A2(i, j);
      }
    }
    for(int k = 0; k < 3; k++)
    {
      Q(3, k) = B1(k) + B2(k);
      Q(k, 3) = B1(k) + B2(k);
    }
    Q(3, 3) = d1 + d2;
    Matrix VtQ = vt * Q;
    Matrix VtQV = VtQ * v;
    double cost = VtQV(0, 0);
    return cost;
  }

  double ComputeCostEdge(edge_descriptor e)
  {
    std::cout << "should not happen" << std::endl;
    return 0.0;
  }

  std::string getMethodasString() const override { return "QEM3D"; }

protected:
  std::map< face_descriptor, Quadric > faces_quadrics;

  std::map< vertex_descriptor, Quadric > vertices_quadrics;

  std::vector< double > getPlaneEquation(face_descriptor f) const
  {
    halfedge_descriptor h1 = halfedge(f, SuperClass::_g);
    halfedge_descriptor h2 = next(h1, SuperClass::_g);

    // get points a the face and dequantize them
    const Point& p1 = get(SuperClass::_pm, source(h1, SuperClass::_g));
    const Point& p2 = get(SuperClass::_pm, target(h1, SuperClass::_g));
    const Point& p3 = get(SuperClass::_pm, target(h2, SuperClass::_g));

    Point p1_tmp = SuperClass::_dequantiz.dequantize_vertex(p1);
    Point p2_tmp = SuperClass::_dequantiz.dequantize_vertex(p2);
    Point p3_tmp = SuperClass::_dequantiz.dequantize_vertex(p3);


    auto x1 = SuperClass::_gt.get_x(p1_tmp);
    auto y1 = SuperClass::_gt.get_y(p1_tmp);
    auto z1 = SuperClass::_gt.get_z(p1_tmp);

    auto x2 = SuperClass::_gt.get_x(p2_tmp);
    auto y2 = SuperClass::_gt.get_y(p2_tmp);
    auto z2 = SuperClass::_gt.get_z(p2_tmp);

    auto x3 = SuperClass::_gt.get_x(p3_tmp);
    auto y3 = SuperClass::_gt.get_y(p3_tmp);
    auto z3 = SuperClass::_gt.get_z(p3_tmp);

    // compute equation of a plane: compute normal vector to the plane
    Eigen::Vector3d v1;
    v1 << x2 - x1, y2 - y1, z2 - z1;
    Eigen::Vector3d v2;
    v2 << x3 - x1, y3 - y1, z3 - z1;
    // cross product computes the normal vector to the plane
    Eigen::Vector3d cp = v1.cross(v2);

    // normalize the vector
    cp.normalize();
    if(cp.hasNaN())
    {
      cp = Eigen::Vector3d(0, 0, 0);
    }
    double a = cp[0];
    double b = cp[1];
    double c = cp[2];
    double d = -(a * x1 + b * y1 + c * z1);
    std::vector< double > plane_eq;
	plane_eq.reserve(4);
    plane_eq.push_back(a);
    plane_eq.push_back(b);
    plane_eq.push_back(c);
    plane_eq.push_back(d);
    return plane_eq;
  }


  virtual Quadric ComputeFaceQuadric(face_descriptor f) const
  {
    // the quadric of a face is a representation of the equation of this place
    std::vector< double > plane_eq = getPlaneEquation(f);
    Eigen::Vector3d n(plane_eq[0], plane_eq[1], plane_eq[2]);

    double d = plane_eq[3];

    Matrix Af(3, 3);
    VectorX Bf(3);
    double Cf = d * d;

    for(size_t i = 0; i < 3; i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        Af(i, j) = n(i) * n(j);
      }
    }
    Bf[0] = d * n(0);
    Bf[1] = d * n(1);
    Bf[2] = d * n(2);

    /* the quadric takes the form as a matrix:
    |a^2 ab ac  ad|
    |ab b^2 bc  bd|
    |ac  bc c^2 cd|
    |ad  bd cd d^2|
    A = |a^2 ab ac |
            |ab b^2 bc |
            |ac  bc c^2|
    B = |ad  bd cd |
    C = d^2
    */


    return std::make_tuple(Af, Bf, Cf);
  }

  virtual Quadric ComputeVertexQuadric(vertex_descriptor v) const
  {
    // The quadric of a vertex is the sum of the quadrics of its adjacent faces.
    boost::iterator_range<
        CGAL::Halfedge_around_target_iterator< HalfedgeGraph > >
        iterator_range = CGAL::halfedges_around_target(v, SuperClass::_g);

    Matrix M(3, 3);
    VectorX V(3);
    double D = 0;
    M << 0, 0, 0, 0, 0, 0, 0, 0, 0;
    V << 0, 0, 0;
    // get every adjacent quadric
    for(auto h_v : iterator_range)
    {
      face_descriptor f = face(h_v, SuperClass::_g);
      // avoid border case
      if(f != boost::graph_traits< HalfedgeGraph >::null_face())
      {
        auto it_curr_quadric = faces_quadrics.find(f);
        if(it_curr_quadric != faces_quadrics.end())
        {
          // if a quadric has been found (no border case)
          double area = ComputeAreaTriangleAfterDequantization(f);
          M = M + area * std::get< 0 >(it_curr_quadric->second);
          V = V + area * std::get< 1 >(it_curr_quadric->second);
          D = D + area * std::get< 2 >(it_curr_quadric->second);
        }
      }
      if(f == boost::graph_traits< HalfedgeGraph >::null_face() ||
         CGAL::is_border(opposite(h_v, SuperClass::_g), SuperClass::_g))
      {
        // border case: create a perpendicular plane
        halfedge_descriptor h_v_opp = opposite(h_v, SuperClass::_g);
        face_descriptor f_opp;
        if(f == boost::graph_traits< HalfedgeGraph >::null_face())
        {
          f_opp = face(h_v_opp, SuperClass::_g);
        }
        else
        {
          f_opp = face(h_v, SuperClass::_g);
        }
        
        std::vector< double > plane_eq_opp = getPlaneEquation(f_opp);
        Eigen::Vector3d n(plane_eq_opp[0], plane_eq_opp[1], plane_eq_opp[2]);

        const Point& p1 = get(SuperClass::_pm, source(h_v_opp, SuperClass::_g));
        const Point& p2 = get(SuperClass::_pm, target(h_v_opp, SuperClass::_g));

        Point p1_tmp = SuperClass::_dequantiz.dequantize_vertex(p1);
        Point p2_tmp = SuperClass::_dequantiz.dequantize_vertex(p2);
        double area = SuperClass::_gt.length2(SuperClass::_gt.sub_p(p2_tmp, p1_tmp)) * std::sqrt(3)/2.; // equilateral triangle area x 2

        auto x1 = SuperClass::_gt.get_x(p1_tmp);
        auto y1 = SuperClass::_gt.get_y(p1_tmp);
        auto z1 = SuperClass::_gt.get_z(p1_tmp);

        auto x2 = SuperClass::_gt.get_x(p2_tmp);
        auto y2 = SuperClass::_gt.get_y(p2_tmp);
        auto z2 = SuperClass::_gt.get_z(p2_tmp);

        Eigen::Vector3d v1;
        v1 << x2 - x1, y2 - y1, z2 - z1;
        Eigen::Vector3d v2;
        v2 << x1 + n(0), y1 + n(1), z1 + n(2);
        // equation of perpendicular plane
        Eigen::Vector3d perpendicular_equation = v1.cross(v2);

        perpendicular_equation.normalize();
        if(perpendicular_equation.hasNaN())
        {
          perpendicular_equation = Eigen::Vector3d(0, 0, 0);
        }
        double a = perpendicular_equation[0];
        double b = perpendicular_equation[1];
        double c = perpendicular_equation[2];
        double d = -(a * x1 + b * y1 + c * z1);

        Matrix Af(3, 3);
        VectorX Bf(3);
        double Cf = d * d;

        for(size_t i = 0; i < 3; i++)
        {
          for(size_t j = 0; j < 3; j++)
          {
            Af(i, j) = perpendicular_equation(i) * perpendicular_equation(j);
          }
        }
        Bf[0] = d * perpendicular_equation(0);
        Bf[1] = d * perpendicular_equation(1);
        Bf[2] = d * perpendicular_equation(2);
        M = M + area * Af;
        V = V + area * Bf;
        D = D + area * Cf;
      }
    }
    return std::make_tuple(M, V, D);
  }
};
} // namespace Filters
} // namespace FEVV
