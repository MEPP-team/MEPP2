// Copyright (c) 2012-2020 University of Lyon and CNRS (France).
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

#include <set>
#include <stack>

#include <Eigen/Dense> // for Matrix::Identity()
#include <cmath> // for std::isnan()

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/boost/graph/iterator.h> // for circulators
#include <CGAL/boost/graph/helpers.h>  // for is_border_edge()

#include "FEVV/Filters/Generic/Manifold/Curvature/eigen_val_vect.hpp"

#include "FEVV/Wrappings/Geometry_traits.h"

#include "FEVV/Tools/Math/MatrixOperations.hpp"
#include "FEVV/Operators/Geometry/AngleOperations.hpp"
#include "FEVV/Operators/Geometry/ClippingAndIntersection.hpp"

#include "FEVV/Operators/Geometry/triangles.hpp"

namespace FEVV {
namespace Filters {

template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
struct v_Curv_cmdm
{
  double KminCurv; ///< The minimum curvature
  double KmaxCurv; ///< The minimum curvature

  typename GeometryTraits::Vector
      VKminCurv; ///< The minimum curvature direction
  typename GeometryTraits::Vector
      VKmaxCurv; ///< The minimum curvature direction
};

template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
double
fabs_cmdm(const v_Curv_cmdm< HalfedgeGraph > &input)
{
  return std::max(::fabs(input.KminCurv), ::fabs(input.KmaxCurv));
}

template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
double
triangle_area_cmdm(
    const HalfedgeGraph &g,
    const PointMap &pm,
    typename boost::graph_traits< HalfedgeGraph >::face_descriptor fd,
    const GeometryTraits &gt)
{
  typedef typename boost::property_traits< PointMap >::value_type Point3d;
  using Vector = typename GeometryTraits::Vector;

  CGAL::Halfedge_around_face_circulator< HalfedgeGraph > p_halfedge(
      halfedge(fd, g), g);
  Point3d p = get(pm, target(*p_halfedge, g));
  Point3d q = get(pm, target(next(*p_halfedge, g), g));
  Point3d r = get(pm, target(next(next(*p_halfedge, g), g), g));

  return FEVV::Operators::Geometry::triangle_area< GeometryTraits >(
      p, q, r, gt);
}

//********************************************
// 1-ring neighborhood curvature for a vertex
//********************************************
template< typename HalfedgeGraph,
          typename PointMap,
          typename FaceNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >,
          typename GraphTraits = boost::graph_traits< HalfedgeGraph >,
          typename VertexDescriptor = typename GraphTraits::vertex_descriptor >
void
principal_curvature_per_vert_cmdm(const HalfedgeGraph &g,
                                  const PointMap &pm,
                                  const FaceNormalMap &f_nm,
                                  VertexDescriptor vd,
                                  double pp_matrix_sum[3][3],
                                  const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Point Point3d;
  using Vector = typename GeometryTraits::Vector;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  double area = 0;

  // iterate over all edges
  CGAL::Halfedge_around_target_circulator< HalfedgeGraph > p_halfedge(vd, g),
      done(p_halfedge);
  do
  {
    // build edge vector and comput its norm
    Point3d p1 = get(pm, target(*p_halfedge, g));
    Point3d p2 = get(pm, source(*p_halfedge, g));
#if 1 // WAITING for AIF stuff - TODO : REMOVE THIS BLOC LATER...
    Vector edge = gt.sub_p(p1, p2);
#else
    Vector edge(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
#endif
    double len_edge = gt.length(edge);
    if(len_edge == 0) // avoid divide by zero
      continue;

    // compute (signed) angle between two incident faces, if exists
    face_descriptor p_facet1 = face(*p_halfedge, g);
    face_descriptor p_facet2 = face(opposite(*p_halfedge, g), g);
    assert(p_facet1 != p_facet2);
    if(p_facet1 == GraphTraits::null_face() ||
       p_facet2 == GraphTraits::null_face())
      continue; // border edge

    Vector normal1 = f_nm[p_facet1];
    Vector normal2 = f_nm[p_facet2];

    if(std::isnan(normal1[0]) || std::isnan(normal1[1]) ||
       std::isnan(normal1[2]) ||
       std::isnan(normal2[0]) || std::isnan(normal2[1]) ||
       std::isnan(normal2[2]))
      continue; // Normal vector is set to be (NaN, NaN, NaN) if 2 edges of
                // the face are collinear -- [YN]

    area += triangle_area_cmdm(g, pm, p_facet1, gt);

    // ---
    Vector cp(normal1[1] * normal2[2] - normal1[2] * normal2[1],
              normal1[2] * normal2[0] - normal1[0] * normal2[2],
              normal1[0] * normal2[1] - normal1[1] * normal2[0]);
    double ps_cp_xedge = cp[0] * edge[0] + cp[1] * edge[1] + cp[2] * edge[2];

    // ---
    double sine =
        ps_cp_xedge /
        len_edge; // = (CGAL::cross_product(normal1,normal2)*edge)/len_edge;
                  // TODO BETTER
    double beta = FEVV::Operators::Geometry::asin(sine);

    // compute edge * edge^t * coeff, and add it to current matrix
    double p_vector_edge[3] = {edge[0], edge[1], edge[2]};
    double pp_matrix[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    FEVV::Math::Matrix::Square::vector_times_transpose_mult(
        p_vector_edge, pp_matrix, beta / len_edge);
    FEVV::Math::Matrix::Square::add(pp_matrix, pp_matrix_sum);
  } while(++p_halfedge != done);

  if(area != 0)
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
        pp_matrix_sum[i][j] /= area;
  else
    std::cout << "Unable to compute curvature." << std::endl;
}

//**********************************************
// geodesic neighborhood curvature for a vertex
//**********************************************
template< typename HalfedgeGraph,
          typename PointMap,
          typename FaceNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >,
          typename GraphTraits = boost::graph_traits< HalfedgeGraph >,
          typename VertexDescriptor = typename GraphTraits::vertex_descriptor >
void
geodes_principal_curvature_per_vert_cmdm(const HalfedgeGraph &g,
                                         const PointMap &pm,
                                         const FaceNormalMap &f_nm,
                                         VertexDescriptor vd,
                                         double pp_matrix_sum[3][3],
                                         double radius,
                                         const GeometryTraits &gt)
{
  typedef typename GeometryTraits::Point Point3d;
  using Vector = typename GeometryTraits::Vector;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  std::set< VertexDescriptor > vertices;
  Point3d o = get(pm, vd);
  std::stack< VertexDescriptor > s;
  s.push(vd);
  vertices.insert(vd);
  int iter = 0;
  while(!s.empty())
  {
    VertexDescriptor v = s.top();
    s.pop();
    Point3d p = get(pm, v);
    CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h(v, g), done(h);
    do
    {
      Point3d p1 = get(pm, target(*h, g));
      Point3d p2 = get(pm, source(*h, g));

#if 1 // WAITING for AIF stuff - TODO : REMOVE THIS BLOC LATER...
      Vector vec = gt.sub_p(p2, p1);
      Vector pm_o = gt.sub_p(p, o);
#else

      Vector vec(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
      Vector PmO(P[0] - O[0], P[1] - O[1], P[2] - O[2]);
#endif
      double ps_vx_pm_o =
          vec[0] * pm_o[0] + vec[1] * pm_o[1] + vec[2] * pm_o[2];

      // ---
      if(v == vd || ps_vx_pm_o >= 0.0) // if (v==pVertex || vec * (P - O) >=
                                       // 0.0) // TODO BETTER ; [>"=" (YN)]
      {
        bool normal1_nan = false, normal2_nan = false;               // [YN]
        bool pt_face1_superpose = false, pt_face2_superpose = false; // [YN]

        bool isect = FEVV::Operators::Geometry::sphere_clip_vector(
            o, radius, p, vec, gt);

        if(!CGAL::is_border_edge(*h, g) && isect) // ["&& isect" (YN)]
        {
          double len_edge =
              std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] +
                        vec[2] * vec[2]); // std::sqrt(vec*vec);  TODO BETTER

          // compute (signed) angle between two incident faces, if exists
          face_descriptor p_facet1 = face(*h, g);
          face_descriptor p_facet2 = face(opposite(*h, g), g);
          assert(p_facet1 != p_facet2);

          std::vector< Point3d > face1_vertices; // [YN]
          BOOST_FOREACH(VertexDescriptor vi, CGAL::vertices_around_face(*h, g))
            face1_vertices.push_back(get(pm, vi));

          if(face1_vertices[0][0] == face1_vertices[1][0] &&
             face1_vertices[1][0] == face1_vertices[2][0] &&
             face1_vertices[0][0] == face1_vertices[2][0] &&
             face1_vertices[0][1] == face1_vertices[1][1] &&
             face1_vertices[1][1] == face1_vertices[2][1] &&
             face1_vertices[0][1] == face1_vertices[2][1] &&
             face1_vertices[0][2] == face1_vertices[1][2] &&
             face1_vertices[1][2] == face1_vertices[2][2] &&
             face1_vertices[0][2] == face1_vertices[2][2])
            pt_face1_superpose = true;

          std::vector< Point3d > face2_vertices; // [YN]
          BOOST_FOREACH(VertexDescriptor vi,
                        CGAL::vertices_around_face(opposite(*h, g), g))
            face2_vertices.push_back(get(pm, vi));

          if(face2_vertices[0][0] == face2_vertices[1][0] &&
             face2_vertices[1][0] == face2_vertices[2][0] &&
             face2_vertices[0][0] == face2_vertices[2][0] &&
             face2_vertices[0][1] == face2_vertices[1][1] &&
             face2_vertices[1][1] == face2_vertices[2][1] &&
             face2_vertices[0][1] == face2_vertices[2][1] &&
             face2_vertices[0][2] == face2_vertices[1][2] &&
             face2_vertices[1][2] == face2_vertices[2][2] &&
             face2_vertices[0][2] == face2_vertices[2][2])
            pt_face2_superpose = true;

          if(p_facet1 != GraphTraits::null_face() &&
             p_facet2 != GraphTraits::null_face())
          {
            Vector normal1 = f_nm[p_facet1];
            Vector normal2 = f_nm[p_facet2];

            // Check if face1 and face2 contain collinear edges [YN]
            if(std::isnan(normal1[0]) || std::isnan(normal1[1]) ||
               std::isnan(normal1[2]))
              normal1_nan = true;
            if(std::isnan(normal2[0]) || std::isnan(normal2[1]) ||
               std::isnan(normal2[2]))
              normal2_nan = true;

            // ---
            double beta = 0;

            if(!pt_face2_superpose && !pt_face1_superpose)
            {
              // Check for common cases encountered in a mesh distorted by a
              // uniform geometric quantization [YN]
              if(!normal1_nan && !normal2_nan)
              {
                // Case 1: the 2 normals are colinear in opposite direction
                double ps = normal1[0] * normal2[0] + normal1[1] * normal2[1] +
                            normal1[2] * normal2[2]; // scalar product
                double acos = ps / std::sqrt((normal1[0] * normal1[0] +
                                              normal1[1] * normal1[1] +
                                              normal1[2] * normal1[2]) *
                                             (normal2[0] * normal2[0] +
                                              normal2[1] * normal2[1] +
                                              normal2[2] * normal2[2]));

                if(acos == -1) // cos=180 degrees => the 2 normals oppose
                  double beta =
                      M_PI /
                      2; // set to maximum curvature achievable (in radian)

                else
                {
                  // Normal and most common case for all types of distortion
                  Vector cp(normal1[1] * normal2[2] - normal1[2] * normal2[1],
                            normal1[2] * normal2[0] - normal1[0] * normal2[2],
                            normal1[0] * normal2[1] - normal1[1] * normal2[0]);
                  double ps_cpx_v =
                      cp[0] * vec[0] + cp[1] * vec[1] + cp[2] * vec[2];

                  double sine =
                      ps_cpx_v /
                      len_edge; // (CGAL::cross_product(normal1, normal2)*vec)
                                // / len_edge; // TODO BETTER
                  beta = FEVV::Operators::Geometry::asin(sine);
                }
              }
              else if((normal1_nan && !normal2_nan) ||
                      (normal2_nan && !normal1_nan))
              {
                // Case 2: face2 contains collinear edges
                if(!normal1_nan)
                {
                  // consider all haledge of the source vertex
                  bool found_face = false;
                  Point3d destination = get(pm, target(*h, g));
                  VertexDescriptor t1 = source(*h, g);

                  std::stack< VertexDescriptor > s2;
                  std::set< VertexDescriptor > s22;
                  s22.insert(t1);
                  s2.push(t1);
                  while(!s2.empty() && !found_face)
                  {
                    VertexDescriptor t = s2.top();
                    s2.pop();
                    CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h2(
                        t, g),
                        done2(h2);
                    do
                    {
                      if(get(pm, source(*h2, g)) == destination)
                      {
                        face_descriptor newface2 = face(*h2, g);
                        if(newface2 != GraphTraits::null_face())
                        {
                          normal2 = f_nm[newface2];
                          if(!std::isnan(normal2[0]) &&
                             !std::isnan(normal2[1]) &&
                             !std::isnan(normal2[2]) && newface2 != p_facet1)
                          {
                            Vector cp(normal1[1] * normal2[2] -
                                          normal1[2] * normal2[1],
                                      normal1[2] * normal2[0] -
                                          normal1[0] * normal2[2],
                                      normal1[0] * normal2[1] -
                                          normal1[1] * normal2[0]);
                            double ps_cpx_v = cp[0] * vec[0] + cp[1] * vec[1] +
                                              cp[2] * vec[2];

                            double sine = ps_cpx_v / len_edge;
                            beta = FEVV::Operators::Geometry::asin(sine);
                            found_face = true;
                            break;
                          }
                          else if(s22.find(source(*h2, g)) == s22.end())
                          {
                            s2.push(source(*h2, g));
                            s22.insert(source(*h2, g));
                          }
                        }
                      }
                      h2++;
                    } while(h2 != done2);
                    destination = get(pm, source(*h, g));
                  }
                }

                // Case 3: face1 contains collinear edges
                if(!normal2_nan)
                {
                  // consider all haledge of the source vertex
                  Point3d destination = get(pm, target(*h, g));
                  bool found_face = false;
                  VertexDescriptor v1 = source(*h, g);

                  std::stack< VertexDescriptor > s1;
                  std::set< VertexDescriptor > s11;
                  s11.insert(v1);
                  s1.push(v1);
                  while(!s1.empty() && !found_face)
                  {
                    VertexDescriptor s = s1.top();
                    s1.pop();
                    CGAL::Halfedge_around_target_circulator< HalfedgeGraph > h1(
                        s, g),
                        done1(h1);
                    do
                    {
                      if(get(pm, source(*h1, g)) == destination)
                      {
                        face_descriptor newface1 = face(*h1, g);
                        if(newface1 != GraphTraits::null_face())
                        {
                          normal1 = f_nm[newface1]; // face or opposite face??
                          if(!std::isnan(normal1[0]) &&
                             !std::isnan(normal1[1]) &&
                             !std::isnan(normal1[2]) && newface1 != p_facet2)
                          {
                            Vector cp(normal1[1] * normal2[2] -
                                          normal1[2] * normal2[1],
                                      normal1[2] * normal2[0] -
                                          normal1[0] * normal2[2],
                                      normal1[0] * normal2[1] -
                                          normal1[1] * normal2[0]);
                            double ps_cpx_v = cp[0] * vec[0] + cp[1] * vec[1] +
                                              cp[2] * vec[2];

                            double sine = ps_cpx_v / len_edge;
                            beta = FEVV::Operators::Geometry::asin(sine);
                            found_face = true;
                            break;
                          }
                          else if(s11.find(source(*h1, g)) == s11.end())
                          {
                            s1.push(source(*h1, g));
                            s11.insert(source(*h1, g));
                          }
                        }
                      }
                      h1++;
                    } while(h1 != done1);
                    destination = get(pm, source(*h, g));
                  }
                }
              }
              else
              {
                // Case 4: 2 colinear faces
                if(normal2_nan &&
                   normal1_nan) // TODO BETTER: if(...) can be removed
                  beta = 0;
              }

              // ---
              // compute edge * edge^t * coeff, and add it to current matrix
              double p_vector_edge[3] = {vec[0], vec[1], vec[2]};
              double pp_matrix[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
              FEVV::Math::Matrix::Square::vector_times_transpose_mult(
                  p_vector_edge, pp_matrix, beta / len_edge);
              FEVV::Math::Matrix::Square::add(pp_matrix, pp_matrix_sum);
            }
          }
        }

        if(!isect)
        {
          VertexDescriptor w = source(
              *h, g); // previously : iterator here, NO !!! -> error from GL
                      // ??? Vertex_iterator w=h->opposite()->vertex();
          assert(w != GraphTraits::null_vertex());
          if(vertices.find(w) == vertices.end())
          {
            vertices.insert(w);
            s.push(w);
          }
        }
      }
    } while(++h != done);

    iter++;
  }

  double area = M_PI * radius * radius;

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      pp_matrix_sum[i][j] /= area;
}


//#define DBG_calculate_curvature_cmdm


/**
 * \brief  Calculate the curvature for a mesh
 *
 * \param  g                   input mesh
 * \param  v_cm                output vertex curvature map
 * \param  pm                  point map
 * \param  f_nm                face normal map
 * \param  isGeod              isGeod if true = geodesic neighborhood, if
 * false = 1-ring neighborhood \param  radius              radius \param
 * MinNrmMinCurvature  output minimum value of the minimum curvature (usefull
 * for color rendering) \param  MaxNrmMinCurvature  output maximum value of
 * the minimum curvature (usefull for color rendering) \param
 * MinNrmMaxCurvature  output minimum value of the maximum curvature (usefull
 * for color rendering) \param  MaxNrmMaxCurvature  output maximum value of
 * the maximum curvature (usefull for color rendering) \param  gt the geometry
 * traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 *
 * \ingroup  GenericManifoldFilters
 */
template< typename HalfedgeGraph,
          typename VertexCurvatureMap,
          typename PointMap,
          typename FaceNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
calculate_curvature_cmdm(const HalfedgeGraph &g,
                         VertexCurvatureMap &v_cm,
                         const PointMap &pm,
                         const FaceNormalMap &f_nm,
                         bool is_geod,
                         double radius,
                         double &min_nrm_min_curvature,
                         double &max_nrm_min_curvature,
                         double &min_nrm_max_curvature,
                         double &max_nrm_max_curvature,
                         const GeometryTraits &gt)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::vertex_descriptor Vertex_Descriptor;
  using Vector = typename GeometryTraits::Vector;

  min_nrm_min_curvature = 100000;
  max_nrm_min_curvature = -100000;

  min_nrm_max_curvature = 100000;
  max_nrm_max_curvature = -100000;

  // ---

  std::cout << "Asking to calculate curvature" << std::endl;
#ifdef DBG_calculate_curvature_cmdm
  std::cout << "Real radius = " << radius << std::endl;
#endif
  auto iterator_pair = vertices(g); // vertices() returns a vertex_iterator pair
  vertex_iterator vi = iterator_pair.first;
  vertex_iterator vi_end = iterator_pair.second;

  // #pragma omp parallel for
  for(vi; vi != vi_end; ++vi)
  {
    bool no_val_pro = false;

    double pp_matrix_sum[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    if(is_geod == true) // geodesic neighborhood
      geodes_principal_curvature_per_vert_cmdm<
          HalfedgeGraph,
          PointMap,
          FaceNormalMap,
          GeometryTraits,
          GraphTraits,
          typename GraphTraits::vertex_descriptor >(g, pm, f_nm, *vi,
                                                    pp_matrix_sum, radius, gt);

    else // 1-ring neighborhood
      principal_curvature_per_vert_cmdm<
          HalfedgeGraph,
          PointMap,
          FaceNormalMap,
          GeometryTraits,
          GraphTraits,
          typename GraphTraits::vertex_descriptor >(g, pm, f_nm, *vi,
                                                    pp_matrix_sum, gt);

    // Eigen values/vectors
    EigenMatrix cov_mat;
    EigenvectorsType vect_pro;
    EigenvalueType valpro;

    cov_mat(0, 0) = pp_matrix_sum[0][0];
    cov_mat(0, 1) = pp_matrix_sum[0][1];
    cov_mat(0, 2) = pp_matrix_sum[0][2];
    cov_mat(1, 0) = pp_matrix_sum[1][0];
    cov_mat(1, 1) = pp_matrix_sum[1][1];
    cov_mat(1, 2) = pp_matrix_sum[1][2];
    cov_mat(2, 0) = pp_matrix_sum[2][0];
    cov_mat(2, 1) = pp_matrix_sum[2][1];
    cov_mat(2, 2) = pp_matrix_sum[2][2];

#ifdef DBG_calculate_curvature_cmdm
    std::cout << "\n\nCovMat\n\n" << CovMat << std::endl;
#endif

    if(pp_matrix_sum[0][1] == 0 && pp_matrix_sum[0][2] == 0 &&
       pp_matrix_sum[1][0] == 0 && pp_matrix_sum[1][2] == 0 &&
       pp_matrix_sum[2][1] == 0 && pp_matrix_sum[2][0] == 0)
    {
      for(int i = 0; i < 3; i++)
        valpro(i) = cov_mat(i, i);
      vect_pro = EigenMatrix::Identity();
    }
    else
    {
      if(eigen_val_vect_compute(cov_mat, vect_pro, valpro) == -1)
      {
        no_val_pro = true;
      }
    }

#ifdef DBG_calculate_curvature_cmdm
    std::cout << "\n\nValpro\n\n" << Valpro << std::endl;
    std::cout << "\n\nVectPro\n\n" << VectPro << std::endl;
#endif

    if(!no_val_pro)
    {
      eigen_val_vect_sort(valpro, vect_pro);
      Vector v_kmax_curv(vect_pro(0, 1), vect_pro(1, 1), vect_pro(2, 1));
      Vector v_kmin_curv(vect_pro(0, 0), vect_pro(1, 0), vect_pro(2, 0));

#ifdef DBG_calculate_curvature_cmdm
      std::cout << "\n\nValpro sorted\n\n" << Valpro << std::endl;
      std::cout << "\n\nVectPro sorted\n\n" << VectPro << std::endl;
      std::cout << "\nVKmaxCurv = " << VKmaxCurv[0] << " " << VKmaxCurv[1]
                << " " << VKmaxCurv[2] << std::endl;
      std::cout << "\nVKminCurv = " << VKminCurv[0] << " " << VKminCurv[1]
                << " " << VKminCurv[2] << std::endl;
#endif

      v_cm[*vi].VKmaxCurv = v_kmax_curv;
      v_cm[*vi].VKminCurv = v_kmin_curv;

      if(valpro(0) > 0)
      {
        v_cm[*vi].KmaxCurv = valpro(0);
        v_cm[*vi].KminCurv = valpro(1);
      }
      else
      {
        v_cm[*vi].KmaxCurv = -valpro(0);
        v_cm[*vi].KminCurv = -valpro(1);
      }
    }
    else
    {
      v_cm[*vi].VKmaxCurv = gt.NULL_VECTOR;
      v_cm[*vi].VKminCurv = gt.NULL_VECTOR;
      v_cm[*vi].KmaxCurv = 0;
      v_cm[*vi].KminCurv = 0;
    }

#ifdef DBG_calculate_curvature_cmdm
    std::cout << "\nv_cm[*vi].KmaxCurv = " << v_cm[*vi].KmaxCurv << std::endl;
    std::cout << "\nv_cm[*vi].KminCurv = " << v_cm[*vi].KminCurv << std::endl;
#endif

    min_nrm_min_curvature =
        std::min< double >(min_nrm_min_curvature, v_cm[*vi].KminCurv);
    max_nrm_min_curvature =
        std::max< double >(max_nrm_min_curvature, v_cm[*vi].KminCurv);

    min_nrm_max_curvature =
        std::min< double >(min_nrm_max_curvature, v_cm[*vi].KmaxCurv);
    max_nrm_max_curvature =
        std::max< double >(max_nrm_max_curvature, v_cm[*vi].KmaxCurv);
  }

  // ---

  /*
  int cpt = 0;
  auto iterator_pair_log = vertices(g); // vertices() returns a
  vertex_iterator pair vertex_iterator vi_log = iterator_pair_log.first;
  vertex_iterator vend_log = iterator_pair_log.second; for (; vi_log !=
  vend_log; ++vi_log)
  {
          printf("cpt: %4d - Kmin: %.3lf - Kmax: %.3lf | VKmin: (%.3lf, %.3lf,
  %.3lf) | VKmax: (%.3lf, %.3lf, %.3lf)\n", cpt, v_cm[*vi_log].KminCurv,
  v_cm[*vi_log].KmaxCurv, v_cm[*vi_log].VKminCurv[0],
  v_cm[*vi_log].VKminCurv[1], v_cm[*vi_log].VKminCurv[2],
  v_cm[*vi_log].VKmaxCurv[0], v_cm[*vi_log].VKmaxCurv[1],
  v_cm[*vi_log].VKmaxCurv[2]); cpt++;
  }
  */
}

/**
 * \brief  Calculate the curvature for a mesh
 *
 * \param  g                   input mesh
 * \param  v_cm                output vertex curvature map
 * \param  pm                  point map
 * \param  f_nm                face normal map
 * \param  isGeod              isGeod if true = geodesic neighborhood, if
 * false = 1-ring neighborhood \param  radius              radius \param
 * MinNrmMinCurvature  output minimum value of the minimum curvature (usefull
 * for color rendering) \param  MaxNrmMinCurvature  output maximum value of
 * the minimum curvature (usefull for color rendering) \param
 * MinNrmMaxCurvature  output minimum value of the maximum curvature (usefull
 * for color rendering) \param  MaxNrmMaxCurvature  output maximum value of
 * the maximum curvature (usefull for color rendering)
 *
 * \sa     the variant that use the geometry traits provided by the user.
 *
 * \ingroup  GenericManifoldFilters
 */
template< typename HalfedgeGraph,
          typename VertexCurvatureMap,
          typename PointMap,
          typename FaceNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
calculate_curvature_cmdm(const HalfedgeGraph &g,
                         VertexCurvatureMap &v_cm,
                         const PointMap &pm,
                         const FaceNormalMap &f_nm,
                         bool is_geod,
                         double radius,
                         double &min_nrm_min_curvature,
                         double &max_nrm_min_curvature,
                         double &min_nrm_max_curvature,
                         double &max_nrm_max_curvature)
{
  GeometryTraits gt(g);
  calculate_curvature_cmdm< HalfedgeGraph,
                            VertexCurvatureMap,
                            PointMap,
                            FaceNormalMap,
                            GeometryTraits >(g,
                                             v_cm,
                                             pm,
                                             f_nm,
                                             is_geod,
                                             radius,
                                             min_nrm_min_curvature,
                                             max_nrm_min_curvature,
                                             min_nrm_max_curvature,
                                             max_nrm_max_curvature,
                                             gt);
}

} // namespace Filters
} // namespace FEVV
