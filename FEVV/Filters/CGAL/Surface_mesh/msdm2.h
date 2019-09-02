// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include "FEVV/DataStructures/DataStructures_cgal_surface_mesh.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

//---------------- specific code above --------------------

//---------------- generic code below ---------------------

// The code below is generic and works with all mesh datastructures
// without any change.

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#include "FEVV/Filters/Generic/Manifold/Curvature/curvature.hpp"
#include "FEVV/Filters/Generic/AABB/get_max_bb_size.hpp"

#include <boost/foreach.hpp>

namespace FEVV {

using Triangle = CGALKernel::Triangle_3;
using Segment = CGALKernel::Segment_3;
using Ray = CGALKernel::Ray_3;
using Vertex_Descriptor = MeshSurface::Vertex_index;

namespace Filters {

template< typename Point3d, typename Vector >
bool
sphere_clip_vector_msdm2(Point3d &o, double r, const Point3d &p, Vector &v)
{

  Vector w = p - o;
  double a = (v * v);
  double b = 2.0 * v * w;
  double c = (w * w) - r * r;
  double delta = b * b - 4 * a * c;

  if(a == 0)
    return true;

  if(delta < 0)
  {
    // Should not happen, but happens sometimes (numerical precision)

    return true;
  }
  double t = (-b + std::sqrt(delta)) / (2.0 * a);
  if(t < 0.0)
  {

    // Should not happen, but happens sometimes (numerical precision)
    return true;
  }
  if(t >= 1.0)
  {
    // Inside the sphere
    return false;
  }

  if(t < 1.e-16)
  {

    t = 0.01;
  }

  v = v * t;

  return true;
}

/**
\brief  Calculates the local curvature statistics per vertex.

\param  pVertex  The considered vertex
\param  TabDistance1  Curvature values from the neighborhood of pVertex
                      regarding the first polyhedron
\param  TabDistance2  Curvature values from the neighborhood of pVertex
                      regarding the second polyhedron
\param  TabPoint1     3D points from the neighborhoodof pVertex regarding the
                      first polyhedron
\param  TabPoint2     3D points from the neighborhood of pVertex regarding the
                      second polyhedron
\param  radius        The radius of the neighborhood.
\param  dim           The max dimension of the Bounding Box Length.


*/
template< typename VertexIterator,
          typename HalfedgeGraph,
          typename PointMap,
          typename MSDM2NearestMap,
          typename MSDM2Map >
void
compute_statistics(const HalfedgeGraph &m_degr,
                   const PointMap &pm_degr,
                   const MSDM2NearestMap &nearest_pmap,
                   MSDM2Map &msdm2_map,
                   const VertexIterator p_vertex,
                   double param,
                   std::vector< double > &tab_distance1,
                   std::vector< double > &tab_distance2,
                   std::vector< CGALPoint > &tab_point1,
                   std::vector< CGALPoint > &tab_point2,
                   double radius,
                   double dim)
{
  double moyenne1 = 0;
  double variance1 = 0;

  double moyenne2 = 0;
  double variance2 = 0;

  double covariance = 0;

  ////gaussian normalisation

  double variance = radius / 2;
  double *tab_wi1 = new double[tab_point1.size()];
  double *tab_wi2 = new double[tab_point1.size()];
  double somme_distance1 = 0;
  double somme_distance2 = 0;
  double somme_wi1 = 0;
  double somme_wi2 = 0;
  for(unsigned int i = 0; i < tab_point1.size(); i++) // MT
  {
    auto distance_pt1 = tab_point1[i] - get(pm_degr, p_vertex);
    double dist_pt1 = sqrt(distance_pt1 * distance_pt1);
    double wi1 = 1 / variance / sqrt(2 * 3.141592) *
                 exp(-(dist_pt1 * dist_pt1) / 2 / variance / variance);

    tab_wi1[i] = wi1;
    somme_wi1 += wi1;
    somme_distance1 += tab_distance1[i] * wi1;

    auto distance_pt2 = tab_point2[i] - get(nearest_pmap, p_vertex);
    double dist_pt2 = sqrt(distance_pt2 * distance_pt2);

    double wi2 = 1 / variance / sqrt(2 * 3.141592) *
                 exp(-(dist_pt2 * dist_pt2) / 2 / variance / variance);
    tab_wi2[i] = wi2;
    somme_wi2 += wi2;
    somme_distance2 += tab_distance2[i] * wi2;
  }
  moyenne1 = somme_distance1 / (double)somme_wi1;
  moyenne2 = somme_distance2 / (double)somme_wi2;


  for(unsigned int i = 0; i < tab_point1.size(); i++) // MT
  {
    variance1 += tab_wi1[i] * pow(tab_distance1[i] - moyenne1, 2);
    variance2 += tab_wi2[i] * pow(tab_distance2[i] - moyenne2, 2);
    covariance += tab_wi2[i] * (tab_distance1[i] - moyenne1) *
                  (tab_distance2[i] - moyenne2);
  }

  variance1 = variance1 / somme_wi1;
  variance2 = variance2 / somme_wi2;
  variance1 = sqrt(variance1);
  variance2 = sqrt(variance2);

  covariance = covariance / somme_wi1;

  /// we then compute the MSDM2_Local value
  // double C1=1; // MT
  // double C2=1; // MT
  // double C3=0.5; // MT
  double fact1 =
      (std::fabs(moyenne1 - moyenne2)) / (std::max(moyenne1, moyenne2) + 1);
  double fact2 =
      (std::fabs(variance1 - variance2)) / (std::max(variance1, variance2) + 1);
  double fact3 = (std::fabs(variance1 * variance2 - covariance)) /
                 (variance1 * variance2 + 1);


  auto tmp_msdm2 = get(msdm2_map, p_vertex);
  tmp_msdm2 += std::pow(
      (std::pow(fact1, 1) + std::pow(fact2, 1) + param * std::pow(fact3, 1)) /
          (2. + param),
      1. / 1.);
  put(msdm2_map, p_vertex, tmp_msdm2);


  delete[] tab_wi1;
  delete[] tab_wi2;
}


/**
\brief	Initialize the matching process

\param	m_PolyDegrad	  	The first Polyhedron.
\param	m_PolyOriginal	  	The second Polyhedron.
\param	[out] _TabMatchedFacet	Facets from m_PolyOriginal on which vertices
from m_PolyDegrad are projected


*/
template< typename MSDM2Map,
          typename HalfedgeGraph,
          typename MSDM2NearestMap,
          typename PointMap,
          typename TagMap,
          typename FaceIterator =
              typename boost::graph_traits< HalfedgeGraph >::face_iterator >
void
matching_multires_init(const HalfedgeGraph &m_poly_degrad,
                       MSDM2Map &msdm2m_degrad,
                       MSDM2NearestMap &msdm2nm_degrad,
                       TagMap &tag_map,
                       const PointMap &pm_degrad,
                       const HalfedgeGraph &m_poly_original,
                       CGAL::SM_Face_index *tab_matched_facet)
{
  using AABB_Primitive =
      typename CGAL::AABB_face_graph_triangle_primitive< HalfedgeGraph >;
  using AABB_Traits = typename CGAL::AABB_traits< CGALKernel, AABB_Primitive >;
  using AABB_Tree = typename CGAL::AABB_tree< AABB_Traits >;
  using Object_and_primitive_id = typename AABB_Tree::Object_and_primitive_id;
  using Point_and_primitive_id = typename AABB_Tree::Point_and_primitive_id;


  auto original_faces = faces(m_poly_original);

  FaceIterator orig_begin = original_faces.first;
  FaceIterator orig_end = original_faces.second;

  AABB_Tree tree(orig_begin, orig_end, m_poly_original);
  tree.accelerate_distance_queries();

  // Searching for the closest point and facet for each vertex

  int ind = 0;
  BOOST_FOREACH(Vertex_Descriptor vi, vertices(m_poly_degrad))
  {
    put(msdm2m_degrad, vi, 0.);
    put(tag_map, vi, ind);
    // computes closest point and primitive id
    Point_and_primitive_id pp =
        tree.closest_point_and_primitive(get(pm_degrad, vi));
    auto f_nearest = pp.second; // closest primitive id

    put(msdm2nm_degrad, vi, pp.first);
    tab_matched_facet[ind] = f_nearest;

    ind++;
  }
}


template< typename PointMap,
          typename HalfedgeGraph,
          typename FaceNormalMap,
          typename VertexMSDM2Map,
          typename GeometryTraits >
double
process_msdm2_multires(const HalfedgeGraph &m_poly_degrad,
                       const PointMap &pm_degrad,
                       const FaceNormalMap &fnm_degrad,
                       FEVV::PMapsContainer &pmaps_degrad,
                       const GeometryTraits &gt_degrad,
                       const HalfedgeGraph &m_poly_original,
                       const PointMap &pm_original,
                       const FaceNormalMap &fnm_original,
                       FEVV::PMapsContainer &pmaps_original,
                       const GeometryTraits &gt_original,
                       const int nb_level,
                       VertexMSDM2Map &msdm2_pmap,
                       const double maxdim,
                       double &msdm2_value)
{
  std::vector< CGALPoint > tab_point1;
  std::vector< CGALPoint > tab_point2;

  double somme_msdm2 = 0;
  int nb_vertex = 0;

  auto curvature_map_degr =
      FEVV::make_vertex_property_map< HalfedgeGraph,
                                      FEVV::Filters::v_Curv< HalfedgeGraph > >(
          m_poly_degrad);
  auto curvature_map_orig =
      FEVV::make_vertex_property_map< HalfedgeGraph,
                                      FEVV::Filters::v_Curv< HalfedgeGraph > >(
          m_poly_original);

  CGAL::SM_Face_index *tab_matched_facet =
      new CGAL::SM_Face_index[size_of_vertices(m_poly_degrad)];


  VertexMSDM2Map msdm2_match_pmap;
  if(has_map(pmaps_degrad, std::string("v:msdm2_match")))
  {
    msdm2_match_pmap =
        boost::any_cast< VertexMSDM2Map >(pmaps_degrad.at("v:msdm2_match"));
  }
  else
  {
    msdm2_match_pmap =
        FEVV::make_vertex_property_map< HalfedgeGraph, double >(m_poly_degrad);
    pmaps_degrad["v:msdm2_match"] = msdm2_match_pmap;
  }

  typedef typename FEVV::Vertex_pmap< HalfedgeGraph, CGALPoint >
      vertex_msdm2_nearest_map;
  vertex_msdm2_nearest_map msdm2_nearest_pmap;

  if(has_map(pmaps_degrad, std::string("v:msdm2_nearest")))
  {
    msdm2_nearest_pmap = boost::any_cast< vertex_msdm2_nearest_map >(
        pmaps_degrad.at("v:msdm2_nearest"));
  }
  else
  {
    msdm2_nearest_pmap =
        FEVV::make_vertex_property_map< HalfedgeGraph, CGALPoint >(
            m_poly_degrad);
    pmaps_degrad["v:msdm2_nearest"] = msdm2_nearest_pmap;
  }

  typedef typename FEVV::Vertex_pmap< HalfedgeGraph, int > vertex_tag_map;
  vertex_tag_map tag_map;

  if(has_map(pmaps_degrad, std::string("v:tag")))
  {
    tag_map = boost::any_cast< vertex_tag_map >(pmaps_degrad.at("v:tag"));
  }
  else
  {
    tag_map =
        FEVV::make_vertex_property_map< HalfedgeGraph, int >(m_poly_degrad);
    pmaps_degrad["v:tag"] = tag_map;
  }

  matching_multires_init(m_poly_degrad,
                         msdm2_pmap,
                         msdm2_nearest_pmap,
                         tag_map,
                         pm_degrad,
                         m_poly_original,
                         tab_matched_facet);

  double radius_curvature = 0.002;
  for(int i = 0; i < nb_level; i++)
  {
    double min_nrm_min_curvature_deg, max_nrm_min_curvature_deg,
        min_nrm_max_curvature_deg, max_nrm_max_curvature_deg;
    double min_nrm_min_curvature_org, max_nrm_min_curvature_org,
        min_nrm_max_curvature_org, max_nrm_max_curvature_org;

    FEVV::Filters::calculate_curvature(m_poly_degrad,
                                       curvature_map_degr,
                                       pm_degrad,
                                       fnm_degrad,
                                       true,
                                       maxdim * radius_curvature,
                                       min_nrm_min_curvature_deg,
                                       max_nrm_min_curvature_deg,
                                       min_nrm_max_curvature_deg,
                                       max_nrm_max_curvature_deg);
    FEVV::Filters::calculate_curvature(m_poly_original,
                                       curvature_map_orig,
                                       pm_original,
                                       fnm_original,
                                       true,
                                       maxdim * radius_curvature,
                                       min_nrm_min_curvature_org,
                                       max_nrm_min_curvature_org,
                                       min_nrm_max_curvature_org,
                                       max_nrm_max_curvature_org);

    kmax_kmean(m_poly_original, curvature_map_orig, maxdim);
    kmax_kmean(m_poly_degrad, curvature_map_degr, maxdim);

    matching_multires_update(m_poly_degrad,
                             pm_degrad,
                             m_poly_original,
                             pm_original,
                             curvature_map_orig,
                             msdm2_match_pmap,
                             tab_matched_facet);

    BOOST_FOREACH(Vertex_Descriptor vi, vertices(m_poly_degrad))
    {
      std::vector< double > tab_distance1;
      std::vector< double > tab_distance2;

      tab_point1.clear();
      tab_point2.clear();

      process_msdm2_per_vertex(m_poly_degrad,
                               pm_degrad,
                               tag_map,
                               curvature_map_degr,
                               msdm2_nearest_pmap,
                               msdm2_match_pmap,
                               gt_degrad,
                               m_poly_original,
                               pm_original,
                               vi,
                               radius_curvature * 5 * maxdim,
                               tab_distance1,
                               tab_distance2,
                               tab_point1,
                               tab_point2);
      compute_statistics(m_poly_degrad,
                         pm_degrad,
                         msdm2_nearest_pmap,
                         msdm2_pmap,
                         vi,
                         0.5,
                         tab_distance1,
                         tab_distance2,
                         tab_point1,
                         tab_point2,
                         radius_curvature * 5 * maxdim,
                         maxdim);
    }


    radius_curvature += 0.001;
  }

  BOOST_FOREACH(Vertex_Descriptor vi, vertices(m_poly_degrad))
  {
    double l_msdm2 = get(msdm2_pmap, vi) / nb_level;
    put(msdm2_pmap, vi, l_msdm2);
    somme_msdm2 += std::pow(l_msdm2, 3);
    nb_vertex++;
  }


  msdm2_value = std::pow(somme_msdm2 / (double)nb_vertex, 1. / 3.);

  delete[] tab_matched_facet;

  return 0;
}

/**
\brief Updates the matching process

\param	m_PolyDegrad	  	The first Polyhedron.
\param	_TabMatchedFacet	Facets from m_PolyOriginal on which vertices
from m_PolyDegrad are projected


*/
template< typename PointMap,
          typename HalfedgeGraph,
          typename CurvatureVertexMap,
          typename MSMD2CurvatureMatchMap >
void
matching_multires_update(const HalfedgeGraph &m_poly_degrad,
                         const PointMap &pm_degr,
                         const HalfedgeGraph &m_poly_orig,
                         const PointMap &pm_orig,
                         const CurvatureVertexMap curv_orig,
                         MSMD2CurvatureMatchMap curvmatch_pmap,
                         CGAL::SM_Face_index *tab_matched_facet)
{

  int ind = 0;
  BOOST_FOREACH(Vertex_Descriptor vi, vertices(m_poly_degrad))
  {
    // Point3d Nearest=pVertex->match; // MT
    auto *f_nearest = &tab_matched_facet[ind];

    ind++;
    /// calculation of the nearest point curvature value using vertices of the
    /// Nearest triangle
    // we use linear interpolation using barycentric coordinates
    auto x1 =
        get(pm_orig, target(halfedge(*f_nearest, m_poly_orig), m_poly_orig));
    auto x2 = get(pm_orig,
                  target(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                         m_poly_orig));
    auto x3 =
        get(pm_orig,
            target(next(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                        m_poly_orig),
                   m_poly_orig));

    double l1 = sqrt((x3 - x2) * (x3 - x2));
    double l2 = sqrt((x1 - x3) * (x1 - x3));
    double l3 = sqrt((x1 - x2) * (x1 - x2));

    auto v1 = x1 - get(pm_degr, vi);
    auto v2 = x2 - get(pm_degr, vi);
    auto v3 = x3 - get(pm_degr, vi);

    double t1 = sqrt(v1 * v1);
    double t2 = sqrt(v2 * v2);
    double t3 = sqrt(v3 * v3);

    double p1 = (l1 + t2 + t3) / 2;
    double p2 = (t1 + l2 + t3) / 2;
    double p3 = (t1 + t2 + l3) / 2;

    double a1 = (p1 * (p1 - l1) * (p1 - t3) * (p1 - t2));
    double a2 = (p2 * (p2 - l2) * (p2 - t3) * (p2 - t1));
    double a3 = (p3 * (p3 - l3) * (p3 - t1) * (p3 - t2));

    if(a1 > 0)
      a1 = sqrt(a1);
    else
      a1 = 0;
    if(a2 > 0)
      a2 = sqrt(a2);
    else
      a2 = 0;
    if(a3 > 0)
      a3 = sqrt(a3);
    else
      a3 = 0;

    double c1 =
        get(curv_orig, target(halfedge(*f_nearest, m_poly_orig), m_poly_orig))
            .KmaxCurv;
    double c2 = get(curv_orig,
                    target(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                           m_poly_orig))
                    .KmaxCurv;
    double c3 =
        get(curv_orig,
            target(next(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                        m_poly_orig),
                   m_poly_orig))
            .KmaxCurv;

    double curvmatch = 0.0;
    if((a1 + a2 + a3) > 0)
      curvmatch = (a1 * c1 + a2 * c2 + a3 * c3) / (a1 + a2 + a3);
    else
      curvmatch = (c1 + c2 + c3) / 3;
    put(curvmatch_pmap, vi, curvmatch);
  }
}

/**
\brief  Computes the local neighborhoods


\param  pVertex  The considered vertex
\param  radius   radius of the neighborhood
\param [out] TabDistance1  Curvature values from the neighborhood of pVertex
                           regarding the first polyhedron
\param [out] TabDistance2  Curvature values from the neighborhood of pVertex
                           regarding the second polyhedron
\param [out] TabPoint1     3D points from the neighborhoodof pVertex regarding
                           the first polyhedron
\param [out] TabPoint2     3D points from the neighborhood of pVertex regarding
                           the second polyhedron
*/
template< typename VertexIterator,
          typename HalfedgeGraph,
          typename PointMap,
          typename TagMap,
          typename CurvatureVertexMap,
          typename MSDM2NearestMap,
          typename MSMD2CurvatureMatchMap,
          typename GeometryTraits >
void
process_msdm2_per_vertex(const HalfedgeGraph &m_degr,
                         const PointMap &pm_degr,
                         const TagMap &tags_pmap,
                         const CurvatureVertexMap &curv_pmap,
                         const MSDM2NearestMap &nearest_pmap,
                         const MSMD2CurvatureMatchMap &curvmatch_pmap,
                         const GeometryTraits &gt_degr,
                         const HalfedgeGraph &m_orig,
                         const PointMap &pm_orig,
                         const VertexIterator &p_vertex,
                         double radius,
                         std::vector< double > &tab_distance1,
                         std::vector< double > &tab_distance2,
                         std::vector< CGALPoint > &tab_point1,
                         std::vector< CGALPoint > &tab_point2)
{

  std::set< int > vertices;

  std::stack< VertexIterator > s;

  auto o = get(pm_degr, p_vertex);

  s.push(p_vertex);
  vertices.insert(get(tags_pmap, p_vertex));


  tab_distance1.push_back(get(curv_pmap, p_vertex).KmaxCurv);
  tab_point1.push_back(get(pm_degr, p_vertex));

  tab_distance2.push_back(get(curvmatch_pmap, p_vertex));
  tab_point2.push_back(get(nearest_pmap, p_vertex));

  int nb_sommet_in_sphere = 0;


  while(!s.empty())
  {
    VertexIterator v_it = s.top();
    s.pop();
    CGALPoint p = get(pm_degr, v_it);
    auto h = halfedge(v_it, m_degr);
    auto h0 = h;
    do
    {
      auto p1 = get(pm_degr, target(h, m_degr));
      auto p2 = get(pm_degr, target(opposite(h, m_degr), m_degr));

      auto p1m = get(nearest_pmap, target(h, m_degr));
      auto p2m = get(nearest_pmap, target(opposite(h, m_degr), m_degr));

      auto v = (p2 - p1);
      auto vm = (p2m - p1m);

      if(v_it == p_vertex || v * (p - o) > 0.0)
      {

        double len_old = std::sqrt(v * v);
        bool isect = sphere_clip_vector_msdm2(o, radius, p, v);
        double len_edge = std::sqrt(v * v);

        nb_sommet_in_sphere++;


        double weighted_curv1, weighted_curv2;
        CGALPoint weighted_p1, weighted_p2;

        bool is_already_integrated = false;
        if(!isect)
        {

          auto w = target(opposite(h, m_degr), m_degr);
          if(vertices.find(get(tags_pmap, w)) == vertices.end())
          {
            vertices.insert(get(tags_pmap, w));
            s.push(w);
          }
          else
            is_already_integrated = true;
        }

        if(is_already_integrated == false)
        {
          if(len_old != 0 && isect)
          {
            weighted_curv1 =
                (1 - len_edge / len_old) *
                    get(curv_pmap, target(h, m_degr)).KmaxCurv +
                len_edge / len_old *
                    get(curv_pmap, target(opposite(h, m_degr), m_degr))
                        .KmaxCurv;
            weighted_p1 = p1 + v;

            weighted_curv2 =
                (1 - len_edge / len_old) *
                    get(curvmatch_pmap, target(h, m_degr)) +
                len_edge / len_old *
                    get(curvmatch_pmap, target(opposite(h, m_degr), m_degr));
            weighted_p2 = p1m + (len_edge / len_old) * vm;
          }
          else
          {
            weighted_curv1 =
                get(curv_pmap, target(opposite(h, m_degr), m_degr)).KmaxCurv;
            weighted_curv2 =
                get(curvmatch_pmap, target(opposite(h, m_degr), m_degr));
            weighted_p1 = p2;
            weighted_p2 = p2m;
          }

          tab_distance1.push_back(weighted_curv1);
          tab_point1.push_back(weighted_p1);

          tab_distance2.push_back(weighted_curv2);
          tab_point2.push_back(weighted_p2);
        }
      }
      h = opposite(next(h, m_degr), m_degr);
    } while(h != h0);
  }
}

/**
\brief	Computes the mean curvature field (kmin+kmax)/2 and normalize it
according to the size of the model

\param	polyhedron_ptr	The polyhedron.
\param	coef		  	The normalization coef.
*/
template< typename HalfedgeGraph, typename CurvatureVertexMap >
void
kmax_kmean(const HalfedgeGraph &mesh,
           CurvatureVertexMap &curv_vmap,
           const double coef)
{

  BOOST_FOREACH(Vertex_Descriptor vi, vertices(mesh))
  {
    auto curv = get(curv_vmap, vi);
    double kmax = curv.KmaxCurv * coef;
    double kmin = std::abs(curv.KminCurv * coef);
    curv.KmaxCurv = (kmax + kmin) / 2.;
    put(curv_vmap, vi, curv);
  }
}

/**
\brief	Computes the multiscale MSDM2 metric

\param	m_PolyDegrad	  	The first Polyhedron.
\param	m_PolyOriginal	  	The second Polyhedron.
\param	NbLevel		Number of scales used
\param	maxdim	The max dimension of the Bounding Box Length.
\param [out] MSDM2Value	The computed value of the MSDM2 metric
*/
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >,
          typename PointMap,
          typename FaceNormalMap,
          typename VertexMSDM2Map >
double
process_msdm2_multires(const HalfedgeGraph &m_poly_degrad,
                       const PointMap &pm_degrad,
                       const FaceNormalMap &fnm_degrad,
                       FEVV::PMapsContainer &pmaps_degrad,
                       const HalfedgeGraph &m_poly_original,
                       const PointMap &pm_original,
                       const FaceNormalMap &fnm_original,
                       FEVV::PMapsContainer &pmaps_original,
                       const int nb_level,
                       VertexMSDM2Map &msdm2_pmap,
                       double &msdm2_value)
{

  GeometryTraits gt_deg(m_poly_degrad);
  GeometryTraits gt_org(m_poly_original);

  return process_msdm2_multires(m_poly_degrad,
                                pm_degrad,
                                fnm_degrad,
                                pmaps_degrad,
                                gt_deg,
                                m_poly_original,
                                pm_original,
                                fnm_original,
                                pmaps_original,
                                gt_org,
                                nb_level,
                                msdm2_pmap,
                                Filters::get_max_bb_size(m_poly_degrad,
                                                    pm_degrad,
                                                    gt_deg),
                                msdm2_value);
}

} // namespace Filters
} // namespace FEVV
