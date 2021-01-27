// Copyright (c) 2012-2020 University of Lyon and CNRS (France).
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

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#include "curvature_cmdm.hpp"
#include <FEVV/Filters/Generic/AABB/get_max_bb_size.hpp>
#include <FEVV/Filters/CGAL/Surface_mesh/CMDM/rgb_to_lab.h>
#include <FEVV/Filters/CGAL/Surface_mesh/CMDM/transfo_lab2lab2000hl.h>

#include <boost/foreach.hpp>

#ifdef _MSC_VER
// disable some warnings on Windows
#pragma warning(push)
#pragma warning(disable : 4244)
#pragma warning(disable : 4267)
// 4244 & 4267: converting type A to type B, possible data loss
#endif

namespace FEVV {

using Triangle = CGALKernel::Triangle_3;
using Segment = CGALKernel::Segment_3;
using Ray = CGALKernel::Ray_3;
using Vertex_Descriptor = MeshSurface::Vertex_index;

namespace Filters {

template< typename Point3d, typename Vector >
bool
sphere_clip_vector_cmdm(Point3d &o, double r, const Point3d &p, Vector &v)
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
\brief  Calculates the local geometry-based features per vertex.

\param  pVertex       The considered vertex
\param  TabDistance1  Curvature values from the neighborhood of pVertex
regarding the first polyhedron
\param  TabDistance2  Curvature values from the neighborhood of pVertex
regarding the second polyhedron
\param  TabPoint1     3D points from the neighborhood of pVertex regarding the
first polyhedron
\param  TabPoint2     3D points from the neighborhood of pVertex regarding the
second polyhedron
\param  radius        The radius of the neighborhood.
*/

template< typename VertexIterator,
          typename HalfedgeGraph,
          typename PointMap,
          typename CMDMNearestMap >
void
compute_geometry_statistics(const HalfedgeGraph &m_degr,
                            const PointMap &pm_degr,
                            const CMDMNearestMap &nearest_pmap,
                            const VertexIterator p_vertex,
                            std::vector< double > &tab_distance1,
                            std::vector< double > &tab_distance2,
                            std::vector< CGALPoint > &tab_point1,
                            std::vector< CGALPoint > &tab_point2,
                            double *tab_f1,
                            double *tab_f2,
                            double *tab_f3,
                            double radius)
{
  double moyenne1 = 0;
  double variance1 = 0;

  double moyenne2 = 0;
  double variance2 = 0;

  double covariance = 0;
  double variance = radius / 3;
  std::vector< double > tab_wi1;
  std::vector< double > tab_wi2;
  double somme_distance1 = 0;
  double somme_distance2 = 0;
  double somme_wi1 = 0;
  double somme_wi2 = 0;
  for(unsigned int i = 0; i < tab_point1.size(); i++)
  {
    auto distance_pt1 = tab_point1[i] - get(pm_degr, p_vertex);
    double dist_pt1 = sqrt(distance_pt1 * distance_pt1);
    double wi1 = 1 / variance / sqrt(2 * 3.141592) *
                 exp(-(dist_pt1 * dist_pt1) / 2 / variance / variance);

    tab_wi1.push_back(wi1);
    somme_wi1 += wi1;
    somme_distance1 += tab_distance1[i] * wi1;

    auto distance_pt2 = tab_point2[i] - get(nearest_pmap, p_vertex);
    double dist_pt2 = sqrt(distance_pt2 * distance_pt2);

    double wi2 = 1 / variance / sqrt(2 * 3.141592) *
                 exp(-(dist_pt2 * dist_pt2) / 2 / variance / variance);
    tab_wi2.push_back(wi2);
    somme_wi2 += wi2;
    somme_distance2 += tab_distance2[i] * wi2;
  }
  moyenne1 = somme_distance1 / (double)somme_wi1;
  moyenne2 = somme_distance2 / (double)somme_wi2;


  for(unsigned int i = 0; i < tab_point1.size(); i++)
  {
    variance1 += tab_wi1[i] * pow(tab_distance1[i] - moyenne1, 2);
    variance2 += tab_wi2[i] * pow(tab_distance2[i] - moyenne2, 2);

    covariance += tab_wi1[i] * (tab_distance1[i] - moyenne1) *
                  (tab_distance2[i] - moyenne2);
  }

  variance1 = variance1 / somme_wi1;
  variance2 = variance2 / somme_wi2;
  variance1 = sqrt(variance1);
  variance2 = sqrt(variance2);
  covariance = covariance / somme_wi1;
  covariance = (covariance < 0) ? 0 : covariance;

  // We then compute the local geometry-based features
  // 1) Curvature comparison
  double f1 =
      (std::fabs(moyenne1 - moyenne2)) / (std::max(moyenne1, moyenne2) + 1);
  // 2) Curvature contrast
  double f2 =
      (std::fabs(variance1 - variance2)) / (std::max(variance1, variance2) + 1);
  // 3) Curvature structure
  double f3 = (std::fabs(variance1 * variance2 - covariance)) /
              (variance1 * variance2 + 1);

  tab_f1[p_vertex] += f1;
  tab_f2[p_vertex] += f2;
  tab_f3[p_vertex] += f3;
}

/**
\brief  Calculates the local color-based features per vertex.

\param  pVertex    The considered vertex
\param  TabColor1  Color values (L,A,B) from the neighborhood of pVertex
regarding the first polyhedron
\param  TabColor2  Color values (L,A,B) from the neighborhood of pVertex
regarding the second polyhedron
\param  TabPoint1  3D points from the neighborhood of pVertex regarding the
first polyhedron
\param  TabPoint2  3D points from the neighborhood of pVertex regarding the
second polyhedron
\param  radius     The radius of the neighborhood.
*/
template< typename VertexIterator,
          typename HalfedgeGraph,
          typename PointMap,
          typename CMDMNearestMap >
void
compute_color_statistics(const HalfedgeGraph &m_degr,
                         const PointMap &pm_degr,
                         const CMDMNearestMap &nearest_pmap,
                         const VertexIterator p_vertex,
                         std::vector< std::array< double, 3 > > &tab_color1,
                         std::vector< std::array< double, 3 > > &tab_color2,
                         std::vector< CGALPoint > &tab_point1,
                         std::vector< CGALPoint > &tab_point2,
                         double *tab_f4,
                         double *tab_f5,
                         double *tab_f6,
                         double *tab_f7,
                         double *tab_f8,
                         double radius)
{
  double variance = radius / 3;
  std::vector< double > tab_wi1;
  std::vector< double > tab_wi2;
  double somme_wi1 = 0;
  double somme_wi2 = 0;

  for(unsigned int i = 0; i < tab_point1.size(); i++)
  {
    auto distance_pt1 = tab_point1[i] - get(pm_degr, p_vertex);
    double dist_pt1 = sqrt(distance_pt1 * distance_pt1);
    double wi1 = 1 / variance / sqrt(2 * 3.141592) *
                 exp(-(dist_pt1 * dist_pt1) / 2 / variance / variance);
    tab_wi1.push_back(wi1);
    somme_wi1 += wi1;

    auto distance_pt2 = tab_point2[i] - get(nearest_pmap, p_vertex);
    double dist_pt2 = sqrt(distance_pt2 * distance_pt2);
    double wi2 = 1 / variance / sqrt(2 * 3.141592) *
                 exp(-(dist_pt2 * dist_pt2) / 2 / variance / variance);
    tab_wi2.push_back(wi2);
    somme_wi2 += wi2;
  }

  std::vector< double > tab_Chr1; // Chroma from the neighborhood of pVertex
                                  // regarding the first polyhedron
  std::vector< double > tab_Chr2; // Chroma from the neighborhood of pVertex
                                  // regarding the second polyhedron
  std::vector< double > tab_dH;   // Hue difference

  double muL1 = 0, muCh1 = 0;
  double variance1_L = 0, sL1 = 0;

  double muL2 = 0, muCh2 = 0;
  double variance2_L = 0, sL2 = 0;

  double dL_sq = 0, dCh_sq = 0, dH_sq = 0, sL12 = 0; // Mixed terms

  double somme1_L = 0, somme1_Chr = 0;
  double somme2_L = 0, somme2_Chr = 0;
  for(unsigned int i = 0; i < tab_point1.size(); i++)
  {
    // Chroma  = sqrt(A^2 +  B^2)
    tab_Chr1.push_back(sqrt(tab_color1[i][1] * tab_color1[i][1] +
                            tab_color1[i][2] * tab_color1[i][2]));
    somme1_L += tab_color1[i][0] * tab_wi1[i];
    somme1_Chr += tab_Chr1[i] * tab_wi1[i];

    tab_Chr2.push_back(sqrt(tab_color2[i][1] * tab_color2[i][1] +
                            tab_color2[i][2] * tab_color2[i][2]));
    somme2_L += tab_color2[i][0] * tab_wi2[i];
    somme2_Chr += tab_Chr2[i] * tab_wi2[i];

    // Hue difference = sqrt[(A1-A2)^2 + (B1-B2)^2 - (Chr1-Chr2)^2]
    double dH = 0;
    dH = std::pow(tab_color1[i][1] - tab_color2[i][1], 2) +
         std::pow(tab_color1[i][2] - tab_color2[i][2], 2) -
         std::pow(tab_Chr1[i] - tab_Chr2[i], 2);
    tab_dH.push_back((dH < 0) ? 0 : sqrt(dH));
  }

  muL1 = somme1_L / somme_wi1;
  muCh1 = somme1_Chr / somme_wi1;

  muL2 = somme2_L / somme_wi2;
  muCh2 = somme2_Chr / somme_wi2;

  // Compute standard deviation and mixed terms
  for(unsigned int i = 0; i < tab_point1.size(); i++)
  {
    variance1_L += tab_wi1[i] * pow(tab_color1[i][0] - muL1, 2);
    variance2_L += tab_wi2[i] * pow(tab_color2[i][0] - muL2, 2);

    dH_sq += tab_wi1[i] * tab_dH[i];
    sL12 += tab_wi1[i] * (tab_color1[i][0] - muL1) *
            (tab_color2[i][0] - muL2); // covariance
  }

  variance1_L = variance1_L / somme_wi1;
  sL1 = (variance1_L < 0) ? 0 : sqrt(variance1_L);

  variance2_L = variance2_L / somme_wi2;
  sL2 = (variance2_L < 0) ? 0 : sqrt(variance2_L);

  sL12 = sL12 / somme_wi1;
  sL12 = (sL12 < 0) ? 0 : sL12;

  dL_sq = (muL1 - muL2) * (muL1 - muL2);
  dCh_sq = (muCh1 - muCh2) * (muCh1 - muCh2);

  dH_sq = dH_sq / somme_wi1;
  dH_sq = dH_sq * dH_sq;

  // We then compute the local color-based features
  double c1 = 0.002, c2 = 0.1, c3 = 0.1, c4 = 0.002,
         c5 = 0.008; // The constants
                     // 1) Lightness comparison
  double f4 = 1 - (1 / (c1 * dL_sq + 1));
  // 2) Lightness-contrast comparison
  double f5 = 1 - ((2 * sL1 * sL2 + c2) / (sL1 * sL1 + sL2 * sL2 + c2));
  // 3) Lightness-structure comparison
  double f6 = 1 - ((sL12 + c3) / (sL1 * sL2 + c3));
  // 4) Chroma comparison
  double f7 = 1 - (1 / (c4 * dCh_sq + 1));
  // 5) Hue comparison
  double f8 = 1 - (1 / (c5 * dH_sq + 1));

  tab_f4[p_vertex] += f4;
  tab_f5[p_vertex] += f5;
  tab_f6[p_vertex] += f6;
  tab_f7[p_vertex] += f7;
  tab_f8[p_vertex] += f8;
}

/**
\brief  Retrieve the color (RGB values) of each vertex of the mesh
and transform it into the working color space LAB2000HL.

\param	m_Poly	The Polyhedron.
\param	init_grid_L, grid_L, init_grid_AB Grids used for the color
transformation LAB -> LAB2000HL. \Return	[out] _LAB2000HL_Vertexcolormap
The vertex color map (in LAB2000HL space) of the mesh.
*/
template< typename HalfedgeGraph,
          typename VertexColorMap = typename FEVV::
              PMap_traits< FEVV::vertex_color_t, MeshSurface >::pmap_type >
VertexColorMap
Get_VertexColorMap(
    const HalfedgeGraph &m_poly,
    FEVV::PMapsContainer &pmaps,
    std::vector< double > &init_grid_L,
    std::vector< double > &grid_L,
    std::vector< std::vector< std::pair< double, double > > > &init_grid_AB)
{
  VertexColorMap v_cm;
  v_cm = get_property_map(FEVV::vertex_color, m_poly, pmaps);
  
  // Transform the color of each vertex into the LAB2000HL space
  // RGB -> LAB -> LAB2000HL
  using Color = typename boost::property_traits< VertexColorMap >::value_type;                       
  VertexColorMap _LAB2000HL_Vertexcolormap;

  BOOST_FOREACH(Vertex_Descriptor vi, vertices(m_poly))
  {
	double *RGB_Vcolor = new double[3];  // vertex color in RGB space
	double *LAB_Vcolor = new double[3];       // vertex color in LAB space
    double *LAB2000HL_Vcolor = new double[3]; // vertex color in LAB2000HL space
	  
	RGB_Vcolor[0] = (v_cm)[vi][0];
	RGB_Vcolor[1] = (v_cm)[vi][1];
	RGB_Vcolor[2] = (v_cm)[vi][2];

    // Color space conversion.
    rgb_to_lab(RGB_Vcolor[0], RGB_Vcolor[1], RGB_Vcolor[2], LAB_Vcolor);
    lab2000hl_conversion(LAB_Vcolor[0],
                         LAB_Vcolor[1],
                         LAB_Vcolor[2],
                         LAB2000HL_Vcolor,
                         init_grid_L,
                         grid_L,
                         init_grid_AB);

    put(_LAB2000HL_Vertexcolormap,
        vi,
        Color(LAB2000HL_Vcolor[0], LAB2000HL_Vcolor[1], LAB2000HL_Vcolor[2])); 

	delete[] RGB_Vcolor;
	delete[] LAB_Vcolor;
	delete[] LAB2000HL_Vcolor;
  }

  return _LAB2000HL_Vertexcolormap;
}

/**
\brief	Initialize the matching process

\param	m_PolyDegrad	  	The first Polyhedron.
\param	m_PolyOriginal	  	The second Polyhedron.
\param	[out] _TabMatchedFacet	Facets from m_PolyOriginal on which vertices
from m_PolyDegrad are projected
*/

template< typename CMDMMap,
          typename HalfedgeGraph,
          typename CMDMNearestMap,
          typename TagMap,
          typename PointMap,
          typename FaceIterator =
              typename boost::graph_traits< HalfedgeGraph >::face_iterator >
void
matching_multires_init_cmdm(
    const HalfedgeGraph &m_poly_degrad,
    CMDMMap &cmdm_degrad,
    CMDMNearestMap &cmdmm_degrad,
    TagMap &tag_map,
    const PointMap &pm_degrad, // Or CGALPoint *tab_pm_degrad,
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

  // CGALPoint *tab_cmdmm_degrad = new CGALPoint[(int)num_vertices(m_poly_degrad)]; - used in parallelization trials 

  //#pragma omp parallel for - program crashes with a critical error
							// detected c0000374 (problem with Parallelism of
							// AABB_tree::closest_point_and_primitive() function)
  for(int i = 0; i < (int)num_vertices(m_poly_degrad); i++) // BOOST_FOREACH(Vertex_Descriptor vi, vertices(m_poly_degrad))
  {
    Vertex_Descriptor vi(i);
    // computes closest point and primitive id
    Point_and_primitive_id pp =
        tree.closest_point_and_primitive(get(pm_degrad, vi));
    auto f_nearest = pp.second; // closest primitive id

    /* formulation used in parallelization trials
    auto *p = &tab_pm_degrad[vi];
    Point_and_primitive_id pp =(tree.closest_point_and_primitive(*p));
    tab_cmdmm_degrad[i] = pp.first;
    */

    put(cmdmm_degrad, vi, pp.first);
    tab_matched_facet[i] = f_nearest;
  }

  int ind = 0;
  BOOST_FOREACH(Vertex_Descriptor vi, vertices(m_poly_degrad))
  {
    put(cmdm_degrad, vi, 0.);
    put(tag_map, vi, ind);
    // put(cmdmm_degrad, vi, tab_cmdmm_degrad[ind]);
    ind++;
  }
}


template< typename PointMap,
          typename HalfedgeGraph,
          typename FaceNormalMap,
          typename VertexCMDMMap >
double
process_CMDM_multires(const HalfedgeGraph &m_poly_degrad,
                       const PointMap &pm_degrad,
                       const FaceNormalMap &fnm_degrad,
                       FEVV::PMapsContainer &pmaps_degrad,
                       const HalfedgeGraph &m_poly_original,
                       const PointMap &pm_original,
                       const FaceNormalMap &fnm_original,
                       FEVV::PMapsContainer &pmaps_original,
                       const int nb_level,
                       VertexCMDMMap &cmdm_pmap,
                       const double maxdim,
                       double &CMDM_value)
{
  auto curvature_map_degr =
      FEVV::make_vertex_property_map< HalfedgeGraph,
                                      FEVV::Filters::v_Curv_cmdm< HalfedgeGraph > >(
          m_poly_degrad);
  auto curvature_map_orig =
      FEVV::make_vertex_property_map< HalfedgeGraph,
                                      FEVV::Filters::v_Curv_cmdm< HalfedgeGraph > >(
          m_poly_original);

  CGAL::SM_Face_index *tab_matched_facet =
      new CGAL::SM_Face_index[size_of_vertices(m_poly_degrad)];


  VertexCMDMMap curv_match_pmap;

  if(has_map(pmaps_degrad, std::string("v:cmdm_match")))
  {
    curv_match_pmap =
        boost::any_cast< VertexCMDMMap >(pmaps_degrad.at("v:cmdm_match"));
  }
  else
  {
    curv_match_pmap =
        FEVV::make_vertex_property_map< HalfedgeGraph, double >(m_poly_degrad);
    pmaps_degrad["v:cmd_match"] = curv_match_pmap;
  }

  typedef typename FEVV::Vertex_pmap< HalfedgeGraph, CGALPoint >
      CMDM_nearest_map;
  CMDM_nearest_map cmdm_nearest_pmap;

  if(has_map(pmaps_degrad, std::string("v:cmdm_nearest")))
  {
    cmdm_nearest_pmap = boost::any_cast< CMDM_nearest_map >(
        pmaps_degrad.at("v:cmdm_nearest"));
  }
  else
  {
    cmdm_nearest_pmap =
        FEVV::make_vertex_property_map< HalfedgeGraph, CGALPoint >(
            m_poly_degrad);
    pmaps_degrad["v:cmdm_nearest"] = cmdm_nearest_pmap;
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

  /***Note that, static arrays were used to store feature values and vertices
  coordinates instead of proprty maps (PointMap, CMDMMap; e.g. CMDMMap
  &feature1), as the latter caused problems when parallelizing the code (notably
  the get() and put() functions): get(pm_degrad, vi) => tab_pm_degrad[i] ****/
  double *tab_f1 = new double[(int)num_vertices(m_poly_degrad)](); // initialized to 0
  double *tab_f2 = new double[(int)num_vertices(m_poly_degrad)]();
  double *tab_f3 = new double[(int)num_vertices(m_poly_degrad)]();
  double *tab_f4 = new double[(int)num_vertices(m_poly_degrad)]();
  double *tab_f5 = new double[(int)num_vertices(m_poly_degrad)]();
  double *tab_f6 = new double[(int)num_vertices(m_poly_degrad)]();
  double *tab_f7 = new double[(int)num_vertices(m_poly_degrad)]();
  double *tab_f8 = new double[(int)num_vertices(m_poly_degrad)]();

  CGALPoint *tab_pm_original =
      new CGALPoint[(int)num_vertices(m_poly_original)];
  CGALPoint *tab_pm_degrad = 
	  new CGALPoint[(int)num_vertices(m_poly_degrad)];
  CGALPoint *tab_cmdm_nearest_pmap =
      new CGALPoint[(int)num_vertices(m_poly_degrad)];

#pragma omp parallel for
  for(int ii = 0; ii < (int)num_vertices(m_poly_original); ii++)
  {
    Vertex_Descriptor vi(ii);
    auto p0 = get(pm_original, vi);
    tab_pm_original[ii] = p0;
  }

#pragma omp parallel for
  for(int ii = 0; ii < (int)num_vertices(m_poly_degrad); ii++)
  {
    Vertex_Descriptor vi(ii);
    auto p1 = get(pm_degrad, vi);
    tab_pm_degrad[ii] = p1;
  }

  matching_multires_init_cmdm(m_poly_degrad,
                              cmdm_pmap,
                              cmdm_nearest_pmap,
                              tag_map,
                              pm_degrad,
                              m_poly_original,
                              tab_matched_facet);

#pragma omp parallel for
  for(int ii = 0; ii < (int)num_vertices(m_poly_degrad); ii++)
  {
    Vertex_Descriptor vi(ii);
    auto p2 = get(cmdm_nearest_pmap, vi);
    tab_cmdm_nearest_pmap[ii] = p2;
  }

  // Initialize, from files, the 1D and 2D grids for LAB to LAB2000HL
  // conversion/interpolation
  static std::vector< double > init_grid_L;
  static std::vector< double > grid_L;
  static std::vector< std::vector< std::pair< double, double > > > init_grid_AB;
  initMatLABCH(init_grid_L, grid_L, init_grid_AB);

  using VertexColorMap = typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                                     MeshSurface >::pmap_type;
  VertexColorMap v_cm_orig = Get_VertexColorMap(
      m_poly_original, pmaps_original, init_grid_L, grid_L, init_grid_AB);
  VertexColorMap v_cm_degr = Get_VertexColorMap(
      m_poly_degrad, pmaps_degrad, init_grid_L, grid_L, init_grid_AB);
  VertexColorMap v_cm_match_pmap; // empty map which will be filled by the color
                                  // of the projected/matched points, in
                                  // matching_multires_update() function

  double radius_curvature = 0.001;
  for(int l = 0; l < nb_level; l++)
  {
    double min_nrm_min_curvature_deg, max_nrm_min_curvature_deg,
        min_nrm_max_curvature_deg, max_nrm_max_curvature_deg;
    double min_nrm_min_curvature_org, max_nrm_min_curvature_org,
        min_nrm_max_curvature_org, max_nrm_max_curvature_org;

    FEVV::Filters::calculate_curvature_cmdm(m_poly_degrad,
                                            curvature_map_degr,
                                            pm_degrad,
                                            fnm_degrad,
                                            true,
                                            maxdim * radius_curvature,
                                            min_nrm_min_curvature_deg,
                                            max_nrm_min_curvature_deg,
                                            min_nrm_max_curvature_deg,
                                            max_nrm_max_curvature_deg);
    FEVV::Filters::calculate_curvature_cmdm(m_poly_original,
                                            curvature_map_orig,
                                            pm_original,
                                            fnm_original,
                                            true,
                                            maxdim * radius_curvature,
                                            min_nrm_min_curvature_org,
                                            max_nrm_min_curvature_org,
                                            min_nrm_max_curvature_org,
                                            max_nrm_max_curvature_org);

    kmax_kmean_cmdm(m_poly_original, curvature_map_orig, maxdim);
    kmax_kmean_cmdm(m_poly_degrad, curvature_map_degr, maxdim);

    matching_multires_update(m_poly_degrad,
                             m_poly_original,
                             tab_pm_original,
                             tab_cmdm_nearest_pmap,
                             curvature_map_orig,
                             v_cm_orig,
                             curv_match_pmap,
                             v_cm_match_pmap,
                             tab_matched_facet);

    //#pragma omp parallel for
    for(int i = 0; i < (int)num_vertices(m_poly_degrad);
        ++i) // for(Vertex_Descriptor vi: vertices(m_poly_degrad))
    {
      // std::cout << "nb thread " << omp_get_num_threads() << std::endl;
      Vertex_Descriptor vi(i);

      std::vector< CGALPoint > tab_point1;
      std::vector< CGALPoint > tab_point2;

      std::vector< double > tab_distance1;
      std::vector< double > tab_distance2;

      std::vector< std::array< double, 3 > > tab_color1;
      std::vector< std::array< double, 3 > > tab_color2;

      process_cmdm_per_vertex(m_poly_degrad,
                              tab_pm_degrad,
                              tag_map,
                              tab_cmdm_nearest_pmap,
                              curvature_map_degr,
                              v_cm_degr,
                              curv_match_pmap,
                              v_cm_match_pmap,
                              vi,
                              radius_curvature * 3 * maxdim,
                              tab_point1,
                              tab_point2,
                              tab_distance1,
                              tab_distance2,
                              tab_color1,
                              tab_color2);

      compute_geometry_statistics(m_poly_degrad,
                                  pm_degrad,
                                  cmdm_nearest_pmap,
                                  vi,
                                  tab_distance1,
                                  tab_distance2,
                                  tab_point1,
                                  tab_point2,
                                  tab_f1,
                                  tab_f2,
                                  tab_f3,
                                  radius_curvature * 3 * maxdim);

      compute_color_statistics(m_poly_degrad,
                               pm_degrad,
                               cmdm_nearest_pmap,
                               vi,
                               tab_color1,
                               tab_color2,
                               tab_point1,
                               tab_point2,
                               tab_f4,
                               tab_f5,
                               tab_f6,
                               tab_f7,
                               tab_f8,
                               radius_curvature * 3 * maxdim);
    }
    radius_curvature += 0.0005;
  }

  double sum_f1 = 0, sum_f2 = 0, sum_f3 = 0, sum_f4 = 0, sum_f5 = 0, sum_f6 = 0,
         sum_f7 = 0, sum_f8 = 0;
  double w2 = 0.091, w5 = 0.22, w6 = 0.032,
         w7 = 0.656; // The recommended weights
  int nb_vertex = (int)num_vertices(m_poly_degrad);

  //#pragma omp parallel for
  for(int i = 0; i < nb_vertex; ++i)
  {
    //#pragma omp atomic
    sum_f1 += tab_f1[i];
    sum_f2 += tab_f2[i];
    sum_f3 += tab_f3[i];
    sum_f4 += tab_f4[i];
    sum_f5 += tab_f5[i];
    sum_f6 += tab_f6[i];
    sum_f7 += tab_f7[i];
    sum_f8 += tab_f8[i];

    // Local distortion
    double LD =
        (w2 * tab_f2[i] + w5 * tab_f5[i] + w6 * tab_f6[i] + w7 * tab_f7[i]) /
        nb_level;
    Vertex_Descriptor vi(i);
    put(cmdm_pmap, vi, LD);
  }

  sum_f1 = sum_f1 / (nb_vertex * nb_level);
  sum_f2 = sum_f2 / (nb_vertex * nb_level);
  sum_f3 = sum_f3 / (nb_vertex * nb_level);
  sum_f4 = sum_f4 / (nb_vertex * nb_level);
  sum_f5 = sum_f5 / (nb_vertex * nb_level);
  sum_f6 = sum_f6 / (nb_vertex * nb_level);
  sum_f7 = sum_f7 / (nb_vertex * nb_level);
  sum_f8 = sum_f8 / (nb_vertex * nb_level);

  CMDM_value = w2 * sum_f2 + w5 * sum_f5 + w6 * sum_f6 + w7 * sum_f7;

  delete[] tab_pm_original;
  delete[] tab_pm_degrad;
  delete[] tab_matched_facet;
  delete[] tab_cmdm_nearest_pmap;
  delete[] tab_f1;
  delete[] tab_f2;
  delete[] tab_f3;
  delete[] tab_f4;
  delete[] tab_f5;
  delete[] tab_f6;
  delete[] tab_f7;
  delete[] tab_f8;

  return 0;
}

/**
\brief Updates the matching process

\param	m_PolyDegrad	  	The first Polyhedron.
\param	_TabMatchedFacet	Facets from m_PolyOriginal on which vertices
from m_PolyDegrad are projected
*/

template< typename HalfedgeGraph,
          typename CurvatureVertexMap,
          typename CMDMCurvatureMatchMap,
          typename VertexColorMap >
void
matching_multires_update(
    const HalfedgeGraph &m_poly_degrad,
    const HalfedgeGraph &m_poly_orig,
    CGALPoint *tab_pm_orig,      // Or const PointMap &pm_orig,
    CGALPoint *tab_cmdmm_degrad, // Or CMDMNearestMap &cmdmm_degrad,
    const CurvatureVertexMap curv_orig,
    const VertexColorMap vcm_orig,
    CMDMCurvatureMatchMap curvmatch_pmap,
    VertexColorMap colmatch_pmap,
    CGAL::SM_Face_index *tab_matched_facet)
{
  using Color = typename boost::property_traits< VertexColorMap >::value_type;
  Color *tab_colmatch_pmap = new Color[(int)num_vertices(m_poly_degrad)];
  double *tab_curvmatch_pmap = new double[(int)num_vertices(m_poly_degrad)];

  //#pragma omp parallel for
  for(int i = 0; i < (int)num_vertices(m_poly_degrad);
      ++i) // BOOST_FOREACH(Vertex_Descriptor vi, vertices(m_poly_degrad))
  {
    // Vertex_Descriptor vi(i);
    int ind = i;
    auto *f_nearest = &tab_matched_facet[ind];

    /// calculation of the nearest point curvature value using vertices of the
    /// Nearest triangle 
	//we use linear interpolation using barycentric coordinates

    /* formulation used in NOT parallelized function
    auto x1 = get(pm_orig, target(halfedge(*f_nearest, m_poly_orig),
    m_poly_orig)); auto x2 = get(pm_orig, target(next(halfedge(*f_nearest,
    m_poly_orig), m_poly_orig), m_poly_orig)); auto x3 =get(pm_orig,
    target(next(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
    m_poly_orig), m_poly_orig));
    */

    auto *x1 =
        &tab_pm_orig[target(halfedge(*f_nearest, m_poly_orig), m_poly_orig)];
    auto *x2 = &tab_pm_orig[target(
        next(halfedge(*f_nearest, m_poly_orig), m_poly_orig), m_poly_orig)];
    auto *x3 = &tab_pm_orig[target(
        next(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig), m_poly_orig),
        m_poly_orig)];

    double l1 = sqrt((*x3 - *x2) * (*x3 - *x2));
    double l2 = sqrt((*x1 - *x3) * (*x1 - *x3));
    double l3 = sqrt((*x1 - *x2) * (*x1 - *x2));

    auto *pt_matched = &tab_cmdmm_degrad[i];
    auto v1 = *x1 - *pt_matched; // x1 - get(cmdmm_degrad, vi);
    auto v2 = *x2 - *pt_matched; // x2 - get(cmdmm_degrad, vi);
    auto v3 = *x3 - *pt_matched; // x3 - get(cmdmm_degrad, vi);

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

    std::array< double, 3 > col1 = {
        get(vcm_orig,
            target(halfedge(*f_nearest, m_poly_orig), m_poly_orig))[0],
        get(vcm_orig,
            target(halfedge(*f_nearest, m_poly_orig), m_poly_orig))[1],
        get(vcm_orig,
            target(halfedge(*f_nearest, m_poly_orig), m_poly_orig))[2]};
    std::array< double, 3 > col2 = {
        get(vcm_orig,
            target(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                   m_poly_orig))[0],
        get(vcm_orig,
            target(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                   m_poly_orig))[1],
        get(vcm_orig,
            target(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                   m_poly_orig))[2]};
    std::array< double, 3 > col3 = {
        get(vcm_orig,
            target(next(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                        m_poly_orig),
                   m_poly_orig))[0],
        get(vcm_orig,
            target(next(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                        m_poly_orig),
                   m_poly_orig))[1],
        get(vcm_orig,
            target(next(next(halfedge(*f_nearest, m_poly_orig), m_poly_orig),
                        m_poly_orig),
                   m_poly_orig))[2]};

    double curvmatch = 0.0;
    std::array< double, 3 > colmatch = {0, 0, 0};
    if((a1 + a2 + a3) > 0)
    {
      curvmatch = (a1 * c1 + a2 * c2 + a3 * c3) / (a1 + a2 + a3);
      colmatch = {(a1 * col1[0] + a2 * col2[0] + a3 * col3[0]) / (a1 + a2 + a3),
                  (a1 * col1[1] + a2 * col2[1] + a3 * col3[1]) / (a1 + a2 + a3),
                  (a1 * col1[2] + a2 * col2[2] + a3 * col3[2]) /
                      (a1 + a2 + a3)};
    }
    else
    {
      curvmatch = (c1 + c2 + c3) / 3;
      colmatch = {(col1[0] + col2[0] + col3[0]) / 3,
                  (col1[1] + col2[1] + col3[1]) / 3,
                  (col1[2] + col2[2] + col3[2]) / 3};
    }

    tab_curvmatch_pmap[i] = curvmatch; // put(curvmatch_pmap, vi, curvmatch);
    tab_colmatch_pmap[i] =
        Color(colmatch[0],
              colmatch[1],
              colmatch[2]); // put(colmatch_pmap, vi, Color(colmatch[0],
                            // colmatch[1], colmatch[2]));
  }

  int index = 0;
  BOOST_FOREACH(Vertex_Descriptor vi, vertices(m_poly_degrad))
  {
    put(curvmatch_pmap, vi, tab_curvmatch_pmap[index]);
    put(colmatch_pmap, vi, tab_colmatch_pmap[index]);
    index++;
  }

  delete[] tab_curvmatch_pmap;
  delete[] tab_colmatch_pmap;
}


/**
\brief  Computes the local neighborhoods

\param  pVertex  The considered vertex
\param  radius   radius of the neighborhood
\param [out] TabPoint1     3D points from the neighborhoodof pVertex regarding
the first polyhedron
\param [out] TabPoint2     3D points from the neighborhood of pVertex regarding
the second polyhedron
\param [out] TabDistance1  Curvature values from the neighborhood of pVertex
regarding the first polyhedron
\param [out] TabDistance2  Curvature values from the neighborhood of pVertex
regarding the second polyhedron
\param [out] TabColor1  Color values (L,A,B) from the neighborhood of pVertex
regarding the first polyhedron
\param [out] TabColor2  Color values (L,A,B) from the neighborhood of pVertex
regarding the second polyhedron
*/

template< typename VertexIterator,
          typename HalfedgeGraph,
          typename TagMap,
          typename CurvatureVertexMap,
          typename CMDMCurvatureMatchMap,
          typename VertexColorMap >
void
process_cmdm_per_vertex(
    const HalfedgeGraph &m_degr,
    CGALPoint *tab_pm_degr, // const PointMap &pm_degr,
    const TagMap &tags_pmap,
    CGALPoint *tab_nearest_pmap, // const CMDMNearestMap &nearest_pmap,
    const CurvatureVertexMap &curv_pmap_degr,
    const VertexColorMap &col_pmap_degr,
    const CMDMCurvatureMatchMap &curvmatch_pmap,
    const VertexColorMap &colmatch_pmap,
    const VertexIterator &p_vertex,
    double radius,
    std::vector< CGALPoint > &tab_point1,
    std::vector< CGALPoint > &tab_point2,
    std::vector< double > &tab_distance1,
    std::vector< double > &tab_distance2,
    std::vector< std::array< double, 3 > > &tab_color1,
    std::vector< std::array< double, 3 > > &tab_color2)
{
  std::set< int > vertices;

  std::stack< VertexIterator > s;

  auto *o = &tab_pm_degr[p_vertex]; // auto o = get(pm_degr, p_vertex);

  s.push(p_vertex);
  vertices.insert(get(tags_pmap, p_vertex));


  tab_point1.push_back(
      tab_pm_degr[p_vertex]); // tab_point1.push_back(get(pm_degr, p_vertex));
  tab_distance1.push_back(get(curv_pmap_degr, p_vertex).KmaxCurv);
  std::array< double, 3 > VertexColor1 = {(col_pmap_degr)[p_vertex][0],
                                          (col_pmap_degr)[p_vertex][1],
                                          (col_pmap_degr)[p_vertex][2]};
  tab_color1.push_back(VertexColor1);

  tab_point2.push_back(tab_nearest_pmap[p_vertex]); // tab_point2.push_back(get(nearest_pmap,
                                                    // p_vertex));
  tab_distance2.push_back(get(curvmatch_pmap, p_vertex));
  std::array< double, 3 > VertexColor2 = {(colmatch_pmap)[p_vertex][0],
                                          (colmatch_pmap)[p_vertex][1],
                                          (colmatch_pmap)[p_vertex][2]};
  tab_color2.push_back(VertexColor2);

  int nb_sommet_in_sphere = 0;

  while(!s.empty())
  {
    VertexIterator v_it = s.top();
    s.pop();
    CGALPoint *p = &tab_pm_degr[v_it]; // CGALPoint p = get(pm_degr, v_it);
    auto h = halfedge(v_it, m_degr);
    auto h0 = h;
    do
    {
      // A halfedge is directed from its source vertex to its target vertex
      auto v_target_of_h = target(h, m_degr);
      auto h_opp = opposite(h, m_degr);
      auto v_target_of_hopp = target(h_opp, m_degr);

      auto *p1 = &tab_pm_degr[v_target_of_h];    // auto p1 = get(pm_degr,
                                                 // target(h, m_degr));
      auto *p2 = &tab_pm_degr[v_target_of_hopp]; // auto p2 = get(pm_degr,
                                                 // target(opposite(h, m_degr),
                                                 // m_degr));

      auto *p1m =
          &tab_nearest_pmap[v_target_of_h]; // auto p1m = get(nearest_pmap,
                                            // target(h, m_degr));
      auto *p2m =
          &tab_nearest_pmap[v_target_of_hopp]; // auto p2m = get(nearest_pmap,
                                               // target(opposite(h, m_degr),
                                               // m_degr));

      auto v = (*p2 - *p1);
      auto vm = (*p2m - *p1m);

      if(v_it == p_vertex || v * (*p - *o) >= 0.0)
      {
        double len_old = std::sqrt(v * v);
        bool isect = sphere_clip_vector_cmdm(*o, radius, *p, v);
        double len_edge = std::sqrt(v * v);

        nb_sommet_in_sphere++;

        CGALPoint weighted_p1, weighted_p2;
        double weighted_curv1, weighted_curv2;
        std::array< double, 3 > weighted_col1, weighted_col2;

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
            weighted_p1 = *p1 + v;
            weighted_curv1 =
                (1 - len_edge / len_old) *
                    get(curv_pmap_degr, target(h, m_degr)).KmaxCurv +
                len_edge / len_old *
                    get(curv_pmap_degr, target(opposite(h, m_degr), m_degr))
                        .KmaxCurv;
            weighted_col1 = {
                (1 - len_edge / len_old) *
                        get(col_pmap_degr, target(h, m_degr))[0] +
                    len_edge / len_old *
                        get(col_pmap_degr,
                            target(opposite(h, m_degr), m_degr))[0],
                (1 - len_edge / len_old) *
                        get(col_pmap_degr, target(h, m_degr))[1] +
                    len_edge / len_old *
                        get(col_pmap_degr,
                            target(opposite(h, m_degr), m_degr))[1],
                (1 - len_edge / len_old) *
                        get(col_pmap_degr, target(h, m_degr))[2] +
                    len_edge / len_old *
                        get(col_pmap_degr,
                            target(opposite(h, m_degr), m_degr))[2]};

            weighted_p2 = *p1m + (len_edge / len_old) * vm;
            weighted_curv2 =
                (1 - len_edge / len_old) *
                    get(curvmatch_pmap, target(h, m_degr)) +
                len_edge / len_old *
                    get(curvmatch_pmap, target(opposite(h, m_degr), m_degr));
            weighted_col2 = {
                (1 - len_edge / len_old) *
                        get(colmatch_pmap, target(h, m_degr))[0] +
                    len_edge / len_old *
                        get(colmatch_pmap,
                            target(opposite(h, m_degr), m_degr))[0],
                (1 - len_edge / len_old) *
                        get(colmatch_pmap, target(h, m_degr))[1] +
                    len_edge / len_old *
                        get(colmatch_pmap,
                            target(opposite(h, m_degr), m_degr))[1],
                (1 - len_edge / len_old) *
                        get(colmatch_pmap, target(h, m_degr))[2] +
                    len_edge / len_old *
                        get(colmatch_pmap,
                            target(opposite(h, m_degr), m_degr))[2]};
          }
          else
          {
            weighted_p1 = *p2;
            weighted_p2 = *p2m;
            weighted_curv1 =
                get(curv_pmap_degr, target(opposite(h, m_degr), m_degr))
                    .KmaxCurv;
            weighted_curv2 =
                get(curvmatch_pmap, target(opposite(h, m_degr), m_degr));
            weighted_col1 = {
                get(col_pmap_degr, target(opposite(h, m_degr), m_degr))[0],
                get(col_pmap_degr, target(opposite(h, m_degr), m_degr))[1],
                get(col_pmap_degr, target(opposite(h, m_degr), m_degr))[2]};
            weighted_col2 = {
                get(colmatch_pmap, target(opposite(h, m_degr), m_degr))[0],
                get(colmatch_pmap, target(opposite(h, m_degr), m_degr))[1],
                get(colmatch_pmap, target(opposite(h, m_degr), m_degr))[2]};
          }

          tab_point1.push_back(weighted_p1);
          tab_distance1.push_back(weighted_curv1);
          tab_color1.push_back(weighted_col1);

          tab_point2.push_back(weighted_p2);
          tab_distance2.push_back(weighted_curv2);
          tab_color2.push_back(weighted_col2);
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
kmax_kmean_cmdm(const HalfedgeGraph &mesh,
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
\brief	Computes the multiscale CMDM metric

\param	m_PolyDegrad	  	The first Polyhedron.
\param	m_PolyOriginal	  	The second Polyhedron.
\param	NbLevel		Number of scales used
\param	maxdim	The max dimension of the Bounding Box Length.
\param [out] CMDMValue	The computed value of the CMDM metric
*/
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >,
          typename PointMap,
          typename FaceNormalMap,
          typename VertexCMDMMap >
double
process_CMDM_multires(const HalfedgeGraph &m_poly_degrad,
                       const PointMap &pm_degrad,
                       const FaceNormalMap &fnm_degrad,
                       FEVV::PMapsContainer &pmaps_degrad,
                       const HalfedgeGraph &m_poly_original,
                       const PointMap &pm_original,
                       const FaceNormalMap &fnm_original,
                       FEVV::PMapsContainer &pmaps_original,
                       const int nb_level,
                       VertexCMDMMap &cmdm_pmap,
                       double &CMDM_value)
{
  GeometryTraits gt_deg(m_poly_degrad);
  GeometryTraits gt_org(m_poly_original);

  return process_CMDM_multires(
      m_poly_degrad,
      pm_degrad,
      fnm_degrad,
      pmaps_degrad,
      m_poly_original,
      pm_original,
      fnm_original,
      pmaps_original,
      nb_level,
      cmdm_pmap,
      Filters::get_max_bb_size(m_poly_degrad, pm_degrad, gt_deg),
      CMDM_value);
}

} // namespace Filters
} // namespace FEVV
