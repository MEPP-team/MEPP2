// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
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

#include <time.h>
#include <boost/foreach.hpp>

#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/genericparametricmodel.h"
#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/contrastmasking.h"
#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/contrastsensitivity.h"
#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/psychometricfunction.h"
#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/threshold.h"
#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/types.h"
#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/visibility.h"

#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/perception/flatcontrastcomputor.h"
#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/perception/flatfrequencycomputor.hpp"
#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/utils/lightsampler.hpp"
#include "FEVV/Filters/Generic/AABB/compute_mesh_bounding_box.hpp"

//---------------------------------------------------------
//                  Filter definition
//---------------------------------------------------------

namespace FEVV {
namespace Filters {

/**
\brief	Compute the visibility of the movement of an halfedge, given a light
source and the displacement
*/
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits,
          typename FaceNormalMap,
          typename HalfEdge,
          typename CamType,
          typename LightType,
          typename DirType >
double
compute_visibility(const GeometryTraits &geom,
                   const PointMap &point_map,
                   const HalfedgeGraph &mesh,
                   const FaceNormalMap &f_nm,
                   const HalfEdge &halfedge,
                   const DirType &displacement,
                   const LightType &light,
                   const CamType &cam,
                   const ScreenParam &screen,
                   const UserParam &user,
                   const SceneParam &scene,
                   const NWHWD16_Threshold &threshold,
                   bool both_faces_move = true,
                   bool data_output = false)
{
  if(face(halfedge, mesh) ==
         boost::graph_traits< HalfedgeGraph >::null_face() ||
     face(opposite(halfedge, mesh), mesh) ==
         boost::graph_traits< HalfedgeGraph >::null_face())
    return 0.0;

  double init_contrast =
      compute_flat_contrast(geom, light, mesh, f_nm, halfedge);
  double init_frequency = compute_flat_frequency(
      cam, point_map, mesh, geom, halfedge, screen, user, scene);

  double visibility = 0.0;
  typename GeometryTraits::Vector n1;
  typename GeometryTraits::Vector n2;
  double new_contrast;
  double new_frequency;
  if(both_faces_move)
  {
    /*
    3 - - 2
      \   | \
        \ |   \
          1 - - 4
  1 is moving
    */
    auto v1 =
        get(point_map, source(halfedge, mesh)); // the vertex that is moving
    auto v2 = get(point_map, target(halfedge, mesh));
    auto v3 = get(point_map, target(next(halfedge, mesh), mesh));
    auto v4 =
        get(point_map, target(next(opposite(halfedge, mesh), mesh), mesh));

    typename GeometryTraits::Point p1(
        geom.get_x(v1), geom.get_y(v1), geom.get_z(v1));
    typename GeometryTraits::Point p2(
        geom.get_x(v2), geom.get_y(v2), geom.get_z(v2));
    typename GeometryTraits::Point p3(
        geom.get_x(v3), geom.get_y(v3), geom.get_z(v3));
    typename GeometryTraits::Point p4(
        geom.get_x(v4), geom.get_y(v4), geom.get_z(v4));
    typename GeometryTraits::Vector _v3(
        geom.get_x(v3), geom.get_y(v3), geom.get_z(v3));
    typename GeometryTraits::Vector _v4(
        geom.get_x(v4), geom.get_y(v4), geom.get_z(v4));

    n1 = geom.unit_normal(p1 + displacement, p2, p3);
    n2 = geom.unit_normal(p1 + displacement, p4, p2);

    new_contrast = compute_flat_contrast(geom, n1, n2, light);
    new_frequency =
        compute_flat_frequency(cam, _v3, _v4, geom, screen, user, scene);
  }
  else
  {
    /*
    1 - - 2
      \   | \
        \ |   \
          3 - - 4
  1 is moving
    */
    auto v1 = get(point_map, target(next(halfedge, mesh), mesh));
    auto v2 = get(point_map, source(halfedge, mesh));
    auto v3 = get(point_map, target(halfedge, mesh));
    auto v4 =
        get(point_map, target(next(opposite(halfedge, mesh), mesh), mesh));

    typename GeometryTraits::Point p1(
        geom.get_x(v1), geom.get_y(v1), geom.get_z(v1));
    typename GeometryTraits::Point p2(
        geom.get_x(v2), geom.get_y(v2), geom.get_z(v2));
    typename GeometryTraits::Point p3(
        geom.get_x(v3), geom.get_y(v3), geom.get_z(v3));
    typename GeometryTraits::Point p4(
        geom.get_x(v4), geom.get_y(v4), geom.get_z(v4));
    typename GeometryTraits::Vector _v1(
        geom.get_x(v1), geom.get_y(v1), geom.get_z(v1));
    typename GeometryTraits::Vector _v4(
        geom.get_x(v4), geom.get_y(v4), geom.get_z(v4));

    n1 = geom.unit_normal(p2, p3, p1 + displacement);
    n2 = geom.unit_normal(p2, p4, p3);

    new_contrast = compute_flat_contrast(geom, n1, n2, light);
    new_frequency = compute_flat_frequency(
        cam, _v1 + displacement, _v4, geom, screen, user, scene);
  }


  double T1 = threshold(NWHWD16_Threshold::InputType(
      init_contrast, std::max(.1, init_frequency)));
  double T2 = threshold(
      NWHWD16_Threshold::InputType(init_contrast, std::max(.1, new_frequency)));

  VisibilityModel visibility_m;

  bool is_ambigous = !(new_contrast < 0.) != !(init_contrast < 0.);

  double dc = is_ambigous ? fabs(init_contrast) + fabs(new_contrast)
                          : fabs(init_contrast - new_contrast);

  visibility = std::max(visibility_m(VisibilityModel::InputType(dc, T1)),
                        visibility_m(VisibilityModel::InputType(dc, T2)));

  return visibility;
}

/**
\brief	Compute the JND value for a vertex, using a given light
*/
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits,
          typename FaceNormalMap,
          typename VertexType,
          typename CamType,
          typename LightType,
          typename DirType >
double
compute_threshold(const GeometryTraits &geom,
                  const PointMap &point_map,
                  const HalfedgeGraph &mesh,
                  const FaceNormalMap &f_nm,
                  const VertexType &vertex,
                  const DirType &dir,
                  const LightType &light,
                  const CamType &cam,
                  const ScreenParam &screen,
                  const UserParam &user,
                  const SceneParam &scene,
                  const NWHWD16_Threshold &threshold,
                  double max_displacement,
                  bool data_output = false)
{
  double tolerence = 0.80;
  double precision = 0.1;
  double int_width = 1.e-8;

  double a = 0.0;
  double b = max_displacement;

  double T = b;
  double visibility = 0.0;

  auto h0 = next(halfedge(vertex, mesh), mesh);
  auto h = h0;

  // Computation of visibility for the maximum displacement
  do
  {
    double visibility_b = 0.0;
    if(face(h, mesh) != boost::graph_traits< HalfedgeGraph >::null_face())
      visibility_b = compute_visibility(geom,
                                        point_map,
                                        mesh,
                                        f_nm,
                                        h,
                                        T * dir,
                                        light,
                                        cam,
                                        screen,
                                        user,
                                        scene,
                                        threshold,
                                        true,
                                        data_output);
    visibility = std::max(visibility, visibility_b);
    if(visibility >= 1.0)
      break;
    if(face(next(h, mesh), mesh) !=
       boost::graph_traits< HalfedgeGraph >::null_face())
      visibility_b = compute_visibility(geom,
                                        point_map,
                                        mesh,
                                        f_nm,
                                        next(h, mesh),
                                        T * dir,
                                        light,
                                        cam,
                                        screen,
                                        user,
                                        scene,
                                        threshold,
                                        false,
                                        data_output);
    visibility = std::max(visibility, visibility_b);
    if(visibility >= 1.0)
      break;
    h = next(opposite(h, mesh), mesh);
  } while(h0 != h);

  // Computation of the visibility until the tolerence is met
  while(fabs(visibility - tolerence) > precision)
  {
    T = (b + a) * 0.5;
    visibility = 0.0;
    h = h0;
    do
    {
      double visibility_b = 0.0;
      if(face(h, mesh) != boost::graph_traits< HalfedgeGraph >::null_face())
        visibility_b = compute_visibility(geom,
                                          point_map,
                                          mesh,
                                          f_nm,
                                          h,
                                          T * dir,
                                          light,
                                          cam,
                                          screen,
                                          user,
                                          scene,
                                          threshold,
                                          true,
                                          data_output);
      visibility = std::max(visibility, visibility_b);
      if(visibility >= 1.0)
        break;
      if(face(next(h, mesh), mesh) !=
         boost::graph_traits< HalfedgeGraph >::null_face())
        visibility_b = compute_visibility(geom,
                                          point_map,
                                          mesh,
                                          f_nm,
                                          next(h, mesh),
                                          T * dir,
                                          light,
                                          cam,
                                          screen,
                                          user,
                                          scene,
                                          threshold,
                                          false,
                                          data_output);
      visibility = std::max(visibility, visibility_b);
      if(visibility >= 1.0)
        break;
      h = next(opposite(h, mesh), mesh);
    } while(h0 != h);

    if(visibility > tolerence)
      b = T;
    else
      a = T;

    if((b - a) < int_width)
      break;
  }

  return T;
}


/**
\brief  Computes the Just Noticeable Distortion metric

\param  halfedge_graph  The mesh whose JND is computed
\param  oint_map  The Point map of g
\param  vertex_normal_map  The vertices normal map of g
\param  geometry_traits  The geometry traits of g
\param  ace_map  The Face map of g
\param  face_normal_map  The faces normal map of g
\param  jnd_map  The JND map of g (a map from vertex to
        double)
\param  screen  The screen parameters used (included in models/types.h)
\param  user  The user parameters used (included in models/types.h)
\param  scene  The scene parameters used (included in models/types.h)
\param  nb_light_sources  The number of randomly generated lights around
        each vertex
\param  use_random  Should we use a randomness to generate the lights ?
*/
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >,
          typename FaceNormalMap,
          typename JNDMap >
void
just_noticeable_distortion_filter(const HalfedgeGraph &halfedge_graph,
                                  PointMap &point_map,
                                  const VertexNormalMap &vertex_normal_map,
                                  const GeometryTraits &geometry_traits,
                                  const FaceNormalMap &face_normal_map,
                                  JNDMap &jnd_map,
                                  const ScreenParam &screen,
                                  const UserParam &user,
                                  const SceneParam &scene,
                                  const int nb_light_sources,
                                  bool data_output = false,
                                  bool use_random = true)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_iterator face_iterator;
  typedef typename boost::property_traits< PointMap >::value_type Point;
  typedef typename boost::property_traits< VertexNormalMap >::value_type Normal;
  typedef
      typename boost::property_traits< FaceNormalMap >::value_type FaceNormal;

  int counter = 0;

  // Preparation of the model defined in G. Nader thesis
  SarkisonCSF csf(SarkisonCSF::ParameterType(-15.13, 0.0096, 0.64));
  DalyMasking masking(DalyMasking::ParameterType(0.0078, 88.29, 1.0, 4.207));

  NWHWD16_Threshold threshold;
  threshold.param().csf = csf;
  threshold.param().masking = masking;

  // Computation of the bounding box of the mesh
  typename boost::property_traits< PointMap >::value_type minAABB, maxAABB;
  Filters::compute_mesh_bounding_box< HalfedgeGraph, PointMap, GeometryTraits >(
      halfedge_graph, point_map, minAABB, maxAABB, geometry_traits);

  typename GeometryTraits::Vector minbb(geometry_traits.get_x(minAABB),
                                        geometry_traits.get_y(minAABB),
                                        geometry_traits.get_z(minAABB));
  typename GeometryTraits::Vector maxbb(geometry_traits.get_x(maxAABB),
                                        geometry_traits.get_y(maxAABB),
                                        geometry_traits.get_z(maxAABB));

  typename GeometryTraits::Vector bbox_c = (maxbb + minbb) * 0.5;
  double bbox_diag_len = geometry_traits.length(maxbb - minbb);

  // Creation of a camera position (should be externalized later)
  typename GeometryTraits::Vector cam =
      bbox_c + typename GeometryTraits::Vector(0., 0., 1.) * bbox_diag_len;
#ifndef FEVV_USE_AIF
  if(data_output)
    std::cout << "Cam:" << cam << std::endl;
#endif

  LightSampler lsampler;
  Eigen::MatrixX3d lights;

  // preparation of some metrics for avancement
  int avance = 0;
  int mesh_size = static_cast< int >(size_of_vertices(halfedge_graph));
  int step = mesh_size * nb_light_sources < 120000 ? 10 : 1;

  clock_t start = clock();

  BOOST_FOREACH(vertex_descriptor vi, vertices(halfedge_graph))
  {
    // Display of estimated computation time and advancement of computation
    if(!data_output && (counter > (avance + step) * (mesh_size / 100)))
    {
      if(avance == 0)
      {
        double estimated_time =
            static_cast< double >(clock() - start) /
            (static_cast< double >(CLOCKS_PER_SEC * 100) / step);
        std::cout << "Estimated time of computation : " << estimated_time
                  << " s for " << mesh_size << " vertices." << std::endl;
      }
      avance += step;
      std::cout << "Progression: " << avance << "%\r" << std::flush;
    }
    counter++;


    Point p = get(point_map, vi);
    Normal n = get(vertex_normal_map, vi);

    halfedge_descriptor h0 = next(halfedge(vi, halfedge_graph), halfedge_graph);
    halfedge_descriptor h = h0;

#ifndef FEVV_USE_AIF
    if(data_output)
    {
      std::cout << "processing vertex #" << counter << " : "
                << geometry_traits.get_x(p) << " " << geometry_traits.get_y(p)
                << " " << geometry_traits.get_z(p) << std::endl;
      std::cout << "\tnormal :" << n[0] << " " << n[1] << " " << n[2]
                << std::endl;
    }
#endif

    Eigen::Vector3d p_n(n[0], n[1], n[2]);
    double jnd = 0.0;

    // Computation of average angle between vertex normal and adjacent faces
    // normals
    double avg_angle = 0.0;
    double valence = 0.0;
    do
    {
      if(face(h, halfedge_graph) !=
         boost::graph_traits< HalfedgeGraph >::null_face())
      {
        auto face_normal = get(face_normal_map, face(h, halfedge_graph));
        Eigen::Vector3d f_n(face_normal[0], face_normal[1], face_normal[2]);
        valence += 1.0;
        avg_angle += std::acos(p_n.normalized().dot(f_n.normalized()));
      }
      h = next(opposite(h, halfedge_graph), halfedge_graph);
    } while(h0 != h);

    avg_angle = avg_angle / valence;
    avg_angle =
        (avg_angle != avg_angle)
            ? M_PI_4 * 0.5
            : M_PI_2 - avg_angle; // (avg_angle != avg_angle) instead of isnan

#ifndef FEVV_USE_AIF
    if(data_output)
      std::cout << "\tMean angle: " << avg_angle << std::endl;
#endif
    // Creation of random light sources around the vertex
    lsampler.sample_to_global(lights,
                              nb_light_sources,
                              p_n,
                              0.80 * avg_angle,
                              0.85 * avg_angle,
                              use_random);

    // initialization of threshold at 10% of the bounding box size
    VectorXd T;
    T.resize(lights.rows());
    T.setConstant(0.1 * bbox_diag_len);

    int count = 0;
    for(unsigned int i = 0; i < lights.rows(); ++i)
    {
      typename GeometryTraits::Vector light(
          lights.row(i)[0], lights.row(i)[1], lights.row(i)[2]);
      if(geometry_traits.dot_product(light, n) > 0.)
      {
        count++;
        T(i) = compute_threshold(geometry_traits,
                                 point_map,
                                 halfedge_graph,
                                 face_normal_map,
                                 vi,
                                 n,
                                 light,
                                 cam,
                                 screen,
                                 user,
                                 scene,
                                 threshold,
                                 0.1 * bbox_diag_len,
                                 data_output);
#ifndef FEVV_USE_AIF
        if(data_output)
          std::cout << "\t" << light << " -> " << T(i) << std::endl;
#endif
      }
    }

#ifndef FEVV_USE_AIF
    if(data_output)
      std::cout << "\tUsed lights :" << count << std::endl;
#endif

    // if no lights were useable to compute the JND set it to 0
    // in the other case keep the lowest threshold
    if(count != 0)
      jnd = T.minCoeff();

    if(jnd != jnd) // instead of isnan
      jnd = 0.0;

    put(jnd_map, vi, jnd);
  }
  clock_t end = clock();
  std::cout << "Average jnd computing time : "
            << static_cast< double >(end - start) /
                   static_cast< double >(CLOCKS_PER_SEC * counter)
            << " s/vertex" << std::endl;
}


/**
\brief	Computes the Just Noticeable Distortion metric

\param  g           The mesh whose JND is computed
\param  pm          The Point map of g
\param  v_nm        The vertices normal map of g
\param  f_nm        The faces normal map of g
\param  jnd_m       The JND map of g (a map from vertex to double)
\param  screen      The screen parameters used (included in models/types.h)
\param  user        The user parameters used (included in models/types.h)
\param  scene       The scene parameters used (included in models/types.h)
\param  lights      The number of randomly generated lights around each vertex
\param  use_random  Should we use a randomness to generate the lights ?
*/
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph >,
          typename FaceNormalMap,
          typename JNDMap >
void
just_noticeable_distortion_filter(const HalfedgeGraph &g,
                                  PointMap &pm,
                                  const VertexNormalMap &v_nm,
                                  const FaceNormalMap &f_nm,
                                  JNDMap &jnd_m,
                                  const ScreenParam &screen,
                                  const UserParam &user,
                                  const SceneParam &scene,
                                  const int lights = 128,
                                  bool data_output = false,
                                  bool use_random = true)
{
  GeometryTraits gt(g);
  just_noticeable_distortion_filter(g,
                                    pm,
                                    v_nm,
                                    gt,
                                    f_nm,
                                    jnd_m,
                                    screen,
                                    user,
                                    scene,
                                    lights,
                                    data_output,
                                    use_random);
}
} // namespace Filters
} // namespace FEVV
