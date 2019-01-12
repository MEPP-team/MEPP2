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

#include <cmath>
#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/types.h"

int
cpd_to_pixel(int screen_resX,
             int screen_resY,
             double screen_size,
             double user_distance,
             double cpd);

double
pixel_to_wDistance(int camera_vportY,
                   double camera_fov,
                   double camera_distance,
                   int nb_pixel);

double
cpd_to_wDistance(int screen_resX,
                 int screen_resY,
                 double screen_size,
                 double user_distance,
                 int camera_vportY,
                 double camera_fov,
                 double camera_distance,
                 double cpd);

double
pixel_to_cpd(int screen_resX,
             int screen_resY,
             double screen_size,
             double user_distance,
             int nb_pixel);

int
wDistance_to_pixel(int camera_vportY,
                   double camera_fov,
                   double camera_distance,
                   double wDistance);

double
wDistance_to_cpd(int screen_resX,
                 int screen_resY,
                 double screen_size,
                 double user_distance,
                 int camera_vportY,
                 double camera_fov,
                 double camera_distance,
                 double wDistance);

template< typename Cam_t, typename GeometryTraits, typename Point_t >
double
compute_flat_frequency(const Cam_t &cam,
                       const Point_t &p1,
                       const Point_t &p2,
                       const GeometryTraits &geom_t,
                       const ScreenParam &screen,
                       const UserParam &user,
                       const SceneParam &scene)
{

  double cdist = geom_t.length((0.5 * (p1 + p2) - cam));

  double length = geom_t.length(p1 - p2);
  // convert length to cpd
  int npx = wDistance_to_pixel(scene.hvport, scene.fov, cdist, length);

  double cpd =
      pixel_to_cpd(screen.hres, screen.vres, screen.diag, user.dist, npx);

  double max_cpd =
      pixel_to_cpd(screen.hres, screen.vres, screen.diag, user.dist, 2);

  return std::min(cpd, max_cpd);
}

template< typename Cam_t,
          typename PointMap,
          typename HalfedgeGraph,
          typename GeometryTraits,
          typename HalfEdge >
double
compute_flat_frequency(const Cam_t &cam,
                       const PointMap &point_map,
                       const HalfedgeGraph &mesh,
                       const GeometryTraits &geom_t,
                       const HalfEdge &halfedge,
                       const ScreenParam &screen,
                       const UserParam &user,
                       const SceneParam &scene)
{
  auto opp = opposite(halfedge, mesh);
  auto v1 = get(point_map, target(next(halfedge, mesh), mesh));
  auto v2 = get(point_map, target(next(opp, mesh), mesh));

  typename GeometryTraits::Vector p1(
      geom_t.get_x(v1), geom_t.get_y(v1), geom_t.get_z(v1));
  typename GeometryTraits::Vector p2(
      geom_t.get_x(v2), geom_t.get_y(v2), geom_t.get_z(v2));

  return compute_flat_frequency(cam, p1, p2, geom_t, screen, user, scene);
}

#include "flatfrequencycomputor.inl"
