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

inline
int
cpd_to_pixel(int screen_resX,
             int screen_resY,
             double screen_size,
             double user_distance,
             double cpd)
{
  double ppcm = sqrt(double(screen_resX * screen_resX) +
                     double(screen_resY * screen_resY)) /
                screen_size;
  double ppd = 0.0175 * user_distance * ppcm;
  int nb_pixel = int(ppd / cpd);

  return nb_pixel;
}

inline
double
pixel_to_wDistance(int camera_vportY,
                   double camera_fov,
                   double camera_distance,
                   int nb_pixel)
{
  double dpp = 2.f * camera_distance * std::tan(0.5 * camera_fov) /
               double(camera_vportY);
  double wDistance = dpp * double(nb_pixel);

  return wDistance;
}

inline
double
cpd_to_wDistance(int screen_resX,
                 int screen_resY,
                 double screen_size,
                 double user_distance,
                 int camera_vportY,
                 double camera_fov,
                 double camera_distance,
                 double cpd)
{
  return pixel_to_wDistance(
      camera_vportY,
      camera_fov,
      camera_distance,
      cpd_to_pixel(screen_resX, screen_resY, screen_size, user_distance, cpd));
}

inline
double
pixel_to_cpd(int screen_resX,
             int screen_resY,
             double screen_size,
             double user_distance,
             int nb_pixel)
{
  double ppcm = std::sqrt(double(screen_resX * screen_resX) +
                          double(screen_resY * screen_resY)) /
                screen_size;
  double ppd = 0.0175 * user_distance * ppcm;
  double cpd = ppd / double(nb_pixel);

  return cpd;
}

inline
int
wDistance_to_pixel(int camera_vportY,
                   double camera_fov,
                   double camera_distance,
                   double wDistance)
{
  double dpp = 2.f * camera_distance * std::tan(0.5 * camera_fov) /
               double(camera_vportY);
  int nb_pixel = int(wDistance / dpp);
  return nb_pixel;
}

inline
double
wDistance_to_cpd(int screen_resX,
                 int screen_resY,
                 double screen_size,
                 double user_distance,
                 int camera_vportY,
                 double camera_fov,
                 double camera_distance,
                 double wDistance)
{
  return pixel_to_cpd(
      screen_resX,
      screen_resY,
      screen_size,
      user_distance,
      wDistance_to_pixel(
          camera_vportY, camera_fov, camera_distance, wDistance));
}
