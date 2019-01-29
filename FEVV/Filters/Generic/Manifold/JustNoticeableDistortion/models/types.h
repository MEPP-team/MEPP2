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

#include <Eigen/Core>
using namespace Eigen;

//-----------------------------------------------------------------------------------

typedef Eigen::Vector3d LightType;
typedef Eigen::Vector3d CamType;

//-----------------------------------------------------------------------------------
/**
\struct ScreenParam
\brief The parameters of the screen
*/
struct ScreenParam
{
  ScreenParam(int h = 0, int v = 0, double d = 0.) : hres(h), vres(v), diag(d)
  {
  }

  int hres;    /** < Horizontal resolution in pixel*/
  int vres;    /** < Vertical resolution in pixel*/
  double diag; /** < Diagonal in cm*/
};

//-----------------------------------------------------------------------------------
/**
\struct UserParam
\brief The parameters of the user
*/
struct UserParam
{
  UserParam(double d = 0) : dist(d) {}

  double dist; /** < Distance between the screen and the user in cm*/
};

//-----------------------------------------------------------------------------------
/**
\struct SceneParam
\brief The parameters of the scene
*/
struct SceneParam
{
  SceneParam(int h = 0, double f = 0.) : hvport(h), fov(f) {}

  int hvport; /** < Width of the scene in pixels*/
  double fov; /** < Field of view of the scene in radian*/
};

//-----------------------------------------------------------------------------------
