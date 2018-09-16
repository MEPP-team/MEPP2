#ifndef JNDModels_TYPES_H
#define JNDModels_TYPES_H

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

#endif // JNDModels_TYPES_H
