#ifndef FEVV_DATASTRUCTURES_OPENMESH_H
#define FEVV_DATASTRUCTURES_OPENMESH_H


#include <OpenMesh/Core/IO/MeshIO.hh>

// force usage of OpenMesh native Point type instead of CGAL Point type
// in point map in properties_PolyMesh_ArrayKernelT.h
#define CGAL_USE_OM_POINTS
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

namespace FEVV {

/* From :
http://www.openmesh.org/media/Documentations/OpenMesh-Doc-Latest/a00012.html

--> The default traits class looks like this:
struct DefaultTraits
{
  typedef Vec3f  Point;
  typedef Vec3f  Normal;
  typedef Vec2f  TexCoord;
  typedef Vec3uc Color;
  VertexTraits    {};
  HalfedgeTraits  {};
  EdgeTraits      {};
  FaceTraits      {};

  VertexAttributes(0);
  HalfedgeAttributes(Attributes::PrevHalfedge);
  EdgeAttributes(0);
  FaceAttributes(0);
};*/

struct MyTraits : public OpenMesh::DefaultTraits
{
  typedef OpenMesh::Vec3d Point;    // use double-values
  typedef OpenMesh::Vec3d Normal;   // use double-values
  typedef OpenMesh::Vec2d TexCoord; // use double-values
};

using MeshOpenMesh = OpenMesh::PolyMesh_ArrayKernelT< MyTraits >;

} // namespace FEVV


#endif // FEVV_DATASTRUCTURES_OPENMESH_H
