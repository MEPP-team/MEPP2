#include "FEVV/DataStructures/DataStructures_cgal_surface_mesh.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"

#include "Testing/Generic/Manifold/test_curvature_filter.inl"


int
main(int argc, const char **argv)
{
  return test_curvature_filter< FEVV::MeshSurface >(argc, argv);
}
