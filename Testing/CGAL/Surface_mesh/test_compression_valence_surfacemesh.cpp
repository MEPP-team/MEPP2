#include "FEVV/DataStructures/DataStructures_cgal_surface_mesh.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"

#include "Testing/Generic/Manifold/test_compression_valence.inl"


int
main(int argc, const char **argv)
{
  return test_compression_valence< FEVV::MeshSurface >(argc, argv);
}
