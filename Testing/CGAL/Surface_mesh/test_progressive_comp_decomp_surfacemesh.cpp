#include "FEVV/DataStructures/DataStructures_cgal_surface_mesh.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"

#include "Testing/CGAL/test_progressive_comp_decomp.inl"

int
main(int argc, const char **argv)
{
  return test_automatic_progressive_compression_decompression< FEVV::MeshSurface >(argc, argv);
}
