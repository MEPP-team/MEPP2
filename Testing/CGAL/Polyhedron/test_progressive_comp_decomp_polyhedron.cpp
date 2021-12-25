#include "FEVV/DataStructures/DataStructures_cgal_polyhedron_3.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/properties_polyhedron_3.h"

#include "Testing/CGAL/test_progressive_comp_decomp.inl"

int
main(int argc, const char **argv)
{
  return test_automatic_progressive_compression_decompression< FEVV::MeshPolyhedron >(argc, argv);
}
