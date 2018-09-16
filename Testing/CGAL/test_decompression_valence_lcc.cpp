#include "FEVV/DataStructures/DataStructures_cgal_linear_cell_complex.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_linear_cell_complex.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"

#include "Testing/Generic/Manifold/test_decompression_valence.inl"


int
main(int argc, const char **argv)
{
  return test_decompression_valence< FEVV::MeshLCC >(argc, argv);
}
