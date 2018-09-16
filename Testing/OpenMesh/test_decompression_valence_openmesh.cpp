#include "FEVV/DataStructures/DataStructures_openmesh.h"

#include "FEVV/Wrappings/Geometry_traits_openmesh.h"
#include "FEVV/Wrappings/properties_openmesh.h"

#include "Testing/Generic/Manifold/test_decompression_valence.inl"


int
main(int argc, const char **argv)
{
  return test_decompression_valence< FEVV::MeshOpenMesh >(argc, argv);
}
