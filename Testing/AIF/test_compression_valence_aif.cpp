#include "FEVV/DataStructures/DataStructures_aif.h"

#include "FEVV/Wrappings/Geometry_traits_aif.h"
#include "FEVV/Wrappings/properties_aif.h"

#include "Testing/Generic/Manifold/test_compression_valence.inl"


int
main(int argc, const char **argv)
{
  return test_compression_valence< FEVV::MeshAIF >(argc, argv);
}
