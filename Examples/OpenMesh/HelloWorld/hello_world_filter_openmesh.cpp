// The following include is generic: it works without any change for
// all mesh datastructures. This should be true for all code to be
// encountered in directories mentioning "Generic" within their pathname.
#include "Examples/Generic/HelloWorld/hello_world_main.hpp"

// The following include defines the FEVV::MeshOpenMesh type that is
// used below when instantiating hello_world_main<>
#include "FEVV/DataStructures/DataStructures_openmesh.h"

// The filter itself (Examples/Generic/HelloWorld/hello_world_filter.hpp),
// as included by hello_world_main.hpp, manipulates both the geometry
// of the mesh as well as some of the mesh properties. The following
// includes respectively allow the filter to use expressions like
// `FEVV::Geometry_traits<...>` and
// `FEVV::Face_pmap<HalfedgeGraph, float> my_property_map;`.
#include "FEVV/Wrappings/Geometry_traits_openmesh.h"
#include "FEVV/Wrappings/properties_openmesh.h"


// Main: load a mesh, apply the filter, write the mesh
int
main(int argc, const char **argv)
{
  return hello_world_main< FEVV::MeshOpenMesh >(argc, argv);
}
