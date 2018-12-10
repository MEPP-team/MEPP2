


#include "Examples/Generic/HelloWorld/hello_world_main.hpp"

// The following include defines the FEVV::MeshPolyhedron type
#include "FEVV/DataStructures/DataStructures_cgal_polyhedron_3.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/properties_polyhedron_3.h"


// Main: load a mesh, apply the filter, write the mesh
int
main(int argc, const char **argv)
{
  return hello_world_main< FEVV::MeshPolyhedron >(argc, argv);
}
