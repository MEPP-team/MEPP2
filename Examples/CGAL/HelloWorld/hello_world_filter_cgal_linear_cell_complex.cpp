


#include "Examples/Generic/HelloWorld/hello_world_main.hpp"

// The following include defines the FEVV::MeshLCC type
#include "FEVV/DataStructures/DataStructures_cgal_linear_cell_complex.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_linear_cell_complex.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"


// Main: load a mesh, apply the filter, write the mesh
int
main(int argc, const char **argv)
{
  return hello_world_main< FEVV::MeshLCC >(argc, argv);
}
