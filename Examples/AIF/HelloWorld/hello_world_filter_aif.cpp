#include "FEVV/DataStructures/AIF/AIFMesh.hpp"

#include "FEVV/Wrappings/Geometry_traits_aif.h"
#include "FEVV/Wrappings/properties_aif.h"

/// \todo Interestingly enough the following inclusion must be placed _after_
/// the AIF specific include files for the file to compile.
/// This behaviour is different from the ones of the equivalent examples
/// using other data structures like OpenMesh or the ones of CGAL.
/// Inquire on this (peculiar?) behavioral difference of AIF.
#include "Examples/Generic/HelloWorld/hello_world_main.hpp"

int
main(int argc, const char **argv)
{
  return hello_world_main< FEVV::MeshAIF >(argc, argv);
}
