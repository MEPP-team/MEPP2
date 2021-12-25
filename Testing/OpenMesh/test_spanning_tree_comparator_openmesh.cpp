#include "FEVV/DataStructures/DataStructures_openmesh.h"

#include "FEVV/Wrappings/Geometry_traits_openmesh.h"
#include "FEVV/Wrappings/properties_openmesh.h"

#include <OpenMesh/Core/IO/MeshIO.hh>

#include "FEVV/Tools/Comparator/SpanningTreeVertexEdgeComparator.hpp" 

#include "FEVV/Tools/IO/FileUtilities.hpp"

#include <string>

using namespace FEVV;
using namespace FEVV::Comparator;

void
testSpanningTreeComparatorSurfaceMesh(const std::string& filename)
{
  typedef FEVV::MeshOpenMesh Mesh;
  /**********************************************************************************************************/
  Mesh m;
  if(!OpenMesh::IO::read_mesh(m, filename))
  {
    std::cout << "Unable to read file " << filename << std::endl;
    return;
  }
  auto pos_pm = get(boost::vertex_point, m);
  /**********************************************************************************************************/
  auto st = get_spanning_tree_comparator(m, pos_pm, false);
}

int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: " << argv[0]
              << " filename; filename being an off file." << std::endl;
    exit(EXIT_FAILURE);
  }

  testSpanningTreeComparatorSurfaceMesh(argv[1]);
  return 0;
}
