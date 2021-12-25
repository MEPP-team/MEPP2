#include "FEVV/DataStructures/DataStructures_cgal_polyhedron_3.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/properties_polyhedron_3.h"

#include "FEVV/Tools/Comparator/SpanningTreeVertexEdgeComparator.hpp" 

#include "FEVV/Tools/IO/FileUtilities.hpp"

#include <fstream>
#include <string>

using namespace FEVV;
using namespace FEVV::Comparator;

void
testSpanningTreeComparatorSurfaceMesh(const std::string& filename)
{
  std::ifstream in(filename);
  if(!in)
  {
    std::cout << "Unable to read file " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  typedef FEVV::MeshPolyhedron Mesh;
  /**********************************************************************************************************/
  Mesh m;
  in >> m;
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
