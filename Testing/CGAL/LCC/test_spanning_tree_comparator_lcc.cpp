#include "FEVV/DataStructures/DataStructures_cgal_linear_cell_complex.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_linear_cell_complex.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"

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

  typedef FEVV::MeshLCC Mesh;
  /**********************************************************************************************************/
  Mesh m;
  try
  {
    CGAL::read_off(in, m);
  }
  catch (const std::length_error &le)
  {
    std::cerr << "[LCC] Exception caught while reading input file " << filename
      << ": " << le.what() << std::endl;
    BOOST_ASSERT_MSG(false, "[LCC] Exception caught while reading input file.");
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
