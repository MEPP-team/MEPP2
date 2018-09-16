#include <CGAL/Surface_mesh.h>
#include <CGAL/Cartesian.h>
#include "Testing/Concepts/GeometryConceptCheck.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"

int
main(int narg, char **argv)
{
  typedef CGAL::Cartesian< double > Kernel;
  typedef CGAL::Surface_mesh< Kernel::Point_3 > Mesh;

  BOOST_CONCEPT_ASSERT((FEVV::Point_3Concept< Mesh >));

  typedef FEVV::Geometry_traits< Mesh > Geometry;
  Mesh mesh;
  Geometry g(mesh);
  // We don't have to call GeometryTraits::verify_assumption() because
  // we known that with a Kernel::Point_3 kernel CGAL::Surface_mesh
  // hardwires an ambiant space dimension of 3:
  // FEVV::GeometryTraits::verify_assumption( g );
  FEVV::GeometryTraits::verify_operator_results(g);

  return 0;
}
