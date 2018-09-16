#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include "Testing/Concepts/GeometryConceptCheck.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"

int
main(int narg, char **argv)
{
  typedef CGAL::Cartesian< double > Kernel;
  typedef CGAL::Polyhedron_3< Kernel, CGAL::Polyhedron_items_with_id_3 > Mesh;

  BOOST_CONCEPT_ASSERT((FEVV::Point_3Concept< Mesh >));

  typedef FEVV::Geometry_traits< Mesh > Geometry;
  Mesh mesh;
  Geometry g(mesh);
  // We don't have to call GeometryTraits::verify_assumption() because
  // we known that a CGAL::Polyhedron_3 hardwires an ambiant space
  // dimension of 3:
  // FEVV::GeometryTraits::verify_assumption( g );
  FEVV::GeometryTraits::verify_operator_results(g);

  return 0;
}
