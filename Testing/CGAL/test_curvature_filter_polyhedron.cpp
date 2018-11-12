#include "FEVV/DataStructures/DataStructures_cgal_polyhedron_3.h"

#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/properties_polyhedron_3.h"

#include "Testing/Generic/Manifold/test_curvature_filter.inl"


int
main(int argc, const char **argv)
{
  return test_curvature_filter< FEVV::MeshPolyhedron >(argc, argv);
}
