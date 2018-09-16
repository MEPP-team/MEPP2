#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/Wrappings/Geometry_traits_aif.h"

#include "Testing/Concepts/GeometryConceptCheck.h"

int
main(int narg, char **argv)
{
  typedef FEVV::DataStructures::AIF::AIFMesh Mesh;
  BOOST_CONCEPT_ASSERT((FEVV::Point_3Concept< Mesh >));

  typedef FEVV::Geometry_traits< Mesh > Geometry;
  Mesh mesh;
  Geometry g(mesh);
  FEVV::GeometryTraits::verify_assumption(g);
  // FEVV::GeometryTraits::verify_operator_results( g );

  return 0;
}
