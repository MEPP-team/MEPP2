#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/Wrappings/Graph_traits_aif.h"
#include "FEVV/Wrappings/Graph_properties_aif.h" // to validate the propertygraph concept

#include <boost/graph/graph_concepts.hpp>
#include "Testing/Concepts/CGALConceptsCheck.h"

typedef FEVV::DataStructures::AIF::AIFMesh Mesh;
// Technical type alias in order for the tests for the different data
// structures to share some code (preparing generic test for the day
// when all data structure are fully compliant).
typedef Mesh G;

int
main(int narg, char **argv)
{
  using namespace boost;

  //////////////////////////////////////////////////////////////////////////
  // Graph already asserted

  // Assertions are sometimes redundant (IncidenceGraph concept assertion
  // already internally assert the Graph concept). We nevertheless do such
  // redudant checks along the concept "inheritance" hieararchy in order
  // for the test to provide the simplest feedback when failing (fail early).

  // boost Graph concepts
  BOOST_CONCEPT_ASSERT((GraphConcept< G >));
  BOOST_CONCEPT_ASSERT((IncidenceGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((BidirectionalGraphConcept< G >));

  BOOST_CONCEPT_ASSERT((VertexListGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((EdgeListGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((VertexAndEdgeListGraphConcept< G >));

  BOOST_CONCEPT_ASSERT((EdgeMutableGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((VertexMutableGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((MutableGraphConcept< G >));

  BOOST_CONCEPT_ASSERT(
      (PropertyGraphConcept< G,
                             typename graph_traits< G >::vertex_descriptor,
                             boost::vertex_point_t >));


  // CGAL concepts
  BOOST_CONCEPT_ASSERT((FEVV::CGALHalfedgeGraph< G >));
  BOOST_CONCEPT_ASSERT((FEVV::CGALMutableHalfedgeGraph< G >));
  BOOST_CONCEPT_ASSERT((FEVV::CGALFaceGraph< G >));
  BOOST_CONCEPT_ASSERT((FEVV::CGALFaceListGraph< G >));
  BOOST_CONCEPT_ASSERT((FEVV::CGALMutableFaceGraph< G >));
  //////////////////////////////////////////////////////////////////////////
  // BOOST_CONCEPT_ASSERT((FEVV::CGALHalfedgeListGraph<    Mesh >));
  // This concept is not compatible with the design of AIF datastructure.
  // It thus should not be asserted.

  //////////////////////////////////////////////////////////////////////////
  // BOOST_CONCEPT_ASSERT((FEVV::MutableHalfedgeFaceListGraph< Mesh >));
  // This concept is not compatible with the design of AIF datastructure.
  // It thus should not be asserted.

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  // The following are not mandatory within the FEVV layer. It is just to
  // to document some technical perks or failures of AIF.

  // boost Graph concepts
  BOOST_CONCEPT_ASSERT((AdjacencyGraphConcept< G >));

  /// \todo: should ReadablePropertyGraphConcept on boost::vertex_index_t be
  /// mandatory for all FEVV datastructures? Test it.
  BOOST_CONCEPT_ASSERT((ReadablePropertyGraphConcept<
                        G,
                        typename graph_traits< G >::vertex_descriptor,
                        boost::vertex_index_t >));

  /// \todo: test it for all FEVV datastructures to see if mandatory or not
  BOOST_CONCEPT_ASSERT((
      LvaluePropertyGraphConcept< G,
                                  typename graph_traits< G >::vertex_descriptor,
                                  boost::vertex_point_t >));

  // The following cannot work (is expected to fail) since boost::vertex_index_t
  // should not be modifiable (as opposed to vertex positions FAIL
  // BOOST_CONCEPT_ASSERT(
  //        (WritablePropertyMapConcept<
  //           typename property_map< G, boost::vertex_index_t >::type,
  //           typename graph_traits< G>::vertex_descriptor>));


#ifndef WIN32
  // BOOST_CONCEPT_ASSERT(
  //   (PropertyGraphConcept< G,
  //                          graph_traits< G >::halfedge_descriptor,
  //                          halfedge_index_t >)
  // );
#endif

  BOOST_CONCEPT_ASSERT((VertexIndexGraphConcept< G >));
  // FAIL BOOST_CONCEPT_ASSERT((EdgeIndexGraphConcept< G>)); // usefull for AIF:
  // yes concept implemented in a future release

  BOOST_CONCEPT_ASSERT((MutableIncidenceGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((MutableBidirectionalGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((MutableEdgeListGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((VertexMutablePropertyGraphConcept< G >));
  // FAIL BOOST_CONCEPT_ASSERT((EdgeMutablePropertyGraphConcept< G>)); // do not
  // attach any property to edge for the time being
  BOOST_CONCEPT_ASSERT((AdjacencyMatrixConcept< G >));

  return 0;
}
