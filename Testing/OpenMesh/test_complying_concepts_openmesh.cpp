// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <boost/graph/graph_concepts.hpp>
#define CGAL_USE_OM_POINTS
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>
#include "Testing/Concepts/CGALConceptsCheck.h"

typedef OpenMesh::PolyMesh_ArrayKernelT< /*MyTraits*/ > Mesh;
// Technical type alias in order for the tests for the different data
// structures to share some code (preparing generic test for the day
// when all data structure are fully compliant).
typedef Mesh G;

int
main(int narg, char **argv)
{
  using namespace boost;
  // Assertions are sometimes redundant (IncidenceGraph concept assertion
  // already internally assert the Graph concept). We nevertheless do such
  // redudant checks along the concept "inheritance" hieararchy in order
  // for the test to provide the simplest feedback when failing (fail early).
  BOOST_CONCEPT_ASSERT((GraphConcept< G >));
  BOOST_CONCEPT_ASSERT((IncidenceGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((BidirectionalGraphConcept< G >));

  // Graph already asserted
  BOOST_CONCEPT_ASSERT((EdgeListGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((VertexListGraphConcept< G >));
  BOOST_CONCEPT_ASSERT((VertexAndEdgeListGraphConcept< G >));

  //////////////////////////////////////////////////////////////////////////
  ////////////// BOOST_CONCEPT_ASSERT( (EdgeMutableGraphConcept< G >) );
  // The following is the unfolding of the FEVV::CGALHalfedgeGraph concept
  // on debugging (of Polyhedron) purposes
  { // Technical scope to avoid variable/type collisions with other tests
    typedef typename graph_traits< G >::edge_descriptor edge_descriptor;
    G g;
    edge_descriptor e;
    // std::pair<edge_descriptor, bool> p;
    // typename graph_traits<G>::vertex_descriptor u, v;

    // The following assertion fails with message:
    //   error: no matching function for call to 'add_edge'
    // UNCOMMENT ME (FIXME) p = add_edge(u, v, g);

    // The following assertion fails because the returned type is not
    // convertible to an edge_descriptor. UNCOMMENT ME (FIXME) remove_edge(u, v,
    // g);

    remove_edge(e, g);

    // The following assertion fails with message:
    //   error: use of undeclared identifier 'clear_vertex(v, g);`
    // UNCOMMENT ME (FIXME) clear_vertex(v, g);
  } // Technical scope

  BOOST_CONCEPT_ASSERT((VertexMutableGraphConcept< G >));

  //////////////////////////////////////////////////////////////////////////
  //////////// BOOST_CONCEPT_ASSERT( (MutableGraphConcept<      G>) );
  // MutableGraph concept doesn't introduce new expressions. It only fails
  // because the EdgeMutableGraph concept from which MutableGraph inherits
  // fails.

  //////////////////////////////////////////////////////////////////////////
  // BOOST_CONCEPT_ASSERT(
  //   (boost::concepts::ReadablePropertyGraph< G,
  //                          graph_traits< G >::vertex_descriptor,
  //                          vertex_point_t >)
  // );
  // The following is the unfolding of the above failing concept assertion of
  // ReadablePropertyGraph on debugging purposes:
  { // Technical scope to avoid variable/type collisions with other tests
    typedef graph_traits< G >::vertex_descriptor X;
    typedef vertex_point_t Property;

    typedef typename property_map< G, Property >::const_type const_Map;
    typename property_traits< const_Map >::value_type pval;
    G g;
    X x;
    const_Map pmap = get(Property(), g);
    // The following assertion fails with message:
    //   error: no viable overloaded '='
    // UNCOMMENT ME (FIXME) pval = get(Property(), g, x);
    ignore_unused_variable_warning(pmap);
  }
  // UNCOMMENT the following FIXMEs once the above ReadablePropertyGraph
  // succeeds.
  // FIXME  BOOST_CONCEPT_ASSERT(
  // FIXME    (PropertyGraphConcept< G,
  // FIXME                           graph_traits< G >::vertex_descriptor,
  // FIXME                           vertex_point_t >)
  // FIXME  );

  BOOST_CONCEPT_ASSERT((FEVV::CGALHalfedgeGraph< G >));
  BOOST_CONCEPT_ASSERT((FEVV::CGALHalfedgeListGraph< G >));
  BOOST_CONCEPT_ASSERT((FEVV::CGALMutableHalfedgeGraph< G >));

  BOOST_CONCEPT_ASSERT((FEVV::CGALFaceGraph< G >));
  BOOST_CONCEPT_ASSERT((FEVV::CGALFaceListGraph< G >));
  BOOST_CONCEPT_ASSERT((FEVV::CGALMutableFaceGraph< G >));

  // Just to test the implementation of the FEVV::MutableHalfedgeFaceListGraph
  // concept as opposed to the concept itself. Indeed this concept is
  // just a syntactic gathering of above already asserted concepts (and
  // as such should never fails provided the above succeeds).
  BOOST_CONCEPT_ASSERT((FEVV::MutableHalfedgeFaceListGraph< G >));

  //////////////////////////////////////////////////////////////////////////
  // The following are not mandatory within the FEVV layer. It is just to
  // to document some technical perks or failures of AIF.
  // BOOST_CONCEPT_ASSERT( (concepts::AdjacencyGraph<  G >) );

  return 0;
}
