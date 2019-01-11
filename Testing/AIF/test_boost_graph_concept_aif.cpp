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
#include "FEVV/Wrappings/Graph_traits_aif.h"
#include "FEVV/DataStructures/AIF/AIFMesh.hpp"

#include <boost/concept_check.hpp>

//------------------------------------------------------------------------------

template< class G >
struct GraphConcept
{
  typedef
      typename boost::graph_traits< G >::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits< G >::edge_descriptor edge_descriptor;
  typedef
      typename boost::graph_traits< G >::directed_category directed_category;
  typedef typename boost::graph_traits< G >::edge_parallel_category
      edge_parallel_category;
  typedef
      typename boost::graph_traits< G >::traversal_category traversal_category;

  void constraints()
  {
    BOOST_CONCEPT_ASSERT(
        (boost::DefaultConstructibleConcept< vertex_descriptor >));
    BOOST_CONCEPT_ASSERT(
        (boost::EqualityComparableConcept< vertex_descriptor >));
    BOOST_CONCEPT_ASSERT((boost::AssignableConcept< vertex_descriptor >));
    BOOST_CONCEPT_ASSERT(
        (boost::DefaultConstructibleConcept< edge_descriptor >));
    BOOST_CONCEPT_ASSERT((boost::EqualityComparableConcept< edge_descriptor >));
    BOOST_CONCEPT_ASSERT((boost::AssignableConcept< edge_descriptor >));
  }
  G g;
};

//------------------------------------------------------------------------------

int
main(int argc, const char **argv)
{
  struct GraphConcept< FEVV::DataStructures::AIF::AIFMesh > s;

  return 0;
}
