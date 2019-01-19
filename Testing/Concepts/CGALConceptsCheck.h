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
#pragma once

//#include <boost/concept_check.hpp>
//#include <boost/concept/assert.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/concept/detail/concept_def.hpp>

namespace FEVV {

using namespace boost;
/**
 * \brief [CGAL's
 * HalfedgeGraph](http://doc.cgal.org/latest/BGL/classHalfedgeGraph.html)
 *        concept (expressed with
 *        [Boost Concept Check
 * Library](http://www.boost.org/doc/libs/1_60_0/libs/concept_check/concept_check.htm))
 * (BCCL)
 */
BOOST_concept(CGALHalfedgeGraph, (G))
{
  typedef typename boost::graph_traits< G >::halfedge_descriptor
      halfedge_descriptor;
  typedef typename boost::graph_traits< G >::edge_descriptor edge_descriptor;
  typedef
      typename boost::graph_traits< G >::vertex_descriptor vertex_descriptor;
  BOOST_CONCEPT_USAGE(CGALHalfedgeGraph)
  {
    BOOST_CONCEPT_ASSERT((DefaultConstructible< halfedge_descriptor >));
    BOOST_CONCEPT_ASSERT((Assignable< halfedge_descriptor >));
    BOOST_CONCEPT_ASSERT((EqualityComparable< halfedge_descriptor >));
    BOOST_CONCEPT_ASSERT((LessThanComparable< halfedge_descriptor >));
    e = edge(h, g);
    h = halfedge(e, g);
    h = halfedge(v, g);
    p = halfedge(u, v, g);
    h = opposite(h, g);
    v = source(h, g);
    v = target(h, g);
    h = next(h, g);
    h = prev(h, g);
    h = boost::graph_traits< G >::null_halfedge();
  }
  halfedge_descriptor h;
  edge_descriptor e;
  vertex_descriptor u, v;
  std::pair< halfedge_descriptor, bool > p;
  G g;
};

/**
 * \brief [CGAL's
 * HalfedgeListGraph](http://doc.cgal.org/latest/BGL/classHalfedgeListGraph.html)
 *        concept (expressed with
 *        [Boost Concept Check
 * Library](http://www.boost.org/doc/libs/1_60_0/libs/concept_check/concept_check.htm))
 * (BCCL)
 */
BOOST_concept(CGALHalfedgeListGraph, (G)) : CGALHalfedgeGraph< G >
{
  typedef
      typename boost::graph_traits< G >::halfedge_iterator halfedge_iterator;
  typedef typename boost::graph_traits< G >::halfedges_size_type
      halfedges_size_type;
  BOOST_CONCEPT_USAGE(CGALHalfedgeListGraph)
  {
    BOOST_CONCEPT_ASSERT((DefaultConstructible< halfedge_iterator >));
    BOOST_CONCEPT_ASSERT((Assignable< halfedge_iterator >));
    BOOST_CONCEPT_ASSERT((EqualityComparable< halfedge_iterator >));
    E = num_halfedges(g);
    p = halfedges(g);
  }
  std::pair< halfedge_iterator, halfedge_iterator > p;
  halfedges_size_type E;
  G g;
};

/**
 * \brief [CGAL's
 * MutableHalfedgeGraph](http://doc.cgal.org/latest/BGL/classMutableHalfedgeGraph.html)
 *        concept (expressed with
 *        [Boost Concept Check
 * Library](http://www.boost.org/doc/libs/1_60_0/libs/concept_check/concept_check.htm))
 * (BCCL)
 */
BOOST_concept(CGALMutableHalfedgeGraph, (G))
{
  typedef typename boost::graph_traits< G >::halfedge_descriptor
      halfedge_descriptor;
  typedef typename boost::graph_traits< G >::edge_descriptor edge_descriptor;
  typedef
      typename boost::graph_traits< G >::vertex_descriptor vertex_descriptor;
  BOOST_CONCEPT_USAGE(CGALMutableHalfedgeGraph)
  {
    v = add_vertex(g);
    remove_vertex(v, g);
    e = add_edge(g);
    remove_edge(e, g);
    set_target(h, v, g);
    set_halfedge(v, h, g);
    set_next(h1, h2, g);
  }
  halfedge_descriptor h, h1, h2;
  edge_descriptor e;
  vertex_descriptor v;
  G g;
};

/**
 * \brief [CGAL's FaceGraph](http://doc.cgal.org/latest/BGL/classFaceGraph.html)
 *        concept (expressed with
 *        [Boost Concept Check
 * Library](http://www.boost.org/doc/libs/1_60_0/libs/concept_check/concept_check.htm))
 * (BCCL)
 */
BOOST_concept(CGALFaceGraph, (G)) : CGALHalfedgeGraph< G >
{
  typedef typename boost::graph_traits< G >::face_descriptor face_descriptor;
  typedef typename boost::graph_traits< G >::halfedge_descriptor
      halfedge_descriptor;
  BOOST_CONCEPT_USAGE(CGALFaceGraph)
  {
    f = face(h, g);
    h = halfedge(f, g);
    (void)degree(f, g);
    f = boost::graph_traits< G >::null_face();
  }
  face_descriptor f;
  halfedge_descriptor h;
  G g;
};

/**
 * \brief [CGAL's
 * FaceListGraph](http://doc.cgal.org/latest/BGL/classFaceListGraph.html)
 *        concept (expressed with
 *        [Boost Concept Check
 * Library](http://www.boost.org/doc/libs/1_60_0/libs/concept_check/concept_check.htm))
 * (BCCL)
 */
BOOST_concept(CGALFaceListGraph, (G)) : CGALFaceGraph< G >
{
  typedef typename boost::graph_traits< G >::face_iterator face_iterator;
  typedef typename boost::graph_traits< G >::faces_size_type faces_size_type;

  BOOST_CONCEPT_USAGE(CGALFaceListGraph)
  {
    p = faces(g);
    F = num_faces(g);
  }
  std::pair< face_iterator, face_iterator > p;
  faces_size_type F;
  G g;
};

/**
 * \brief [CGAL's
 * MutableFaceGraph](http://doc.cgal.org/latest/BGL/classMutableFaceGraph.html)
 *        concept (expressed with
 *        [Boost Concept Check
 * Library](http://www.boost.org/doc/libs/1_60_0/libs/concept_check/concept_check.htm))
 * (BCCL)
 */
BOOST_concept(CGALMutableFaceGraph, (G))
    : CGALFaceGraph< G >, CGALMutableHalfedgeGraph< G >
{
  typedef typename boost::graph_traits< G >::face_descriptor face_descriptor;
  typedef typename boost::graph_traits< G >::halfedge_descriptor
      halfedge_descriptor;
  BOOST_CONCEPT_USAGE(CGALMutableFaceGraph)
  {
    f = add_face(g);
    remove_face(f, g);
    set_face(h, f, g);
    set_halfedge(f, h, g);
  }
  face_descriptor f;
  halfedge_descriptor h;
  G g;
};

/**
 * \brief A notational convenience concept gathering CGAL's MutableFaceGraph,
 *        HaledgeListGraph and FaceListGraph concepts (expressed with
 *        [Boost Concept Check
 * Library](http://www.boost.org/doc/libs/1_60_0/libs/concept_check/concept_check.htm))
 * (BCCL)
 */
BOOST_concept(MutableHalfedgeFaceListGraph, (G))
    : CGALMutableFaceGraph< G >, CGALHalfedgeListGraph< G >,
      CGALFaceListGraph< G >{};

} // namespace FEVV

