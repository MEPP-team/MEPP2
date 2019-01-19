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

#include "FEVV/Wrappings/Geometry_traits_aif.h"       // for k_ring testing
#include "FEVV/Operators/Generic/Manifold/k_ring.hpp" // for k_ring testing

#include <algorithm>

// must be the last include to avoid side effect of disabling NDEBUG
#undef NDEBUG
#include <cassert>

using namespace FEVV; // for k_ring testing

//------------------------------------------------------------------------------

/*
 * Display container or range content.
 */
template< typename T >
void
display_container(const T &container, const char *title)
{
  std::cout << title << " =" << std::endl;
  for(auto element : container)
    std::cout << "  " << element << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

const bool VERBOSE = false;

/*
 *  Returns true if 2 containers are equal by rotation.
 */
template< typename T >
bool
rotationally_equal(const T &v1, T v2, bool test_in_reverse_order_too = true)
{
  bool are_equal = false;

  while(true)
  {
    for(int i = 0; i < v1.size(); i++)
    {
      if(VERBOSE)
      {
        std::cout << "rotate number = " << i << std::endl;
        display_container(v1, "v1");
        display_container(v2, "v2");
      }

      if(v1 == v2)
      {
        if(VERBOSE)
          std::cout << "v1 == v2" << std::endl;
        are_equal = true;
        break;
      }
      else
      {
        if(VERBOSE)
          std::cout << "v1 != v2" << std::endl;
      }

      // rotate v2 and test again
      std::rotate(v2.begin(), v2.begin() + 1, v2.end());
    }

    if(are_equal == false && test_in_reverse_order_too)
    {
      if(VERBOSE)
        std::cout << "reverse v2" << std::endl;
      std::reverse(v2.begin(), v2.end());
      test_in_reverse_order_too = false;
    }
    else
      break;
  }

  return are_equal;
}

//------------------------------------------------------------------------------

int
main(int argc, const char **argv)
{
  typedef FEVV::DataStructures::AIF::AIFMesh AIFMesh;
  typedef FEVV::DataStructures::AIF::AIFVertex AIFVertex;
  typedef FEVV::DataStructures::AIF::AIFEdge AIFEdge;
  typedef FEVV::DataStructures::AIF::AIFFace AIFFace;
  typedef FEVV::DataStructures::AIF::AIFTopologyHelpers helpers;
  typedef helpers::AIFHalfEdge AIFHalfEdge;
  typedef helpers::smart_ptr_mesh smart_ptr_mesh;
  typedef helpers::vertex_descriptor vertex_descriptor;
  typedef helpers::edge_descriptor edge_descriptor;
  typedef helpers::halfedge_descriptor halfedge_descriptor;
  typedef helpers::face_descriptor face_descriptor;


  // build mesh
  /*
               v2 ----------- v6
              /|\     e14      \
             / | \              \
            /  |  \              \
           /   |   \              \
          /    |    \              \e0
       e4/     |e15  \e5            \
        /      |      \              \
       /       |       \              \
      /        |        \              \           v9
     /         |         \              \         / \
    /          |          \              \       /   \
   /     f4    |    f3     \      f2      \   e3/     \e1
  /            |            \              \   /   f5  \
 /     e13     |      e9     \              \ /         \
v1 ---------- v0 ----------- v3             v7 --------- v8
 \             |              |             /     e7
  \            |              |            /
   \           |              |           /
    \          |              |          /
     \   f0    |      f1      |         /
      \        |              |        /
       \       |              |       /            v11 (isolated)
     e12\      |e11         e6|    e2/            X
         \     |              |     /
          \    |              |    /
           \   |              |   /
            \  |              |  /
             \ |              | /
              \|      e10     |/      e8
               v5 ----------- v4 ----------- v10
  */

  smart_ptr_mesh m = AIFMesh::New();

  // add vertices
  vertex_descriptor v0 = helpers::add_vertex(m);
  vertex_descriptor v1 = helpers::add_vertex(m);
  vertex_descriptor v2 = helpers::add_vertex(m);
  vertex_descriptor v3 = helpers::add_vertex(m);
  vertex_descriptor v4 = helpers::add_vertex(m);
  vertex_descriptor v5 = helpers::add_vertex(m);
  vertex_descriptor v6 = helpers::add_vertex(m);
  vertex_descriptor v7 = helpers::add_vertex(m);
  vertex_descriptor v8 = helpers::add_vertex(m);
  vertex_descriptor v9 = helpers::add_vertex(m);
  vertex_descriptor v10 = helpers::add_vertex(m);
  vertex_descriptor v11 = helpers::add_vertex(m);

  // add edges
  edge_descriptor e0 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v6, e0, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v7, e0, helpers::vertex_pos::SECOND);
  edge_descriptor e1 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v8, e1, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v9, e1, helpers::vertex_pos::SECOND);
  edge_descriptor e2 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v4, e2, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v7, e2, helpers::vertex_pos::SECOND);
  edge_descriptor e3 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v7, e3, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v9, e3, helpers::vertex_pos::SECOND);
  edge_descriptor e4 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v1, e4, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v2, e4, helpers::vertex_pos::SECOND);
  edge_descriptor e5 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v2, e5, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v3, e5, helpers::vertex_pos::SECOND);
  edge_descriptor e6 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v3, e6, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v4, e6, helpers::vertex_pos::SECOND);
  edge_descriptor e7 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v7, e7, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v8, e7, helpers::vertex_pos::SECOND);
  edge_descriptor e8 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v4, e8, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v10, e8, helpers::vertex_pos::SECOND);
  edge_descriptor e9 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v0, e9, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v3, e9, helpers::vertex_pos::SECOND);
  edge_descriptor e10 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v4, e10, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v5, e10, helpers::vertex_pos::SECOND);
  edge_descriptor e11 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v5, e11, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v0, e11, helpers::vertex_pos::SECOND);
  edge_descriptor e12 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v1, e12, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v5, e12, helpers::vertex_pos::SECOND);
  edge_descriptor e13 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v1, e13, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v0, e13, helpers::vertex_pos::SECOND);
  edge_descriptor e14 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v6, e14, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v2, e14, helpers::vertex_pos::SECOND);
  edge_descriptor e15 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v0, e15, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v2, e15, helpers::vertex_pos::SECOND);

  // add faces [counter-clockwise oriented]
  face_descriptor f0 = helpers::add_face(m);
  helpers::link_edge_and_face(e11, f0);
  helpers::link_edge_and_face(e13, f0);
  helpers::link_edge_and_face(e12, f0);
  face_descriptor f1 = helpers::add_face(m);
  helpers::link_edge_and_face(e6, f1);
  helpers::link_edge_and_face(e9, f1);
  helpers::link_edge_and_face(e11, f1);
  helpers::link_edge_and_face(e10, f1);
  face_descriptor f2 = helpers::add_face(m);
  helpers::link_edge_and_face(e2, f2);
  helpers::link_edge_and_face(e0, f2);
  helpers::link_edge_and_face(e14, f2);
  helpers::link_edge_and_face(e5, f2);
  helpers::link_edge_and_face(e6, f2);
  face_descriptor f3 = helpers::add_face(m);
  helpers::link_edge_and_face(e9, f3);
  helpers::link_edge_and_face(e5, f3);
  helpers::link_edge_and_face(e15, f3);
  face_descriptor f4 = helpers::add_face(m);
  helpers::link_edge_and_face(e4, f4);
  helpers::link_edge_and_face(e13, f4);
  helpers::link_edge_and_face(e15, f4);
  face_descriptor f5 = helpers::add_face(m);
  helpers::link_edge_and_face(e3, f5);
  helpers::link_edge_and_face(e7, f5);
  helpers::link_edge_and_face(e1, f5);

  // display mesh informations
  std::cout << "Input mesh:\n";
  m->Print();

  //---------------------- Tests ---------------------

  // unordered one ring edges of vertex

  std::cout << "testing  get_unordered_one_ring_edges(v0)" << std::endl;
  // display_container(helpers::get_unordered_one_ring_edges(v0),
  // "get_unordered_one_ring_edges(v0)");
  // display_container(std::set<edge_descriptor>{e4, e5, e12, e6, e10},
  // "std::set<edge_descriptor>{e4, e5, e12, e6, e10}");
  assert((helpers::get_unordered_one_ring_edges(v0) ==
          std::set< edge_descriptor >{e4, e5, e12, e6, e10}));

  std::cout << "testing  get_unordered_one_ring_edges(v1)" << std::endl;
  assert((helpers::get_unordered_one_ring_edges(v1) ==
          std::set< edge_descriptor >{e11, e15}));

  std::cout << "testing  get_unordered_one_ring_edges(v2)" << std::endl;
  assert((helpers::get_unordered_one_ring_edges(v2) ==
          std::set< edge_descriptor >{e9, e13, e6, e2, e0}));

  std::cout << "testing  get_unordered_one_ring_edges(v11)" << std::endl;
  assert((helpers::get_unordered_one_ring_edges(v11) ==
          std::set< edge_descriptor >{}));


  // ordered one ring edges of vertex

  std::cout << "testing  get_ordered_one_ring_edges(v0)" << std::endl;
  assert(
      rotationally_equal(helpers::get_ordered_one_ring_edges(v0),
                         std::vector< edge_descriptor >{e5, e6, e10, e12, e4}));

  // display_container(helpers::get_ordered_one_ring_edges(v1),
  // "get_ordered_one_ring_edges(v1)");
  std::cout << "testing  get_ordered_one_ring_edges(v1)" << std::endl;
  assert(rotationally_equal(helpers::get_ordered_one_ring_edges(v1),
                            std::vector< edge_descriptor >{e15, e11}));

  // display_container(helpers::get_ordered_one_ring_edges(v2),
  // "get_ordered_one_ring_edges(v2)");
  std::cout << "testing  get_ordered_one_ring_edges(v2)" << std::endl;
  assert(
      rotationally_equal(helpers::get_ordered_one_ring_edges(v2),
                         std::vector< edge_descriptor >{e0, e2, e6, e9, e13}));

  std::cout << "testing  get_ordered_one_ring_edges(v11)" << std::endl;
  assert((helpers::get_ordered_one_ring_edges(v11) ==
          std::vector< edge_descriptor >{}));


  // ordered one ring of vertices of vertex

  std::cout << "testing  get_ordered_one_ring_vertices(v0)" << std::endl;
  assert(
      rotationally_equal(helpers::get_ordered_one_ring_vertices(v0),
                         std::vector< vertex_descriptor >{v2, v3, v4, v5, v1}));

  std::cout << "testing  get_ordered_one_ring_vertices(v1)" << std::endl;
  assert(rotationally_equal(helpers::get_ordered_one_ring_vertices(v1),
                            std::vector< vertex_descriptor >{v2, v0, v5}));

  std::cout << "testing  get_ordered_one_ring_vertices(v2)" << std::endl;
  assert(rotationally_equal(
      helpers::get_ordered_one_ring_vertices(v2),
      std::vector< vertex_descriptor >{v6, v7, v4, v3, v0, v1}));

  std::cout << "testing  get_ordered_one_ring_vertices(v11)" << std::endl;
  assert((helpers::get_ordered_one_ring_vertices(v11) ==
          std::vector< vertex_descriptor >{}));


  // ordered one ring of adjacent vertices of vertex

  std::cout << "testing  get_ordered_one_ring_of_adjacent_vertices(v0)"
            << std::endl;
  assert(
      rotationally_equal(helpers::get_ordered_one_ring_of_adjacent_vertices(v0),
                         std::vector< vertex_descriptor >{v1, v2, v3, v5}));

  std::cout << "testing  get_ordered_one_ring_of_adjacent_vertices(v1)"
            << std::endl;
  assert(
      rotationally_equal(helpers::get_ordered_one_ring_of_adjacent_vertices(v1),
                         std::vector< vertex_descriptor >{v2, v0, v5}));

  std::cout << "testing  get_ordered_one_ring_of_adjacent_vertices(v2)"
            << std::endl;
  assert(
      rotationally_equal(helpers::get_ordered_one_ring_of_adjacent_vertices(v2),
                         std::vector< vertex_descriptor >{v6, v3, v0, v1}));

  std::cout << "testing  get_ordered_one_ring_of_adjacent_vertices(v11)"
            << std::endl;
  assert((helpers::get_ordered_one_ring_of_adjacent_vertices(v11) ==
          std::vector< vertex_descriptor >{}));

  // test one-ring-vertices cache

  std::cout << "testing  one-ring-vertices cache" << std::endl;
  auto one_r1 = helpers::get_ordered_one_ring_vertices(v2);
  helpers::remove_face(f4, m);
  auto one_r2 = helpers::get_ordered_one_ring_vertices(v2);
  assert(one_r1 != one_r2);

  // test extract_k_ring [k_ring.h]

  std::cout << "testing  extract_k_ring" << std::endl;
  std::vector< vertex_descriptor > qv;
  Operators::extract_k_ring(v0, *m, 0, qv);
  assert((qv == std::vector< vertex_descriptor >{v0}));
  qv.clear();
  Operators::extract_k_ring(
      v0,
      *m,
      1,
      qv); // be careful face f4 has been removed before this current test
  assert((qv == std::vector< vertex_descriptor >{v0, v3, v5, v1, v2}));
  qv.clear();
  Operators::extract_k_ring(v0, *m, 2, qv);
  assert((qv == std::vector< vertex_descriptor >{v0, v3, v5, v1, v2, v4, v6}));
  helpers::remove_edge(e8, m); // remove e8 to ensure 2-manifoldness (else
                               // extract_k_ring will throw an exception)
  qv.clear();
  Operators::extract_k_ring(v0, *m, 3, qv);
  assert(
      (qv == std::vector< vertex_descriptor >{v0, v3, v5, v1, v2, v4, v6, v7}));
  // cannot go further away, because v7 is not a manifold vertex (cut-vertex)

  qv.clear();
  Operators::extract_k_ring(v8, *m, 1, qv);
  assert((qv == std::vector< vertex_descriptor >{v8, v9, v7}));
  qv.clear();
  Operators::extract_k_ring(v9, *m, 1, qv);
  assert((qv == std::vector< vertex_descriptor >{v9, v8, v7}));

  // test extract_vertex_star [k_ring.h]

  std::cout << "testing  extract_vertex_star" << std::endl;
  std::vector< halfedge_descriptor > qh;
  Operators::extract_vertex_star(v0, *m, qh);
  assert((qh == std::vector< halfedge_descriptor >{
                    halfedge_descriptor(v3, v0, f1, e9, false),
                    halfedge_descriptor(v5, v0, f0, e11, false),
                    halfedge_descriptor(
                        v0, v1, f0, e13, false), // f4 has been removed
                    halfedge_descriptor(v2, v0, f3, e15, false)}));
  qh.clear();

  Operators::extract_vertex_star(v6, *m, qh);
  assert((qh == std::vector< halfedge_descriptor >{
                    halfedge_descriptor(v7, v6, f2, e0, false),
                    halfedge_descriptor(v6, v2, f2, e14, false)}));
  qh.clear();

  std::cout << "Test passed.\n";
  return 0; // all tests passed
}
