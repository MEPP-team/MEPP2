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

#include <algorithm>

// must be the last include to avoid side effect of disabling NDEBUG
#undef NDEBUG
#include <cassert>

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
               v2
              /|\
             / | \            v5
            /  |  \          /|
           /   |   \        / |
          /    |    \      /  |
       e0/     |e7   \e1  /   |
        /      |      \  /e9  |
       /       |       \/  f4 |
      /        |       /\     |e8
     /         |  f2  /  \    |
    /          |     /    \   |
   /     f3    |    /      \  |
  /            |   /        \ |
 /     e6      | /    e4     \|
v1 ---------- v0 ----------- v3
 \             |             /
  \            |            /
   \           |           /
    \          |          /
     \   f0    |  f1     /
      \        |        /
       \       |       /
      e3\      |e5    / e2
         \     |     /
          \    |    /
           \   |   /
            \  |  /
             \ | /
              \|/
               v4
  */

  smart_ptr_mesh m = AIFMesh::New();

  // add vertices
  vertex_descriptor v0 = helpers::add_vertex(m);
  vertex_descriptor v1 = helpers::add_vertex(m);
  vertex_descriptor v2 = helpers::add_vertex(m);
  vertex_descriptor v3 = helpers::add_vertex(m);
  vertex_descriptor v4 = helpers::add_vertex(m);
  vertex_descriptor v5 = helpers::add_vertex(m);

  // add edges
  edge_descriptor e0 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v1, e0, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v2, e0, helpers::vertex_pos::SECOND);
  edge_descriptor e1 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v2, e1, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v3, e1, helpers::vertex_pos::SECOND);
  edge_descriptor e2 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v3, e2, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v4, e2, helpers::vertex_pos::SECOND);
  edge_descriptor e3 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v1, e3, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v4, e3, helpers::vertex_pos::SECOND);
  edge_descriptor e4 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v0, e4, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v3, e4, helpers::vertex_pos::SECOND);
  edge_descriptor e5 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v4, e5, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v0, e5, helpers::vertex_pos::SECOND);
  edge_descriptor e6 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v1, e6, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v0, e6, helpers::vertex_pos::SECOND);
  edge_descriptor e7 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v0, e7, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v2, e7, helpers::vertex_pos::SECOND);
  edge_descriptor e8 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v3, e8, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v5, e8, helpers::vertex_pos::SECOND);
  edge_descriptor e9 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v0, e9, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v5, e9, helpers::vertex_pos::SECOND);

  // add faces [counter-clockwise oriented]
  face_descriptor f0 = helpers::add_face(m);
  helpers::link_edge_and_face(e5, f0);
  helpers::link_edge_and_face(e6, f0);
  helpers::link_edge_and_face(e3, f0);
  face_descriptor f1 = helpers::add_face(m);
  helpers::link_edge_and_face(e2, f1);
  helpers::link_edge_and_face(e4, f1);
  helpers::link_edge_and_face(e5, f1);
  face_descriptor f2 = helpers::add_face(m);
  helpers::link_edge_and_face(e4, f2);
  helpers::link_edge_and_face(e1, f2);
  helpers::link_edge_and_face(e7, f2);
  face_descriptor f3 = helpers::add_face(m);
  helpers::link_edge_and_face(e0, f3);
  helpers::link_edge_and_face(e6, f3);
  helpers::link_edge_and_face(e7, f3);
  face_descriptor f4 = helpers::add_face(m);
  helpers::link_edge_and_face(e4, f4);
  helpers::link_edge_and_face(e8, f4);
  helpers::link_edge_and_face(e9, f4);

  // display mesh informations
  std::cout << "Input mesh:\n";
  m->Print();

  //---------------------- Tests ---------------------

  // unordered one ring edges of vertex
  assert(helpers::is_e_one_ring_connected(
             e0, helpers::get_unordered_one_ring_edges(v0)) == true);
  assert(helpers::is_e_one_ring_connected(
             e1, helpers::get_unordered_one_ring_edges(v0)) == false);
  assert(helpers::is_e_one_ring_connected(
             e1, helpers::get_unordered_one_ring_edges(v0), true) == true);
  assert(helpers::is_e_one_ring_connected(
             e2, helpers::get_unordered_one_ring_edges(v0)) == false);
  assert(helpers::is_e_one_ring_connected(
             e2, helpers::get_unordered_one_ring_edges(v0), true) == true);
  assert(helpers::is_e_one_ring_connected(
             e3, helpers::get_unordered_one_ring_edges(v0)) == true);
  assert(helpers::is_e_one_ring_connected(
             e8, helpers::get_unordered_one_ring_edges(v0)) == false);
  assert(helpers::is_e_one_ring_connected(
             e8, helpers::get_unordered_one_ring_edges(v0), true) == false);
  std::cout << "testing  get_ordered_one_ring_edges(v0)" << std::endl;
  display_container(helpers::get_ordered_one_ring_edges(v0),
                    "helpers::get_ordered_one_ring_edges(v0)");
  display_container(std::vector< edge_descriptor >{e2, e1, e0, e3, e8},
                    "std::vector<edge_descriptor>{e2, e1, e0, e3, e8}");
  // auto c = helpers::get_ordered_one_ring_edges(v0);
  assert(rotationally_equal(
      helpers::get_ordered_one_ring_edges(v0),
      std::vector< edge_descriptor >{e1, e8, e2, e3, e0}) // ok
  );
  // ordered one ring of vertices of vertex

  std::cout << "testing  get_ordered_one_ring_vertices(v0)" << std::endl;
  assert(rotationally_equal(
      helpers::get_ordered_one_ring_vertices(v0),
      std::vector< vertex_descriptor >{v3, v5, v3, v4, v1, v2}) // ok
  );

  std::cout << "testing  get_ordered_one_ring_vertices(v1)" << std::endl;
  assert(rotationally_equal(helpers::get_ordered_one_ring_vertices(v1),
                            std::vector< vertex_descriptor >{v2, v0, v4}));

  std::cout << "testing  get_ordered_one_ring_vertices(v2)" << std::endl;
  assert(rotationally_equal(helpers::get_ordered_one_ring_vertices(v2),
                            std::vector< vertex_descriptor >{v3, v0, v1}));

  // ordered one ring of adjacent vertices of vertex

  std::cout << "testing  get_ordered_one_ring_of_adjacent_vertices(v0)"
            << std::endl;
  assert(rotationally_equal(
      helpers::get_ordered_one_ring_of_adjacent_vertices(v0),
      std::vector< vertex_descriptor >{v3, v5, v3, v4, v1, v2}) // ok
  );

  //---------------------- Tests even more complex ---------------------
  vertex_descriptor v6 = helpers::add_vertex(m);

  edge_descriptor e10 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v0, e10, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v6, e10, helpers::vertex_pos::SECOND);
  edge_descriptor e11 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v1, e11, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v6, e11, helpers::vertex_pos::SECOND);

  face_descriptor f5 = helpers::add_face(m);
  helpers::link_edge_and_face(e6, f5);
  helpers::link_edge_and_face(e11, f5);
  helpers::link_edge_and_face(e10, f5);
  // auto c = helpers::get_ordered_one_ring_edges(v0);
  assert(rotationally_equal(
      helpers::get_ordered_one_ring_edges(v0),
      std::vector< edge_descriptor >{e1, e8, e2, e3, e11, e0}) // ok
  );
  // ordered one ring of vertices of vertex

  std::cout << "testing  get_ordered_one_ring_vertices(v0)" << std::endl;
  assert(rotationally_equal(
      helpers::get_ordered_one_ring_vertices(v0),
      std::vector< vertex_descriptor >{v3, v5, v3, v4, v1, v6, v1, v2}) // ok
  );
 
  // bellow the code is a work in progress, please let it here as it is
#if 0 // dangling one-ring polylines not yet managed
  //---------------------- Tests even more complex ---------------------
  vertex_descriptor v7 = helpers::add_vertex(m);

  edge_descriptor e12 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v0, e12, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v7, e12, helpers::vertex_pos::SECOND);

  edge_descriptor e13 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v5, e13, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v7, e13, helpers::vertex_pos::SECOND);

  face_descriptor f6 = helpers::add_face(m);
  helpers::link_edge_and_face(e9, f6);
  helpers::link_edge_and_face(e13, f6);
  helpers::link_edge_and_face(e12, f6);
  auto c = helpers::get_ordered_one_ring_edges(v0);
  assert(rotationally_equal(helpers::get_ordered_one_ring_edges(v0),
std::vector<edge_descriptor>{e1, e8, e13, e2, e3, e11, e0}) // ok
    );
  // ordered one ring of vertices of vertex

  std::cout << "testing  get_ordered_one_ring_vertices(v0)" << std::endl;
  assert(rotationally_equal(helpers::get_ordered_one_ring_vertices(v0),
std::vector<vertex_descriptor>{v3, v5, v7, v5, v3, v4, v1, v6, v1, v2}) // ok
    );
#else
  edge_descriptor e12 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v5, e12, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v6, e12, helpers::vertex_pos::SECOND);

  face_descriptor f6 = helpers::add_face(m);
  helpers::link_edge_and_face(e9, f6);
  helpers::link_edge_and_face(e12, f6);
  helpers::link_edge_and_face(e10, f6);
  //auto c = helpers::get_ordered_one_ring_edges(v0);
  assert(rotationally_equal(helpers::get_ordered_one_ring_edges(v0), std::vector<edge_descriptor>{e11, e12, e3, e2, e8, e1, e0}) // ok
    ||
    rotationally_equal(helpers::get_ordered_one_ring_edges(v0), std::vector<edge_descriptor>{e12, e8, e2, e3, e11, e0, e1}) // ok
    ||
    rotationally_equal(helpers::get_ordered_one_ring_edges(v0), std::vector<edge_descriptor>{e12, e8, e1, e0, e11, e3, e2}) // ok
    ||
    rotationally_equal(helpers::get_ordered_one_ring_edges(v0), std::vector<edge_descriptor>{e0, e1, e8, e2, e3, e11, e12}) // ok
  );
  // ordered one ring of vertices of vertex

  //std::cout << "testing  get_ordered_one_ring_vertices(v0)" << std::endl;
  //assert(rotationally_equal(helpers::get_ordered_one_ring_vertices(v0), std::vector<vertex_descriptor>{v6, v5, v3, v2, v1, v4}) // not idea of expected result
  //  );
#endif
 
  return 0; // all tests passed
}
