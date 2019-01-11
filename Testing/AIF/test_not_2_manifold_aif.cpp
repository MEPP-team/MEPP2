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
#include "FEVV/DataStructures/AIF/AIFMeshHelpers.h"

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

//------------------------------------------------------------------------------

int
main(int argc, const char **argv)
{
  typedef FEVV::DataStructures::AIF::AIFMesh AIFMesh;
  typedef FEVV::DataStructures::AIF::AIFVertex AIFVertex;
  typedef FEVV::DataStructures::AIF::AIFEdge AIFEdge;
  typedef FEVV::DataStructures::AIF::AIFFace AIFFace;
  typedef FEVV::DataStructures::AIF::AIFTopologyHelpers TopoHelpers;
  typedef TopoHelpers::AIFHalfEdge AIFHalfEdge;
  typedef TopoHelpers::smart_ptr_mesh smart_ptr_mesh;
  typedef TopoHelpers::vertex_descriptor vertex_descriptor;
  typedef TopoHelpers::edge_descriptor edge_descriptor;
  typedef TopoHelpers::face_descriptor face_descriptor;


  // build mesh
  /*                 e14 (dangling)
               v2 ----------- v6
              /|
             / |
            /  |
           /   | e15
          /    |
       e4/     |  v11
        /      |  /\
       /   f4  |/   \
      /      / |      \
     /     /   |       \ e16
    /    /e17  |        \
   /   /       |    f5   \
  /  /         |          \
 / /   e13     |      e9   \
v1 ---------- v0 ----------- v3             v7 --------- v8
 \             |\           /                      e7 (isolated)
  \            |  \   f6   /
   \           |e19\      / e18
    \          |    \    /
     \   f0    |     \  /
      \        |      \/
       \       |      v12               v9 (isolated)
     e12\      |e11                          X
         \     |
          \    |
           \   |
            \  |
             \ |
              \|     e10 (dangling)   e8 (dangling)
               v5 ----------- v4 ----------- v10
  */

  smart_ptr_mesh m = AIFMesh::New();

  // add vertices
  vertex_descriptor v0 = TopoHelpers::add_vertex(m);
  vertex_descriptor v1 = TopoHelpers::add_vertex(m);
  vertex_descriptor v2 = TopoHelpers::add_vertex(m);
  vertex_descriptor v3 = TopoHelpers::add_vertex(m);
  vertex_descriptor v4 = TopoHelpers::add_vertex(m);
  vertex_descriptor v5 = TopoHelpers::add_vertex(m);
  vertex_descriptor v6 = TopoHelpers::add_vertex(m);
  vertex_descriptor v7 = TopoHelpers::add_vertex(m);
  vertex_descriptor v8 = TopoHelpers::add_vertex(m);
  vertex_descriptor v9 = TopoHelpers::add_vertex(m);
  vertex_descriptor v10 = TopoHelpers::add_vertex(m);
  vertex_descriptor v11 = TopoHelpers::add_vertex(m);
  vertex_descriptor v12 = TopoHelpers::add_vertex(m);

  // add edges
  edge_descriptor e4 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v1, e4, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v2, e4, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e7 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v7, e7, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v8, e7, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e8 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v4, e8, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v10, e8, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e9 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v0, e9, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v3, e9, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e10 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v4, e10, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v5, e10, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e11 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v5, e11, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v0, e11, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e12 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v1, e12, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v5, e12, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e13 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v1, e13, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v0, e13, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e14 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v6, e14, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v2, e14, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e15 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v0, e15, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v2, e15, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e16 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v3, e16, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v11, e16, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e17 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v11, e17, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v1, e17, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e18 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v12, e18, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v3, e18, TopoHelpers::vertex_pos::SECOND);

  edge_descriptor e19 = TopoHelpers::add_edge(m);
  TopoHelpers::link_vertex_and_edge(v0, e19, TopoHelpers::vertex_pos::FIRST);
  TopoHelpers::link_vertex_and_edge(v12, e19, TopoHelpers::vertex_pos::SECOND);

  // add faces
  face_descriptor f0 = TopoHelpers::add_face(m);
  TopoHelpers::link_edge_and_face(e11, f0);
  TopoHelpers::link_edge_and_face(e12, f0);
  TopoHelpers::link_edge_and_face(e13, f0);

  face_descriptor f4 = TopoHelpers::add_face(m);
  TopoHelpers::link_edge_and_face(e4, f4);
  TopoHelpers::link_edge_and_face(e13, f4);
  TopoHelpers::link_edge_and_face(e15, f4);

  face_descriptor f5 = TopoHelpers::add_face(m);
  TopoHelpers::link_edge_and_face(e9, f5);
  TopoHelpers::link_edge_and_face(e13, f5);
  TopoHelpers::link_edge_and_face(e17, f5);
  TopoHelpers::link_edge_and_face(e16, f5);

  face_descriptor f6 = TopoHelpers::add_face(m);
  TopoHelpers::link_edge_and_face(e18, f6);
  TopoHelpers::link_edge_and_face(e9, f6);
  TopoHelpers::link_edge_and_face(e19, f6);

  // display mesh informations
  std::cout << "Input mesh:\n";
  m->Print();

  // 2-manifoldness of vertex:

  std::cout << "testing  is_2_manifold_vertex(v0)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v0) == false));

  std::cout << "testing  is_2_manifold_vertex(v1)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v1) == false));

  std::cout << "testing  is_2_manifold_vertex(v2)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v2) == false));

  std::cout << "testing  is_2_manifold_vertex(v3)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v3) == true));

  std::cout << "testing  is_2_manifold_vertex(v4)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v4) == false));

  std::cout << "testing  is_2_manifold_vertex(v5)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v5) == false));

  std::cout << "testing  is_2_manifold_vertex(v6)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v6) == false));

  std::cout << "testing  is_2_manifold_vertex(v7)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v7) == false));

  std::cout << "testing  is_2_manifold_vertex(v8)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v8) == false));

  std::cout << "testing  is_2_manifold_vertex(v9)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v9) == false));

  std::cout << "testing  is_2_manifold_vertex(v10)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v10) == false));

  std::cout << "testing  is_2_manifold_vertex(v11)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v11) == true));

  std::cout << "testing  is_2_manifold_vertex(v12)" << std::endl;
  assert((TopoHelpers::is_2_manifold_vertex(v12) == true));

  // 2-manifoldness of edge:

  std::cout << "testing  is_2_manifold_edge(e4)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e4) == false));

  std::cout << "testing  is_2_manifold_edge(e7)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e7) == false));

  std::cout << "testing  is_2_manifold_edge(e8)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e8) == false));

  std::cout << "testing  is_2_manifold_edge(e9)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e9) == false));

  std::cout << "testing  is_2_manifold_edge(e10)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e10) == false));

  std::cout << "testing  is_2_manifold_edge(e11)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e11) == false));

  std::cout << "testing  is_2_manifold_edge(e12)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e12) == false));

  std::cout << "testing  is_2_manifold_edge(e13)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e13) == false));

  std::cout << "testing  is_2_manifold_edge(e14)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e14) == false));

  std::cout << "testing  is_2_manifold_edge(e15)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e15) == false));

  std::cout << "testing  is_2_manifold_edge(e16)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e16) == true));

  std::cout << "testing  is_2_manifold_edge(e17)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e17) == false));

  std::cout << "testing  is_2_manifold_edge(e18)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e18) == true));

  std::cout << "testing  is_2_manifold_edge(e19)" << std::endl;
  assert((TopoHelpers::is_2_manifold_edge(e19) == false));

  // 2-manifoldness of vertex' one-ring:

  std::cout << "testing  is_one_ring_2_manifold(v0)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v0) == false));

  std::cout << "testing  is_one_ring_2_manifold(v1)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v1) == false));

  std::cout << "testing  is_one_ring_2_manifold(v2)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v2) == false));

  std::cout << "testing  is_one_ring_2_manifold(v3)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v3) == false));

  std::cout << "testing  is_one_ring_2_manifold(v4)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v4) == false));

  std::cout << "testing  is_one_ring_2_manifold(v5)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v5) == false));

  std::cout << "testing  is_one_ring_2_manifold(v6)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v6) == false));

  std::cout << "testing  is_one_ring_2_manifold(v7)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v7) == false));

  std::cout << "testing  is_one_ring_2_manifold(v8)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v8) == false));

  std::cout << "testing  is_one_ring_2_manifold(v9)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v9) == false));

  std::cout << "testing  is_one_ring_2_manifold(v10)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v10) == false));

  std::cout << "testing  is_one_ring_2_manifold(v11)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v11) == false));

  std::cout << "testing  is_one_ring_2_manifold(v12)" << std::endl;
  assert((TopoHelpers::is_one_ring_2_manifold(v12) == true));

  std::cout << "Test passed.\n";
  return 0; // all tests passed
}
