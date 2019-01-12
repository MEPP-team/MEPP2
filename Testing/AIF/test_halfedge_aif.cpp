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

// must be the last include to avoid side effect of disabling NDEBUG
#undef NDEBUG
#include <cassert>

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
  typedef helpers::face_descriptor face_descriptor;


  // build mesh
  //
  //         v3
  //        /|\
  //     e3/ | \e2
  //   v0 /f0|f1\ v2
  //      \e4|  /
  //     e0\ | /e1
  //        \|/
  //         v1
  //
  smart_ptr_mesh m = AIFMesh::New();

  // add vertices
  vertex_descriptor v0 = helpers::add_vertex(m);
  vertex_descriptor v1 = helpers::add_vertex(m);
  vertex_descriptor v2 = helpers::add_vertex(m);
  vertex_descriptor v3 = helpers::add_vertex(m);

  // add edges
  edge_descriptor e0 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v0, e0, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v1, e0, helpers::vertex_pos::SECOND);
  edge_descriptor e1 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v1, e1, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v2, e1, helpers::vertex_pos::SECOND);
  edge_descriptor e2 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v2, e2, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v3, e2, helpers::vertex_pos::SECOND);
  edge_descriptor e3 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v3, e3, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v0, e3, helpers::vertex_pos::SECOND);
  edge_descriptor e4 = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v1, e4, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v3, e4, helpers::vertex_pos::SECOND);

  // add faces
  face_descriptor f0 = helpers::add_face(m);
  helpers::link_edge_and_face(e0, f0);
  helpers::link_edge_and_face(e4, f0);
  helpers::link_edge_and_face(e3, f0);
  face_descriptor f1 = helpers::add_face(m);
  helpers::link_edge_and_face(e1, f1);
  helpers::link_edge_and_face(e4, f1);
  helpers::link_edge_and_face(e2, f1);

  // display mesh informations
  std::cout << "Input mesh:\n";
  m->Print();

  // define some halfedges
  AIFHalfEdge h1(v0, v1, f0);
  AIFHalfEdge h2(v1, v3, f0);
  AIFHalfEdge h3(v1, v3, f1);

  try
  {
    AIFHalfEdge hbad(v0, v2, f0); // invalid halfedge, no edge between vertices
    assert(false);
  }
  catch(const std::runtime_error &e)
  {
    assert(e.what() ==
           std::string("In AIFHalfEdge::AIFHalfEdge(s, t, f): invalid "
                       "halfedge, no edge between the two vertices."));
  }

  try
  {
    AIFHalfEdge hbad(v0, v1, f1); // invalid halfedge, v0v1 not incident to f1
    assert(false);
  }
  catch(const std::runtime_error &e)
  {
    assert(e.what() ==
           std::string("In AIFHalfEdge::AIFHalfEdge(s, t, f): invalid "
                       "halfedge, not incident to the face."));
  }

  try
  {
    AIFHalfEdge hbad(v3, v2, f0); // invalid halfedge, v3v2 not incident to f0
    assert(false);
  }
  catch(const std::runtime_error &e)
  {
    assert(e.what() ==
           std::string("In AIFHalfEdge::AIFHalfEdge(s, t, f): invalid "
                       "halfedge, not incident to the face."));
  }

  // test 1
  AIFHalfEdge n;
  n = next(h1, *m);
  assert(n.get_source() == v1 && n.get_target() == v3);

  // test 2
  n = next(h2, *m);
  assert(n.get_source() == v3 && n.get_target() == v0);

  // test 3
  n = next(h3, *m);
  assert(n.get_source() == v3 && n.get_target() == v2);

  std::cout << "Test passed.\n";
  return 0; // all tests passed
}
