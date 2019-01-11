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

  //---------------------- Tests ---------------------

  // Beware!
  // Comparing v1 and v2 is comparing shared pointers
  // comparing *v1 and *v2 is comparing two AIFVertex objects,
  // using AIFVertex class comparison operators.

  // vertex operators
  assert(*v2 == *v2);
  assert(*v2 <= *v2);
  assert(*v2 >= *v2);
  assert(!(*v2 != *v2));
  assert(!(*v2 < *v2));
  assert(!(*v2 > *v2));

  assert(*v2 != *v3);
  assert(!(*v2 == *v3));
  assert((*v2 < *v3) != (*v2 > *v3));
  assert((*v2 <= *v3) != (*v2 >= *v3));

  assert(v2 == v2); // compare smart pointers, ok for vertices !
  assert(v2 != v3); // compare smart pointers, ok for vertices !

  // face operators
  assert(*f0 == *f0);
  assert(*f0 <= *f0);
  assert(*f0 >= *f0);
  assert(!(*f0 != *f0));
  assert(!(*f0 < *f0));
  assert(!(*f0 > *f0));

  assert(*f0 != *f1);
  assert(!(*f0 == *f1));
  assert((*f0 < *f1) != (*f0 > *f1));
  assert((*f0 <= *f1) != (*f0 >= *f1));

  assert(f0 == f0); // compare smart pointers, ok for faces !
  assert(f0 != f1); // compare smart pointers, ok for faces !

  // edge operators
  assert(*e3 == *e3);
  assert(*e3 <= *e3);
  assert(*e3 >= *e3);
  assert(!(*e3 != *e3));
  assert(!(*e3 < *e3));
  assert(!(*e3 > *e3));

  assert(*e3 != *e1);
  assert(!(*e3 == *e1));
  assert((*e3 < *e1) != (*e3 > *e1));
  assert((*e3 <= *e1) != (*e3 >= *e1));

  // Edge comparison is not based on direct pointer comparison
  // but on vertices comparison, so distinct edge objects can be
  // equal.
  edge_descriptor e2twin = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v2, e2twin, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v3, e2twin, helpers::vertex_pos::SECOND);
  edge_descriptor e2reversetwin = helpers::add_edge(m);
  helpers::link_vertex_and_edge(v3, e2reversetwin, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(v2, e2reversetwin, helpers::vertex_pos::SECOND);
  assert(*e2twin == *e2);
  assert(*e2reversetwin == *e2);
  assert(*e2reversetwin == *e2twin);
  assert(*e2twin != *e3);

  assert(e2twin != e2); // compare smart pointers, NOT ok for edges !


  std::cout << "Test passed.\n";
  return 0; // all tests passed
}
