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
#include "FEVV/DataStructures/AIF/AIFMesh.hpp"

#include <cassert>

using namespace FEVV::DataStructures::AIF;
typedef AIFTopologyHelpers helpers;

/*
        Example topology used in this test

               v2 ----------- v6
              /|\              \
             / | \              \
            /  |  \              \
           /   |   \              \
          /    |    \              \
         /     |     \              \
        /      |      \              \
       /       |       \              \
      /        |        \              \		   v9
     /         |         \              \         / \
    /          |          \              \		   /   \
   /     f4    |    f3     \      f2      \		  /     \
  /            |            \              \   /   f5  \
 /             |             \              \ /         \
v1 ---------- v0 ----------- v3             v7 --------- v8
 \             |              |			        /
  \            |              |			       /
   \           |			        |			      /
    \          |			        |		       /
     \   f0    |	  f1        |		      /
      \        |			        |		     /
       \       |			        |		    /			 v11
(isolated) \   |			        |      /
         \     |			        |     /
          \    |			        |    /
           \   |			        |   /
            \  |			        |  /
             \ |			        | /
              \|			        |/
                              v5 ----------- v4 ----------- v10
*/

helpers::face_descriptor
add_face_to_mesh(const std::vector< helpers::vertex_descriptor > &vertices,
                 helpers::smart_ptr_mesh mesh)
{
  typedef std::vector< helpers::vertex_descriptor >::const_iterator iterator;

  helpers::vertex_descriptor prev_vertex = vertices.back(), current_vertex;
  helpers::edge_descriptor current_edge;
  helpers::face_descriptor current_face = helpers::add_face(mesh);

  for(iterator it = vertices.cbegin(); it != vertices.cend(); ++it)
  {
    current_vertex = *it;

    // check if the 2 vertices ain't already link by an edge
    helpers::edge_descriptor current_edge;
    if((current_edge = helpers::common_edge(prev_vertex, current_vertex)) ==
       helpers::null_edge())
    {
      // create the edge in the mesh
      current_edge = helpers::add_edge(mesh);
      // link the edge and the vertices
      helpers::link_vertex_and_edge(
          prev_vertex, current_edge, helpers::vertex_pos::FIRST);
      helpers::link_vertex_and_edge(
          current_vertex, current_edge, helpers::vertex_pos::SECOND);
    }
    // link the edge and the face
    helpers::link_edge_and_face(current_edge, current_face);

    prev_vertex = current_vertex;
  }

  return current_face;
}

unsigned int nbAssert = 0;
void
test_assert(bool expr)
{
  if(!expr)
    throw std::runtime_error("Test Failed");
  nbAssert++;
}

int
main(int narg, char **argv)
{
  helpers::smart_ptr_mesh mesh = helpers::mesh_type::New();

  //// Vertices creation
  const unsigned int nb_vertices = 12;
  helpers::vertex_descriptor vertices[nb_vertices];

  for(unsigned int i = 0; i < nb_vertices; ++i)
    vertices[i] = helpers::add_vertex(mesh);

  //// Dangling Edge creation
  helpers::edge_descriptor dedge = helpers::add_edge(mesh);
  helpers::link_vertex_and_edge(vertices[4], dedge, helpers::vertex_pos::FIRST);
  helpers::link_vertex_and_edge(
      vertices[10], dedge, helpers::vertex_pos::SECOND);

  //// Faces creation
  const unsigned int nb_faces = 6;
  helpers::face_descriptor faces[nb_faces];

  std::vector< helpers::vertex_descriptor > face0 = {
      vertices[0], vertices[5], vertices[1]};
  faces[0] = add_face_to_mesh(face0, mesh);

  std::vector< helpers::vertex_descriptor > face1 = {
      vertices[0], vertices[3], vertices[4], vertices[5]};
  faces[1] = add_face_to_mesh(face1, mesh);

  std::vector< helpers::vertex_descriptor > face2 = {
      vertices[6], vertices[7], vertices[4], vertices[3], vertices[2]};
  faces[2] = add_face_to_mesh(face2, mesh);

  std::vector< helpers::vertex_descriptor > face3 = {
      vertices[0], vertices[2], vertices[3]};
  faces[3] = add_face_to_mesh(face3, mesh);

  std::vector< helpers::vertex_descriptor > face4 = {
      vertices[0], vertices[1], vertices[2]};
  faces[4] = add_face_to_mesh(face4, mesh);

  std::vector< helpers::vertex_descriptor > face5 = {
      vertices[9], vertices[8], vertices[7]};
  faces[5] = add_face_to_mesh(face5, mesh);

  //// Testing
  try
  {
    test_assert(helpers::degree(vertices[0]) == 4);
    test_assert(helpers::degree(vertices[3]) == 3);
    test_assert(helpers::degree(vertices[6]) == 2);

    test_assert(!helpers::is_isolated_vertex(vertices[0]));
    test_assert(helpers::is_isolated_vertex(vertices[11]));

    test_assert(!helpers::is_cut_vertex(vertices[1]));
    test_assert(!helpers::is_cut_vertex(vertices[3]));
    // test_assert(helpers::is_cut_vertex(vertices[7]));

    test_assert(!helpers::is_surface_border_vertex(vertices[0]));
    test_assert(!helpers::is_surface_border_vertex(vertices[3]));
    test_assert(helpers::is_surface_border_vertex(vertices[2]));
    test_assert(helpers::is_surface_border_vertex(vertices[5]));

    test_assert(helpers::is_surface_interior_vertex(vertices[0]));
    test_assert(helpers::is_surface_interior_vertex(vertices[3]));
    test_assert(!helpers::is_surface_interior_vertex(vertices[2]));
    test_assert(!helpers::is_surface_interior_vertex(vertices[5]));

    test_assert(!helpers::are_incident(vertices[0], dedge));
    test_assert(helpers::are_incident(vertices[4], dedge));

    helpers::edge_descriptor common_edge =
        helpers::common_edge(vertices[0], vertices[7]);
    test_assert(common_edge == helpers::null_edge());

    common_edge = helpers::common_edge(vertices[0], vertices[3]);
    test_assert(common_edge != helpers::null_edge());

    test_assert(helpers::are_incident(vertices[0], common_edge));
    test_assert(helpers::are_incident(vertices[3], common_edge));

    test_assert(helpers::are_incident(vertices[9], faces[5]));
    test_assert(helpers::are_incident(vertices[1], faces[4]));
    test_assert(!helpers::are_incident(vertices[5], faces[3]));

    test_assert(helpers::are_adjacent(vertices[0], vertices[3]));
    test_assert(helpers::are_adjacent(vertices[0], vertices[2]));
    test_assert(helpers::are_adjacent(vertices[2], vertices[0]));
    test_assert(!helpers::are_adjacent(vertices[1], vertices[7]));
    test_assert(!helpers::are_adjacent(vertices[1], vertices[1])); // same
                                                                   // vertex

    typedef boost::iterator_range<
        helpers::edge_container_in_vertex::const_iterator >
        ve_range;
    ve_range r1 = helpers::incident_edges(vertices[0]);

    std::vector< unsigned int > vv_vec = {1, 2, 3, 5};
    for(std::vector< unsigned int >::iterator it = vv_vec.begin();
        it != vv_vec.end();
        ++it)
    {
      test_assert(std::find(r1.begin(),
                            r1.end(),
                            helpers::common_edge(vertices[0], vertices[*it])) !=
                  r1.end());
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////

    test_assert(
        helpers::degree(helpers::common_edge(vertices[0], vertices[3])) == 2);
    test_assert(
        helpers::degree(helpers::common_edge(vertices[2], vertices[6])) == 1);
    test_assert(
        helpers::degree(helpers::common_edge(vertices[4], vertices[10])) == 0);

    test_assert(!helpers::is_surface_border_edge(
        helpers::common_edge(vertices[2], vertices[3])));
    test_assert(!helpers::is_surface_border_edge(
        helpers::common_edge(vertices[0], vertices[1])));
    test_assert(!helpers::is_surface_border_edge(
        helpers::common_edge(vertices[4], vertices[10])));
    test_assert(helpers::is_surface_border_edge(
        helpers::common_edge(vertices[4], vertices[7])));

    test_assert(helpers::is_surface_interior_edge(
        helpers::common_edge(vertices[2], vertices[3])));
    test_assert(helpers::is_surface_interior_edge(
        helpers::common_edge(vertices[0], vertices[1])));
    test_assert(!helpers::is_surface_interior_edge(
        helpers::common_edge(vertices[4], vertices[10])));
    test_assert(!helpers::is_surface_interior_edge(
        helpers::common_edge(vertices[4], vertices[7])));

    test_assert(!helpers::is_dangling_edge(
        helpers::common_edge(vertices[2], vertices[3])));
    test_assert(!helpers::is_dangling_edge(
        helpers::common_edge(vertices[0], vertices[1])));
    test_assert(helpers::is_dangling_edge(
        helpers::common_edge(vertices[4], vertices[10])));
    test_assert(!helpers::is_dangling_edge(
        helpers::common_edge(vertices[4], vertices[7])));

    test_assert(!helpers::are_adjacent(
        helpers::common_edge(vertices[4], vertices[7]),
        helpers::common_edge(vertices[4], vertices[7]))); // same edge
    test_assert(
        !helpers::are_adjacent(helpers::common_edge(vertices[7], vertices[8]),
                               helpers::common_edge(vertices[0], vertices[1])));
    test_assert(
        helpers::are_adjacent(helpers::common_edge(vertices[0], vertices[5]),
                              helpers::common_edge(vertices[4], vertices[5])));

    test_assert(helpers::are_incident(
        helpers::common_edge(vertices[2], vertices[3]), faces[2]));
    test_assert(helpers::are_incident(
        helpers::common_edge(vertices[2], vertices[3]), faces[3]));
    test_assert(!helpers::are_incident(
        helpers::common_edge(vertices[2], vertices[3]), faces[4]));

    test_assert(
        helpers::vertex_position(helpers::common_edge(vertices[2], vertices[3]),
                                 vertices[2]) != helpers::vertex_pos::FIRST);
    test_assert(
        helpers::vertex_position(helpers::common_edge(vertices[2], vertices[3]),
                                 vertices[3]) == helpers::vertex_pos::FIRST);
    test_assert(
        helpers::vertex_position(helpers::common_edge(vertices[2], vertices[3]),
                                 vertices[2]) == helpers::vertex_pos::SECOND);
    test_assert(
        helpers::vertex_position(helpers::common_edge(vertices[2], vertices[3]),
                                 vertices[3]) != helpers::vertex_pos::SECOND);

    test_assert(
        helpers::opposite_vertex(helpers::common_edge(vertices[3], vertices[4]),
                                 vertices[3]) == vertices[4]);
    test_assert(
        helpers::opposite_vertex(helpers::common_edge(vertices[3], vertices[4]),
                                 vertices[4]) == vertices[3]);

    helpers::swap_vertices(helpers::common_edge(vertices[3], vertices[4]));

    test_assert(
        helpers::opposite_vertex(helpers::common_edge(vertices[3], vertices[4]),
                                 vertices[3]) == vertices[4]);
    test_assert(
        helpers::opposite_vertex(helpers::common_edge(vertices[3], vertices[4]),
                                 vertices[4]) == vertices[3]);

    test_assert(
        helpers::vertex_position(helpers::common_edge(vertices[3], vertices[4]),
                                 vertices[4]) == helpers::vertex_pos::FIRST);
    test_assert(
        helpers::vertex_position(helpers::common_edge(vertices[3], vertices[4]),
                                 vertices[3]) == helpers::vertex_pos::SECOND);

    helpers::swap_vertices(helpers::common_edge(vertices[3], vertices[4]));

    test_assert(
        helpers::vertex_position(helpers::common_edge(vertices[3], vertices[4]),
                                 vertices[4]) == helpers::vertex_pos::SECOND);
    test_assert(
        helpers::vertex_position(helpers::common_edge(vertices[3], vertices[4]),
                                 vertices[3]) == helpers::vertex_pos::FIRST);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    test_assert(helpers::degree(faces[0]) == 3);
    test_assert(helpers::degree(faces[1]) == 4);
    test_assert(helpers::degree(faces[2]) == 5);

    test_assert(!helpers::are_adjacent(faces[0], faces[0])); // same face
    test_assert(helpers::are_adjacent(faces[0], faces[1]));
    test_assert(helpers::are_adjacent(faces[1], faces[2]));
    test_assert(!helpers::are_adjacent(faces[2], faces[5]));
    test_assert(helpers::are_adjacent(faces[2], faces[3]));

    ////////////////////////////////////////////////////////////////////////////////////////////////

    test_assert(helpers::num_vertices(mesh) == 12);
    test_assert(helpers::num_edges(mesh) == 16);
    test_assert(helpers::num_faces(mesh) == 6);

    typedef boost::iterator_range< helpers::vertex_container::const_iterator >
        v_range;
    typedef boost::iterator_range< helpers::edge_container::const_iterator >
        e_range;
    typedef boost::iterator_range< helpers::face_container::const_iterator >
        f_range;
    typedef helpers::vertex_container_in_edge ve_pair;
    typedef boost::iterator_range<
        helpers::face_container_in_edge::const_iterator >
        fe_range;
    e_range r2 = helpers::edges(mesh);

    for(e_range::iterator it = r2.begin(); it != r2.end(); ++it)
    {
      helpers::unlink_all_faces(*it);
      helpers::unlink_all_vertices(*it);
    }

    test_assert(helpers::num_vertices(mesh) == 12);
    test_assert(helpers::num_edges(mesh) == 16);
    test_assert(helpers::num_faces(mesh) == 6);

    v_range r3 = helpers::vertices(mesh);
    e_range r4 = helpers::edges(mesh);
    f_range r5 = helpers::faces(mesh);

    for(v_range::iterator it = r3.begin(); it != r3.end(); ++it)
      test_assert(helpers::degree(*it) == 0);

    for(e_range::iterator it = r4.begin(); it != r4.end(); ++it)
    {
      test_assert(helpers::is_degenerated_edge(*it));
      test_assert(helpers::degree(*it) == 0);
    }

    for(f_range::iterator it = r5.begin(); it != r5.end(); ++it)
      test_assert(helpers::degree(*it) == 0);

    std::vector< helpers::vertex_descriptor > vertices_to_delete;
    std::copy(
        r3.begin(),
        r3.end(),
        std::back_inserter(vertices_to_delete)); // the pointers need to be
                                                 // copied to keep iterator valid
    // std::vector<helpers::vertex_descriptor> verticesToDelete(r3.size()); //
    // C++ recommends this solution std::copy(r3.begin(), r3.end(),
    // verticesToDelete.begin()); // C++ recommends this solution
    for(std::vector< helpers::vertex_descriptor >::const_iterator it =
            vertices_to_delete.cbegin();
        it != vertices_to_delete.cend();
        ++it)
      helpers::remove_vertex(*it, mesh);

    std::vector< helpers::edge_descriptor > edges_to_delete;
    std::copy(
        r4.begin(),
        r4.end(),
        std::back_inserter(edges_to_delete)); // the pointers need to be copied
                                              // to keep iterator valid
    // std::vector<helpers::edge_descriptor> edgesToDelete(r4.size()); // C++
    // recommends this solution std::copy(r4.begin(), r4.end(),
    // edgesToDelete.begin()); // C++ recommends this solution
    for(std::vector< helpers::edge_descriptor >::const_iterator it =
            edges_to_delete.cbegin();
        it != edges_to_delete.cend();
        ++it)
      helpers::remove_edge(*it, mesh);

    std::vector< helpers::face_descriptor > faces_to_delete;
    std::copy(
        r5.begin(),
        r5.end(),
        std::back_inserter(faces_to_delete)); // the pointers need to be copied
                                              // to keep iterator valid
    // std::vector<helpers::face_descriptor> facesToDelete(r5.size()); // C++
    // recommends this solution std::copy(r5.begin(), r5.end(),
    // facesToDelete.begin()); // C++ recommends this solution
    for(std::vector< helpers::face_descriptor >::const_iterator it =
            faces_to_delete.cbegin();
        it != faces_to_delete.cend();
        ++it)
      helpers::remove_face(*it, mesh);

    test_assert(helpers::num_vertices(mesh) == 0);
    test_assert(helpers::num_edges(mesh) == 0);
    test_assert(helpers::num_faces(mesh) == 0);
  }
  catch(std::runtime_error &e)
  {
    std::cout << e.what() << " (at assert #" << nbAssert << ")" << std::endl;
    return -1;
  }

  std::cout << "Test Successful" << std::endl;
  return 0;
}
