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

// must be the last include to avoid side effect of disabling NDEBUG
#undef NDEBUG
#include <cassert>

//----------------------------------------------------------

template< typename T >
void
print_property_map(FEVV::DataStructures::AIF::PropertyMap< T > *pm)
{
  std::cout << "[";
  for(int i = 0; i < pm->size(); i++)
  {
    std::cout << " " << get(*pm, i);
    if(i != pm->size() - 1)
      std::cout << ",";
  }
  std::cout << "]\n";
}

//----------------------------------------------------------

int
main(int argc, const char **argv)
{
  typedef FEVV::DataStructures::AIF::AIFMesh AIFMesh;
  typedef FEVV::DataStructures::AIF::AIFVertex AIFVertex;
  typedef FEVV::DataStructures::AIF::AIFEdge AIFEdge;
  typedef FEVV::DataStructures::AIF::AIFFace AIFFace;
  typedef AIFMesh::Point Point;
  typedef FEVV::DataStructures::AIF::AIFTopologyHelpers AIFHelpers;
  typedef AIFHelpers::vertex_descriptor vertex_descriptor;
  typedef AIFHelpers::edge_descriptor edge_descriptor;
  typedef AIFHelpers::face_descriptor face_descriptor;
  typedef AIFHelpers::smart_ptr_mesh smart_ptr_mesh;
  typedef FEVV::DataStructures::AIF::AIFPropertiesHelpers AIFPropHelpers;


  {
    smart_ptr_mesh m = AIFMesh::New();

    // create mesh with some vertices
    vertex_descriptor v0 = AIFHelpers::add_vertex(m);
    vertex_descriptor v1 = AIFHelpers::add_vertex(m);
    vertex_descriptor v2 = AIFHelpers::add_vertex(m);
    vertex_descriptor v3 = AIFHelpers::add_vertex(m);

    // populate vertex coordinates property map
    AIFPropHelpers::set_point(m, v0, 0.1, 0.2, 0.3);
    AIFPropHelpers::set_point(m, v1, 1.1, 1.2, 1.3);
    AIFPropHelpers::set_point(m, v2, 2.1, 2.2, 2.3);
    AIFPropHelpers::set_point(m, v3, 3.1, 3.2, 3.3);

    // display mesh informations to check vertex coordinates are ok
    std::cout << "Input mesh:\n";
    m->Print();

    // set_point()/get_point() test
    {
      auto v_range = AIFHelpers::vertices(m);
      auto v_iter = v_range.begin();
      auto v_end = v_range.end();
      for(; v_iter != v_end; ++v_iter)
      {
        Point &p = AIFPropHelpers::get_point(m, *v_iter);
        assert(((*v_iter) == v0 && p == Point(0.1, 0.2, 0.3)) ||
               ((*v_iter) == v1 && p == Point(1.1, 1.2, 1.3)) ||
               ((*v_iter) == v2 && p == Point(2.1, 2.2, 2.3)) ||
               ((*v_iter) == v3 && p == Point(3.1, 3.2, 3.3)));
      }
    }

    // add a new property map to the vertices and populate it
    FEVV::DataStructures::AIF::PropertyMap< Point > *test_pm =
        m->AddPropertyMap< AIFVertex::ptr, Point >("v:test");
    if(!m->isPropertyMap< AIFVertex::ptr >("v:test"))
      throw std::runtime_error("Failed to create property map.");

    (*test_pm)[v0] = Point(10.1, 10.2, 10.3);
    (*test_pm)[v1] = Point(11.1, 11.2, 11.3);
    (*test_pm)[v2] = Point(12.1, 12.2, 12.3);

    // test property map #1
    {
      assert(test_pm->size() == 3);

      auto v_range = AIFHelpers::vertices(m);
      auto v_iter = v_range.begin();
      auto v_end = v_range.end();
      for(; v_iter != v_end; ++v_iter)
      {
        Point &p = (*test_pm)[*v_iter];
        assert(((*v_iter) == v0 && p == Point(10.1, 10.2, 10.3)) ||
               ((*v_iter) == v1 && p == Point(11.1, 11.2, 11.3)) ||
               ((*v_iter) == v2 && p == Point(12.1, 12.2, 12.3)) ||
               ((*v_iter) == v3));
      }
    }

    // display property map content
    auto v_range = AIFHelpers::vertices(m);
    auto v_iter = v_range.begin();
    auto v_end = v_range.end();
    for(; v_iter != v_end; ++v_iter)
    {
      Point &p = (*test_pm)[*v_iter];
      std::cout << "[" << (*v_iter)->GetIndex() << "]: " << p[0] << ", " << p[1]
                << ", " << p[2] << '\n';
    }

    // remove property maps
    m->RemovePropertyMap< AIFVertex::ptr >("v:test");
  }

  //--------------------------------------------------------------------------
  // test adding-removing vertices/edges/faces, and effect upon property_maps
  //--------------------------------------------------------------------------
  {
    smart_ptr_mesh m = AIFMesh::New();

    // test property maps on vertices

    FEVV::DataStructures::AIF::PropertyMap< std::string > *vpm_str =
        m->AddPropertyMap< AIFVertex::ptr, std::string >("v:test_str");
    FEVV::DataStructures::AIF::PropertyMap< int > *vpm_int =
        m->AddPropertyMap< AIFVertex::ptr, int >("v:test_int");

    vertex_descriptor v0 = AIFHelpers::add_vertex(m);
    vertex_descriptor vtmp0 = AIFHelpers::add_vertex(m);
    vertex_descriptor v1 = AIFHelpers::add_vertex(m);
    vertex_descriptor v2 = AIFHelpers::add_vertex(m);
    (*vpm_str)[v0] = "v0";
    (*vpm_str)[v2] = "v2";
    (*vpm_str)[v1] = "v1";
    (*vpm_int)[v2] = 2;
    (*vpm_int)[v1] = 1;
    (*vpm_int)[v0] = 0;
    (*vpm_int)[vtmp0] = 9999;
    // vertices = [  v0 , vtmp0,  v1 ,  v2  ]
    // spm_str  = [ "v0",      , "v1", "v2" ]
    // spm_int  = [   0 , 9999 ,   1 ,   2  ]
    assert(vpm_str->size() == 4);
    assert(vpm_int->size() == 4);

    vertex_descriptor v3 = AIFHelpers::add_vertex(m);
    // vertices = [  v0 , vtmp0,  v1 ,  v2 , v3 ]
    // spm_str  = [ "v0",      , "v1", "v2" ]
    // spm_int  = [   0 , 9999 ,   1 ,   2  ]
    assert(vpm_str->size() == 4);
    assert(vpm_int->size() == 4);

    (*vpm_str)[v3] = "v3";
    // vertices = [  v0 , vtmp0,  v1 ,  v2 ,  v3  ]
    // spm_str  = [ "v0",      , "v1", "v2", "v3" ]
    // spm_int  = [   0 , 9999 ,   1 ,   2  ]
    assert(vpm_str->size() == 5);
    assert(vpm_int->size() == 4); // (*vpm_int)[v3] not created

    AIFHelpers::remove_vertex(vtmp0, m);
    // vertices = [  v0 ,  v3 ,  v1 ,  v2  ]
    // spm_str  = [ "v0", "v3", "v1", "v2", "v3" ]
    // spm_int  = [   0 , 9999,   1 ,   2  ]
    // the last element of spm_int is copied in v3 slot, which looks like a bug;
    // let say that after a remove operation, one property map slot is undefined
    // if no value has been set inside before;
    // deletion doesn't free space in property maps because of
    // boost::vector_property_maps
    assert(vpm_int->size() == 4);
    assert(vpm_str->size() == 5);
    assert(get(*vpm_str, 1) == "v3"); // right value (set before)
    assert(get(*vpm_int, 1) == 9999); // undefined value (never set before)

    vertex_descriptor v4 = AIFHelpers::add_vertex(m);
    (*vpm_str)[v4] = "v4";
    (*vpm_int)[v4] = 4;
    // vertices = [  v0 ,  v3 ,  v1 ,  v2 ,  v4  ]
    // spm_str  = [ "v0", "v3", "v1", "v2", "v4" ]
    // spm_int  = [   0 , 9999,   1 ,   2 ,   4  ]
    assert(vpm_int->size() == 5);
    assert(vpm_str->size() == 5);
    assert(get(*vpm_str, 0) == "v0");
    assert(get(*vpm_str, 1) == "v3");
    assert(get(*vpm_str, 2) == "v1");
    assert(get(*vpm_str, 3) == "v2");
    assert(get(*vpm_str, 4) == "v4");
    assert(get(*vpm_int, 0) == 0);
    assert(get(*vpm_int, 1) == 9999);
    assert(get(*vpm_int, 2) == 1);
    assert(get(*vpm_int, 3) == 2);
    assert(get(*vpm_int, 4) == 4);
    assert((*vpm_str)[v0] == "v0");
    assert((*vpm_str)[v1] == "v1");
    assert((*vpm_str)[v2] == "v2");
    assert((*vpm_str)[v3] == "v3");
    assert((*vpm_str)[v4] == "v4");
    assert((*vpm_int)[v0] == 0);
    assert((*vpm_int)[v1] == 1);
    assert((*vpm_int)[v2] == 2);
    assert((*vpm_int)[v3] == 9999); // undefined
    assert((*vpm_int)[v4] == 4);

    // test property maps on edges

    FEVV::DataStructures::AIF::PropertyMap< std::string > *epm_str =
        m->AddPropertyMap< AIFEdge::ptr, std::string >("e:test_str");

    edge_descriptor e0 = AIFHelpers::add_edge(m);
    edge_descriptor e1 = AIFHelpers::add_edge(m);
    edge_descriptor e2 = AIFHelpers::add_edge(m);
    // edges   = [  e0 ,  e1 ,  e2  ]
    // spm_str = [  ]
    AIFHelpers::remove_edge(e0, m);
    // edges   = [  e2 ,  e1  ]
    // spm_str = [  ]
    edge_descriptor e3 = AIFHelpers::add_edge(m);
    edge_descriptor e4 = AIFHelpers::add_edge(m);
    edge_descriptor e5 = AIFHelpers::add_edge(m);
    (*epm_str)[e4] = "e4";
    (*epm_str)[e1] = "e1";
    (*epm_str)[e2] = "e2";
    (*epm_str)[e5] = "e5";
    (*epm_str)[e3] = "e3";
    // edges   = [  e2 ,  e1 ,  e3 ,  e4 ,  e5  ]
    // spm_str = [ "e2", "e1", "e3", "e4", "e5" ]
    AIFHelpers::remove_edge(e5, m);
    // edges   = [  e2 ,  e1 ,  e3 ,  e4  ]
    // spm_str = [ "e2", "e1", "e3", "e4", "e5" ]
    assert(AIFHelpers::num_edges(m) == 4);
    assert(epm_str->size() == 5);
    assert(get(*epm_str, 0) == "e2");
    assert(get(*epm_str, 1) == "e1");
    assert(get(*epm_str, 2) == "e3");
    assert(get(*epm_str, 3) == "e4");
    assert(get(*epm_str, 4) == "e5");
    assert((*epm_str)[e1] == "e1");
    assert((*epm_str)[e2] == "e2");
    assert((*epm_str)[e3] == "e3");
    assert((*epm_str)[e4] == "e4");

    // test property maps on faces

    FEVV::DataStructures::AIF::PropertyMap< int > *fpm =
        m->AddPropertyMap< AIFFace::ptr, int >("f:test");

    face_descriptor f0 = AIFHelpers::add_face(m);
    (*fpm)[f0] = 0;
    face_descriptor f1 = AIFHelpers::add_face(m);
    (*fpm)[f1] = 1;
    face_descriptor f2 = AIFHelpers::add_face(m);
    (*fpm)[f2] = 2;
    face_descriptor f3 = AIFHelpers::add_face(m);
    (*fpm)[f3] = 3;
    face_descriptor f4 = AIFHelpers::add_face(m);
    (*fpm)[f4] = 4;
    // faces   = [ f0, f1, f2, f3, f4 ]
    // fpm     = [  0,  1,  2,  3,  4 ]
    assert(AIFHelpers::num_faces(m) == 5);
    assert(fpm->size() == 5);
    AIFHelpers::remove_face(f1, m);
    // faces   = [ f0, f4, f2, f3 ]
    // fpm     = [  0,  4,  2,  3,  4 ]
    assert(AIFHelpers::num_faces(m) == 4);
    assert(fpm->size() == 5);
    face_descriptor f5 = AIFHelpers::add_face(m);
    (*fpm)[f5] = 5;
    // faces   = [ f0, f4, f2, f3, f5 ]
    // fpm     = [  0,  4,  2,  3,  5 ]
    assert(AIFHelpers::num_faces(m) == 5);
    assert(fpm->size() == 5);
    AIFHelpers::remove_face(f3, m);
    // faces   = [ f0, f4, f2, f5 ]
    // fpm     = [  0,  4,  2,  5,  5 ]
    assert(AIFHelpers::num_faces(m) == 4);
    assert(fpm->size() == 5);
    face_descriptor f6 = AIFHelpers::add_face(m);
    (*fpm)[f6] = 6;
    face_descriptor f7 = AIFHelpers::add_face(m);
    (*fpm)[f7] = 7;
    // faces   = [ f0, f4, f2, f5, f6, f7 ]
    // fpm     = [  0,  4,  2,  5,  6,  7 ]
    AIFHelpers::remove_face(f7, m);
    AIFHelpers::remove_face(f6, m);
    // faces   = [ f0, f4, f2, f5 ]
    // fpm     = [  0,  4,  2,  5,  6,  7 ]
    face_descriptor f8 = AIFHelpers::add_face(m);
    (*fpm)[f8] = 8;
    // faces   = [ f0, f4, f2, f5, f8 ]
    // fpm     = [  0,  4,  2,  5,  8,  7 ]
    AIFHelpers::remove_face(f0, m);
    // faces   = [ f8, f4, f2, f5 ]
    // fpm     = [  7!,  4,  2,  5,  8,  7 ]
    assert(AIFHelpers::num_edges(m) == 4);
    assert(fpm->size() == 6);
    assert(get(*fpm, 0) == 8);
    assert(get(*fpm, 1) == 4);
    assert(get(*fpm, 2) == 2);
    assert(get(*fpm, 3) == 5);
    assert(get(*fpm, 4) == 8);
    assert(get(*fpm, 5) == 7);
    assert((*fpm)[f2] == 2);
    assert((*fpm)[f4] == 4);
    assert((*fpm)[f5] == 5);
    assert((*fpm)[f8] == 8);
  }


  //--------------------------------------------------------------------------
  // test Associative property maps
  //--------------------------------------------------------------------------
  {
    smart_ptr_mesh m = AIFMesh::New();

    // create and populate an associative property map
    FEVV::DataStructures::AIF::AssocPropertyMap< char, float > *pm =
        m->AddAssocPropertyMap< char, float >("assoc:test");

    (*pm)['a'] = 42.3f;
    (*pm)['b'] = 42.5f;
    assert((*pm)['a'] == 42.3f);
    assert((*pm)['b'] == 42.5f);

    m->RemoveAssocPropertyMap("assoc:test");

    try
    {
      pm = m->GetAssocPropertyMap< char, float >("assoc:test");
      // the previous line must trigger an exception
      // so we must not pass here
      assert(false);
    }
    catch(std::runtime_error)
    {
    }
  }


  std::cout << "Test passed.\n";
  return 0;
}
