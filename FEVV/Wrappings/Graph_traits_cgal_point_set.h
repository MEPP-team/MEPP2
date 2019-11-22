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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <utility> // for std::make_pair

#include "FEVV/DataStructures/DataStructures_cgal_point_set.h"


#if defined(BOOST_MSVC)
#pragma warning(push)
#pragma warning(disable : 4267)
#endif


namespace boost {

template<>
struct graph_traits< FEVV::CGALPointSet >
{
private:
  typedef  FEVV::CGALPointSetPoint                           Point;

public:
  // Graph
  typedef  FEVV::CGALPointSet::Index                         vertex_descriptor;
  typedef  Point                                             vertex_property_type;

  // HalfedgeGraph
  typedef  int                                               halfedge_descriptor;
    //TODO-elo-fix  halfedge_descriptor needed by SimpleViewer::draw_or_redraw_mesh()
    //              remove when draw_redraw_mesh() is fixed

   // FaceGraph
  typedef  int                                               face_descriptor;
    //TODO-elo-fix  face_decriptor needed par SimpleViewer::internal_createMesh()
    //              remove when internal_createMesh() is fixed
  
  // VertexListGraph
  typedef  FEVV::CGALPointSet::iterator                      vertex_iterator;
  typedef  std::size_t                                       vertices_size_type;

  // nulls
  static vertex_descriptor    null_vertex()   { return -1; }
};

template<>
struct graph_traits< const FEVV::CGALPointSet >
    : public graph_traits< FEVV::CGALPointSet >
{ };

} // namespace boost


namespace CGAL {

  // note:
  //   FEVV::CGALPointSet is a typedef to CGAL::Point_set_3<...> ;
  //   vertices(FEVV::CGALPointSet) must be declared in the same
  //   namespace as CGALPointSet real underlying class for the
  //   name lookup mecanism to work properly ;
  //   see https://en.cppreference.com/w/cpp/language/adl

// Essential free functions specialization for CGALPointSet.

// See http://www.boost.org/doc/libs/1_61_0/libs/graph/doc/graph_concepts.html
// for BGL concepts description.

// See http://doc.cgal.org/latest/BGL/group__PkgBGLConcepts.html
// for CGAL-BGL concepts description.

// BGL VertexListGraph
//!
//! \brief  Returns the iterator range of the vertices of the mesh.
//!

inline
std::pair< typename boost::graph_traits<
                    FEVV::CGALPointSet >::vertex_iterator,
           typename boost::graph_traits<
                    FEVV::CGALPointSet >::vertex_iterator >
vertices(FEVV::CGALPointSet &ps)
{
  return std::make_pair(ps.begin(), ps.end());
}

inline
std::pair< FEVV::CGALPointSet::const_iterator,
           FEVV::CGALPointSet::const_iterator >
vertices(const FEVV::CGALPointSet &ps)
{
  return std::make_pair(ps.begin(), ps.end());
}


// BGL VertexListGraph
//!
//! \brief  Returns an upper bound of the number of vertices of the mesh.
//!
inline
typename boost::graph_traits< FEVV::CGALPointSet >::vertices_size_type
num_vertices(const FEVV::CGALPointSet &ps)
{
  return ps.number_of_points();
}


// BGL VertexMutableGraph Concept
//!
//! \brief  Adds a new vertex to the graph without initializing the
//! connectivity.
//!
inline
typename boost::graph_traits< FEVV::CGALPointSet >::vertex_descriptor
add_vertex(FEVV::CGALPointSet &ps)
{
  FEVV::CGALPointSet::iterator new_point_it = ps.insert();

  // return index of new point
  return new_point_it - ps.begin();
}


// CGAL MutableHalfedgeGraph Concept
//!
//! \brief  Removes v from the mesh.
//!
inline
void
remove_vertex(typename boost::graph_traits<
                  FEVV::CGALPointSet >::vertex_descriptor v,
              FEVV::CGALPointSet &ps)
{
  ps.remove(v);

  // collect_garbage() must be called to really remove the point from container
  // and ensure that the point index is equivalent to the non-removed-point
  // iterator
  ps.collect_garbage();
}


} // namespace CGAL
