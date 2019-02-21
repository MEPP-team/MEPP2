// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

#pragma once

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_3.h>

#include "boolops_properties.h"


template< class Refs, class T, class P >
class EnrichedVertex : public CGAL::HalfedgeDS_vertex_base< Refs, T, P >
{
public:
  // life cycle
  EnrichedVertex() {}

  // repeat mandatory constructors
  EnrichedVertex(const P &pt)
      : CGAL::HalfedgeDS_vertex_base< Refs, T, P >(pt)
  {
  }

  // vertex properties
  VertexId Label;
};


template< class Refs, class Tprev, class Tvertex, class Tface >
class EnrichedHalfedge
    : public CGAL::HalfedgeDS_halfedge_base< Refs, Tprev, Tvertex, Tface >
{
public:
  // life cycle
  EnrichedHalfedge() {}

  // halfedge properties
  HalfedgeId Label;
};


template< class Refs, class T >
class EnrichedFacet : public CGAL::HalfedgeDS_face_base< Refs, T >
{
public:
  // life cycle
  EnrichedFacet() {}

  // face properties
  /*! \brief true if the facet belongs to the result*/
  bool IsExt;
  /*! \brief true if the facet has been processed*/
  bool IsOK;
  /*! \brief An Id for the facet*/
  FacetId Label;
};


struct EnrichedItems : public CGAL::Polyhedron_items_3
{
  // wrap vertex
  template< class Refs, class Traits >
  struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Normal;
    typedef EnrichedVertex< Refs, CGAL::Tag_true, Point > Vertex;
  };

  // wrap halfedge
  template< class Refs, class Traits >
  struct Halfedge_wrapper
  {
    typedef typename Traits::Vector_3 Normal;
    typedef EnrichedHalfedge< Refs, CGAL::Tag_true, CGAL::Tag_true, CGAL::Tag_true >
        Halfedge;
  };

  // wrap face
  template< class Refs, class Traits >
  struct Face_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Normal;
    typedef typename Traits::Plane_3 Plane;
    typedef EnrichedFacet< Refs, CGAL::Tag_true > Face;
  };
};


using CGALKernel = CGAL::Cartesian< double >;
using EnrichedPolyhedron =
    CGAL::Polyhedron_3< CGALKernel, EnrichedItems >;
