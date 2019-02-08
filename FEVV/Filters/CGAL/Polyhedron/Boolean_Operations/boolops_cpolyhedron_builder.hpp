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

/*!
 * \file CPolyhedron_from_polygon_builder_3.h
 * \brief An incremental builder to build a polyhedron
 * \author Cyril Leconte
 */

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include "boolops_enriched_polyhedron.hpp"


typedef typename EnrichedPolyhedron::Halfedge_handle Halfedge_handle;
typedef typename EnrichedPolyhedron::Facet_handle    Facet_handle;
typedef typename EnrichedPolyhedron::Point_3         Point3d;

/*!
 * \class CPolyhedron_from_polygon_builder_3
 * \brief A polyhedron incremental builder
 */
template <class HDS>
class CPolyhedron_from_polygon_builder_3 : public CGAL::Modifier_base<HDS> {

public:
  /*!
   * \typedef typename Point_3
   * \brief 3d point of an halfedge data structure
   */
  typedef typename HDS::Traits::Point_3                  Point_3;

  /*!
   * \typedef typename Indices
   * \brief A list of indices (unsigned long) to describe a facet
   */
  typedef typename std::vector<unsigned long>                Indices;

  /*!
   * \typedef typename Builder
   * \brief The polyhedron incremental builder
   */
  typedef typename CGAL::Polyhedron_incremental_builder_3<HDS>      Builder;

private:
  // Member variables
  /*! \brief List of the vertices*/
  std::vector<Point_3>                          m_Sorted_vertices;
  /*! \brief List of the facets*/
  std::vector<Indices>                          m_Facets_indices;

public:
  // Constructors
  /*!
   * \brief Constructor
   */
  CPolyhedron_from_polygon_builder_3() {}
  
  /*!
   * \brief Adds a triangular facet from a facet of a polyhedron
   * \param f : The facet handle
   * \param invert : must be true if the orientation of the facet must be inverted
   */
  void add_triangle(Facet_handle &f, bool invert)
  {
    //initially, the label of the vertices is 0xFFFFFFFF. if a vertex is added to the result, the tag is set to
    //the number of vertices added.
    
    //creation of a list of indices
    Indices  vi;

    //adding the first vertex to the result and adding its label to the list of indices.
    //if the vertex is already added (its label is not 0xFFFFFFFF) we only need to add
    //this label to "vi" without adding the vertex.
    Halfedge_handle he = f->facet_begin();
    if(he->vertex()->Label == 0xFFFFFFFF) add_vertex(he->vertex()->point(), he->vertex()->Label);
    vi.push_back(he->vertex()->Label);

    //the order of the two other vertices depends on the orientation of the facet
    if(!invert)
    {
      if(he->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->vertex()->point(), he->next()->vertex()->Label);
      vi.push_back(he->next()->vertex()->Label);
      if(he->next()->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->next()->vertex()->point(), he->next()->next()->vertex()->Label);
      vi.push_back(he->next()->next()->vertex()->Label);
    }
    else
    {
      if(he->next()->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->next()->vertex()->point(), he->next()->next()->vertex()->Label);
      vi.push_back(he->next()->next()->vertex()->Label);
      if(he->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->vertex()->point(), he->next()->vertex()->Label);
      vi.push_back(he->next()->vertex()->Label);
    }

    //finally, "vi" is added to the list of the facets
    m_Facets_indices.push_back(vi);
  }
  
  /*!
   * \brief Adds a list of triangular facet from an intersected facet of a polyhedron
   * \param T : The list of triangle to add. Each triangle is described as a list of three indices
   * \param he : First halfedge handle of the facet
   */
  void add_triangle(std::vector< std::vector< unsigned long > > &T, Halfedge_handle &he)
  {
    //For each triangle of the vector T...
    for(unsigned int i = 0;i != T.size();++i)
    {
      //we verify that the indices are valid.
      //if one of the indices equals to 0xFFFFFFFF, 0xFFFFFFFE or 0XFFFFFFFD, it means that
      //the corresponding vertex is respectively the first, the second or the third vertex 
      //of the facet.
      for(unsigned int j = 0 ; j != 3 ; ++j)
      {
        switch (T[i][j])
        {
        case 0xFFFFFFFF:
          if(he->vertex()->Label != 0xFFFFFFFF)
          {
            T[i][j] = he->vertex()->Label;
          }
          else
          {
            T[i][j] = static_cast< unsigned long >(m_Sorted_vertices.size());
            he->vertex()->Label = T[i][j];
            m_Sorted_vertices.push_back(he->vertex()->point());
          }
          break;
        case 0xFFFFFFFE:
          if(he->next()->vertex()->Label != 0xFFFFFFFF)
          {
            T[i][j] = he->next()->vertex()->Label;
          }
          else
          {
            T[i][j] = static_cast< unsigned long >(m_Sorted_vertices.size());
            he->next()->vertex()->Label = T[i][j];
            m_Sorted_vertices.push_back(he->next()->vertex()->point());
          }
          break;
        case 0xFFFFFFFD:
          if(he->next()->next()->vertex()->Label != 0xFFFFFFFF)
          {
            T[i][j] = he->next()->next()->vertex()->Label;
          }
          else
          {
            T[i][j] = static_cast< unsigned long >(m_Sorted_vertices.size());
            he->next()->next()->vertex()->Label = T[i][j];
            m_Sorted_vertices.push_back(he->next()->next()->vertex()->point());
          }
          break;
        }
      }

      //finally, the facet is added to the list of the facets
      m_Facets_indices.push_back(T[i]);
    }
  }
  
  /*!
   * \brief Adds a Vertex
   * \param p : The point to add
   * \param l : The corresponding label
   */
  void add_vertex(Point3d p, unsigned long &l) // MT: suppression référence
  {
    //The value of the label is updated
    l = static_cast< unsigned long >(m_Sorted_vertices.size());
    //The vertex is added
    m_Sorted_vertices.push_back(p);
  }
  
  /*!
   * \brief this method builds the polyhedron, using the vertices and the facets stored
   * \param hds : The halfedge data structure
   */
  void operator()(HDS& hds)
  {
    Builder B(hds, true);
    B.begin_surface(3,1);
    add_vertices(B);  
    add_facets(B);
    B.end_surface();
  }
  
private:
  
  /*!
   * \brief Used to build the vertices of the polyhedron
   * \param B : The builder
   */
  void add_vertices(Builder& B)
  {
    for(int i = 0; i != (int)this->m_Sorted_vertices.size(); i++)
    {
      B.add_vertex(this->m_Sorted_vertices[i]);
    }
  }
  
  /*!
   * \brief Used to build the facets of the polyhedron
   * \param B : The builder
   */
  void add_facets(Builder &B)
  {
    for(int i = 0; i != (int)this->m_Facets_indices.size(); i++)
    {
      B.begin_facet();
      for(int j = 0; j != (int)this->m_Facets_indices[i].size(); j++)
      {
        B.add_vertex_to_facet(this->m_Facets_indices[i][j]);
      }
      B.end_facet();
    }
  }
};
