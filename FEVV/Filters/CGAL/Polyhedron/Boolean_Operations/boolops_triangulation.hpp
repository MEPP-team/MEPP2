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
 * \file Boolean_Operations_triangulation.h
 * \brief Creates a Constrained Delaunay triangulation to subdivide a facet
 * \author Cyril Leconte
 */

#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

/*!
 * \class Enriched_vertex_base
 * \brief Enriches the vertices of a triangulation
 */
template <class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class Enriched_vertex_base : public Vb
{
public:
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
    typedef Enriched_vertex_base<Gt,Vb2> Other;
  };

private:
  /*! \brief An Id for the vertex*/
  unsigned long m_Label;

public:
  /*! \brief Accessor
   * \param Label : The value to assign*/
  void set_Label(unsigned long Label) {m_Label = Label;}
  /*! \brief Accessor
   * \return The Label of the vertex*/
  unsigned long get_Label() {return m_Label;}
};
  
/*!
 * \class Enriched_face_base
 * \brief Enriches the faces of a triangulation
 */
template <class Gt, class Fb = CGAL::Constrained_triangulation_face_base_2<Gt> >
class Enriched_face_base : public Fb
{
public:
  /*!
   * \typedef typename Vertex_handle
   * \brief Handle for a vertex of a triangulation
   */
  typedef typename Fb::Triangulation_data_structure::Vertex_handle          Vertex_handle;

  /*!
   * \typedef typename Face_handle
   * \brief Handle for a face of a triangulation
   */
  typedef typename Fb::Triangulation_data_structure::Face_handle            Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Enriched_face_base<Gt,Fb2>            Other;
  };

private:
  /*! \brief True if the vertex has been processed*/
  bool m_OK;
  /*! \brief True if the vertex belongs to the result*/
  bool m_Ext;

public:
  Enriched_face_base() : Fb() {}
  Enriched_face_base(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2) : Fb(v0,v1,v2) {}
  Enriched_face_base(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
    Face_handle n0, Face_handle n1, Face_handle n2) : Fb(v0,v1,v2,n0,n1,n2) {}
  Enriched_face_base(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Face_handle n0, Face_handle n1, Face_handle n2,
          bool c0, bool c1, bool c2 ) : Fb(v0,v1,v2,n0,n1,n2) {}
  /*! \brief Accessor
   * \param Ext : The value to assign*/
  void set_Ext(bool Ext) {m_Ext = Ext;}
  /*! \brief Accessor
   * \return true if the triangle of the triangulation belongs to the result*/
  bool get_Ext() {return m_Ext;}
  /*! \brief Accessor
   * \param OK : The value to assign*/
  void set_OK(bool OK) {m_OK = OK;}
  /*! \brief Accessor
   * \return true if the parameter m_Ext has been determined*/
  bool get_OK() {return m_OK;}
};
  
/*!
 * \class Triangulation
 * \brief To subdivide a facet. (the kernel K must be exact)
 */
template <class K>
class Triangulation
{
  /*!
   * \typedef typename Point_3
   * \brief 3d point using exact number type
   */
  typedef typename K::Point_3                            Point_3;
  
  /*!
   * \typedef typename Tri_vb
   * \brief Vertex base
   */
        typedef /*typename*/ Enriched_vertex_base<K>                    Tri_vb;
  
  /*!
   * \typedef typename Tri_fb
   * \brief Face base
   */
        typedef /*typename*/ Enriched_face_base<K>                      Tri_fb;
  
  /*!
   * \typedef typename Tri_DS
   * \brief Data structure of the triangulation
   */
  typedef typename CGAL::Triangulation_data_structure_2<Tri_vb,Tri_fb>      Tri_DS;
  
  /*!
   * \typedef typename Itag
   * \brief No intersection tag
   */
  typedef typename CGAL::No_intersection_tag                    Itag;
  
  /*!
   * \typedef typename Constrained_Delaunay_tri
   * \brief 2d constrained Delaunay triangulation
   */
  typedef typename CGAL::Constrained_Delaunay_triangulation_2<K, Tri_DS, Itag>  Constrained_Delaunay_tri;
  
  /*!
   * \typedef typename Vertex_handle_tri
   * \brief Vertex handle for the triangulation
   */
  typedef typename Constrained_Delaunay_tri::Vertex_handle            Vertex_handle_tri;
  
  /*!
   * \typedef typename Face_handle_tri
   * \brief Face handle for the triangulation
   */
  typedef typename Constrained_Delaunay_tri::Face_handle              Face_handle_tri;
  
  /*!
   * \typedef typename Point_tri
   * \brief 2d point for the triangulation
   */
  typedef typename Constrained_Delaunay_tri::Point                Point_tri;
  
  /*!
   * \typedef typename Face_iterator_tri
   * \brief Iterator for the faces of the triangulation
   */
  typedef typename Constrained_Delaunay_tri::Face_iterator            Face_iterator_tri;

public:
  /*!
   * \brief Constructor
   * \param he : A halfedge incident to the facet
   * \param norm_dir : The vector directing the normal of the facet
   */
  Triangulation(Halfedge_handle &he, Vector_exact &norm_dir)
  {
    //find the longest coordinate of the normal vector, and its sign
    double x = to_double(norm_dir.x());
    double y = to_double(norm_dir.y());
    double z = to_double(norm_dir.z());
    double absx = std::abs(x);
    double absy = std::abs(y);
    double absz = std::abs(z);

    //this information is stored using a code :
    //0 : The coordinate X is the longest, and positive
    //1 : The coordinate Y is the longest, and positive
    //2 : The coordinate Z is the longest, and positive
    //3 : The coordinate X is the longest, and negative
    //4 : The coordinate Y is the longest, and negative
    //5 : The coordinate Z is the longest, and negative
    if (absx >= absy && absx >= absz) max_coordinate = (x>0)?0:3;
    else if (absy >= absx && absy >= absz) max_coordinate = (y>0)?1:4;
    else if (absz >= absx && absz >= absy) max_coordinate = (z>0)?2:5;

    //we add the three vertices of the facet to the triangulation
    //The Label of these vertices is set for the corresponding point in the triangulation
    v1 = add_new_pt(point_to_exact(he->vertex()->point()), he->vertex()->Label);
    v2 = add_new_pt(point_to_exact(he->next()->vertex()->point()), he->next()->vertex()->Label);
    v3 = add_new_pt(point_to_exact(he->next()->next()->vertex()->point()), he->next()->next()->vertex()->Label);

    //if the vertices does not have an Id (Label = OxFFFFFFFF), the labels
    //of the points in the triangulation is set as follows : 
    //0xFFFFFFFF for the first point
    //0xFFFFFFFE for the second point
    //0xFFFFFFFD for the third point
    if(v2->get_Label() == 0xFFFFFFFF) v2->set_Label(0xFFFFFFFE);
    if(v3->get_Label() == 0xFFFFFFFF) v3->set_Label(0xFFFFFFFD);
  }
  
  /*!
   * \brief Compute the orthogonal projection of a point to the plane defined by the longest coordinate of the normal vector
   * \param p : the point (in 3d)
   * \return The projection as a 2d point
   */
  Point_tri get_minvar_point_2(Point_3 &p)
  {
    switch(max_coordinate)
    {
    case 0:
      return Point_tri(p.y(),p.z());
      break;
    case 1:
      return Point_tri(p.z(),p.x());
      break;
    case 2:
      return Point_tri(p.x(),p.y());
      break;
    case 3:
      return Point_tri(p.z(),p.y());
      break;
    case 4:
      return Point_tri(p.x(),p.z());
      break;
    case 5:
      return Point_tri(p.y(),p.x());
      break;
    default:
      return Point_tri(p.y(),p.z());
    }
  }
  
  /*!
   * \brief Adds a point in the triangulation
   * \param p : The point (in 3d : the projection in 2d is done automatically) 
   * \param Label : The label of the point
   * \return The Vertex_handle of the point added 
   */
        Vertex_handle_tri add_new_pt(Point_3 p, unsigned long &Label)   // MT: suppression r�f�rence
  {
    //if the point is not a new one, we verify that the point has not already been added
    if(Label != 0xFFFFFFFF)
      for(unsigned int i = 0;i != pts_point.size();++i)
        if(Label == pts_point[i])
          //if the point is already in the triangulation, we return its handle
          return pts_vertex[i];
    Vertex_handle_tri v;
    v = ct.insert(get_minvar_point_2(p));
    v->set_Label(Label);
    pts_point.push_back(Label);
    pts_vertex.push_back(v);
    return v;
  }
  
  /*!
   * \brief Adds a constrained segment in the triangulation
   * \param p1 : The first point (in 3d)
   * \param p2 : The second point (in 3d)
   * \param Label1 : The label of the first point
   * \param Label2 : The label of the second point
   */
  void add_segment(Point_3 &p1, Point_3 &p2, unsigned long &Label1, unsigned long &Label2)
  {
    // we add the two points in the triangulation and store their handles in c1 and c2
    c1 = add_new_pt(p1, Label1);
    c2 = add_new_pt(p2, Label2);
    // we set a constrained segment between these two points
    ct.insert_constraint(c1, c2);
    // if an other segment is added, we will overwrite c1 and c2
    // what is important is to memorize the handles of the last segment added
  }
  
  /*!
   * \brief Gets the triangles of the triangulation that belongs to the result
   * and deduce for the three neighboring facets if they belong to the result or not
   * \param inv_triangles : must be true if the orientation of the triangles must be inverted
   * \param IsExt : Pointer on a three-case boolean table. 
   * \return The list of the triangles belonging to the result.
   * each triangle is defined by a list of three labels
   */
  std::vector< std::vector<unsigned long> > get_triangles(bool inv_triangles, bool *IsExt)
  {
    //init
    IsExt[0] = false;
    IsExt[1] = false;
    IsExt[2] = false;
    std::vector<std::vector<unsigned long> > tris;
    for(Face_iterator_tri fi = ct.faces_begin();fi != ct.faces_end();fi++)
      fi->set_OK(false);

    //the constrained segments are oriented, we search the triangle (c1, c2, X), (X, c1, c2) or (c2, X, c1)
    //where c1 and c2 are the two points related to the last constrained segment added, and X another point
    //this triangle belongs to the result (thanks to the orientation of the segments)
    Face_handle_tri f, f2 = c1->face();
                int i=0; // MT
    do {
      f = f2;
      f->has_vertex(c1,i);
      f2 = f->neighbor(f->ccw(i));
    } while( ! ( f->has_vertex(c2) && f2->has_vertex(c2) ) );

    //dans le cas particulier ou la frontiere se trouve exactement sur un bord de la triangulation,
    //et que ce triangle n'appartient pas a la triangulation, on d�marrera avec l'autre triangle
    //incluant le segment c1, c2 (et donc, n'appartenant pas au r�sultat

    //if the segment is exactly on the border of the triangulation, the triangle could be outside the triangulation
    //in that case, we will search the other triangle including the points c1 and c2
    //this triangle does not belong to the result
    if(f->has_vertex(ct.infinite_vertex()))
    {
      f = f2;
      f->set_Ext(false);
    }
    else
    {
      f->set_Ext(true);
    }

    std::stack<Face_handle_tri> sfh;
    f->set_OK(true);
    sfh.push(f);

    //we decide for all the triangles, if they belongs to the result, starting from the first triangle f,
    //by moving on the triangulation using the connectivity between the triangles.
    //If a constrained segment is crossed, the value of the tag "isext" is inverted
    while(!sfh.empty())
    {
      f = sfh.top();
      sfh.pop();

      if(f->get_Ext())
      {
        std::vector<unsigned long> tri;
        int i;
        tri.push_back(f->vertex(0)->get_Label());

        //verify if the neighboring facets belongs to the result or not
        if(f->has_vertex(v1,i) && f->neighbor(f->ccw(i))->has_vertex(ct.infinite_vertex())) IsExt[0] = true;
        if(f->has_vertex(v2,i) && f->neighbor(f->ccw(i))->has_vertex(ct.infinite_vertex())) IsExt[1] = true;
        if(f->has_vertex(v3,i) && f->neighbor(f->ccw(i))->has_vertex(ct.infinite_vertex())) IsExt[2] = true;
        
        if(inv_triangles)
        {
          tri.push_back(f->vertex(2)->get_Label());
          tri.push_back(f->vertex(1)->get_Label());
        }
        else
        {
          tri.push_back(f->vertex(1)->get_Label());
          tri.push_back(f->vertex(2)->get_Label());
        }
        tris.push_back(tri);
      }
      for(i = 0;i!=3;i++)
      {
        if(!(f->neighbor(i)->get_OK() || f->neighbor(i)->has_vertex(ct.infinite_vertex())))
        {
          f->neighbor(i)->set_OK(true);
          f->neighbor(i)->set_Ext((f->is_constrained(i))?!f->get_Ext():f->get_Ext());
          sfh.push(f->neighbor(i));
        }
      }
    }
    return tris;
  }
  
  /*!
   * \brief Gets all the triangles of the triangulation
   * \param inv_triangles : must be true if the orientation of the triangles must be inverted
   * \return The list of the triangles belonging to the result.
   * each triangle is defined by a list of three labels
   */
  std::vector< std::vector<unsigned long> > get_all_triangles(bool inv_triangles)
  {
    std::vector< std::vector<unsigned long> > tris;
    for(Face_iterator_tri f = ct.faces_begin();f != ct.faces_end();f++)
    {
      std::vector<unsigned long> tri;
      tri.push_back(f->vertex(0)->get_Label());
      if(inv_triangles)
      {
        tri.push_back(f->vertex(2)->get_Label());
        tri.push_back(f->vertex(1)->get_Label());
      }
      else
      {
        tri.push_back(f->vertex(1)->get_Label());
        tri.push_back(f->vertex(2)->get_Label());
      }
      tris.push_back(tri);
    }
    return tris;
  }
  
private:
  /*! \brief The triangulation*/
  Constrained_Delaunay_tri ct;
  /*! \brief List of the id of the points added in the triangulation*/
  std::vector<InterId> pts_point;
  /*! \brief List of the handles of the points added in the triangulation*/
  std::vector<Vertex_handle_tri> pts_vertex;
  /*! \brief Handle of the point corresponding to the first vertex of the facet*/
  Vertex_handle_tri v1;
  /*! \brief Handle of the point corresponding to the second vertex of the facet*/
  Vertex_handle_tri v2;
  /*! \brief Handle of the point corresponding to the third vertex of the facet*/
  Vertex_handle_tri v3;
  /*! \brief Handle of the point corresponding to the first vertex of the last segment added*/
  Vertex_handle_tri c1;
  /*! \brief Handle of the point corresponding to the second vertex of the last segment added*/
  Vertex_handle_tri c2;
  /*! \brief Code identifying the plane where the triangulation is done \n
   * 0 : Plane (y, z) \n
   * 1 : Plane (z, x) \n
   * 2 : Plane (x, y) \n
   * 3 : Plane (z, y) \n
   * 4 : Plane (x, z) \n
   * 5 : Plane (y, x)
   */
  int max_coordinate;
};
  
