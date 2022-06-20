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
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"

// DBG  #undef NDEBUG // disable some unwanted assertions in Debug mode
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/internal/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h> // for add_face(vr,g)

#include "FEVV/Types/Mesh_vector_representation.h"


namespace FEVV {
namespace Filters {

/**
 * Create and return a 256x256 "no tex" image.
 */
inline
FEVV::Types::Image* make_no_texture_image(void)
{
  FEVV::Types::Image *cimg = new FEVV::Types::Image(256, 256, 1, 3);

  // draw gradient

  // loop over pixels
  for(int row = 0; row < (*cimg).height(); row++)
  {
    for(int col = 0; col < (*cimg).width(); col++)
    {
      // (*cimg)(x, y, z, channel)
      (*cimg)(col, row, 0, 0) = col; // R = x = col
      (*cimg)(col, row, 0, 1) = row; // G = y = row
      (*cimg)(col, row, 0, 2) = 0;   // B = constant
    }
  }

  // draw text

  const unsigned char gray[] = { 128, 128, 128 };

  int size1 = 72;
  (*cimg).draw_text(20, 128-(size1/2), "NO TEX", gray, 0, 1.0, size1);

  int size2 = 32;
  (*cimg).draw_text(     2,         0, "(0;0)", gray, 0, 1.0, size2);
  (*cimg).draw_text(255-50,         0, "(1;0)", gray, 0, 1.0, size2);
  (*cimg).draw_text(     2, 255-size2, "(0;1)", gray, 0, 1.0, size2);
  (*cimg).draw_text(255-50, 255-size2, "(1;1)", gray, 0, 1.0, size2);

  return cimg;
}

/**
 * Helper type for vertex descriptors vector.
 */
template< typename BoostGraph >
using VertDescVect = std::vector<
    typename boost::graph_traits< BoostGraph >::vertex_descriptor >;


/**
 * Helper type for face descriptors vector.
 */
template< typename FaceGraph >
using FaceDescVect = std::vector<
    typename boost::graph_traits< FaceGraph >::face_descriptor >;


/**
 * Implement the named parameters idiom for mesh_from_vector_representation()
 *
 * \param  use_corner_texture_coord   force use of corner texture
 * \param  vd_target                  a vector to store the vertex descriptors
 * \param  fd_target                  a vector to store the face descriptors
 */
template< typename FaceGraph >
struct MeshFromVectorReprParameters
{
  typedef  VertDescVect< FaceGraph >  VDVect;
  typedef  FaceDescVect< FaceGraph >  FDVect;

  MeshFromVectorReprParameters& use_corner_texcoord(
      bool _use_corner_texcoord)
  {
    m_use_corner_texcoord = _use_corner_texcoord;
    return *this;
  }

  MeshFromVectorReprParameters& vd_target(VDVect *_vd_target)
  {
    m_vd_target = _vd_target;
    return *this;
  }

  MeshFromVectorReprParameters& fd_target(FDVect *_fd_target)
  {
    m_fd_target = _fd_target;
    return *this;
  }

  bool m_use_corner_texcoord;
  VDVect *m_vd_target = nullptr;
  FDVect *m_fd_target = nullptr;
};


/**
 * \brief  Build the mesh from the given vector representation.
 *
 * \param  g           the output mesh
 * \param  pmaps       the output property maps bag of the mesh
 * \param  dup_v_nbr   number of duplicated vertices of the mesh
 * \param  mvr         the vector representation structure (input)
 * \param  params      the named parameters
 * \param  gt          the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 */
template< typename HalfedgeGraph,
          typename coordP_type,
          typename coordN_type,
          typename coordT_type,
          typename coordC_type,
          typename index_type,
          typename GeometryTraits >
void
mesh_from_vector_representation(
    HalfedgeGraph                                       &g,
    FEVV::PMapsContainer                                &pmaps,
    unsigned int                                        &dup_v_nbr,
    FEVV::Types::MVR < coordP_type,
                       coordN_type,
                       coordT_type,
                       coordC_type,
                       index_type >                     &mvr,
    MeshFromVectorReprParameters< HalfedgeGraph > const &params,
    const GeometryTraits                                &/*gt*/)
{
  typedef boost::graph_traits< HalfedgeGraph >       GraphTraits;
  typedef typename GraphTraits::vertex_descriptor    vertex_descriptor;
  //typedef typename GraphTraits::edge_descriptor      edge_descriptor;
  typedef typename GraphTraits::face_descriptor      face_descriptor;
  typedef typename GraphTraits::halfedge_descriptor  halfedge_descriptor;
  //typedef typename GraphTraits::face_iterator        face_iterator;
  typedef typename GeometryTraits::Point             Point;
  typedef typename GeometryTraits::Vector            Vector;

  ///////////////////////// CREATE NEEDED PROPERTY MAPS
  //////////////////////////////

  // vertex-color must be defined for all vertices or none
  if(mvr.vertex_color_coords.size() != mvr.points_coords.size())
    mvr.vertex_color_coords.clear();
  for(auto &color : mvr.vertex_color_coords)
  {
    if(color.size() != 3)
    {
      std::cout << __func__
                << "(): found vertex color with size != 3. Disabling "
                   "vertex color for all vertices."
                << std::endl;
      mvr.vertex_color_coords.clear();
      break;
    }
  }

  // face-color must be defined for all faces or none
  if(mvr.face_color_coords.size() != mvr.faces_indices.size())
    mvr.face_color_coords.clear();
  for(auto &color : mvr.face_color_coords)
  {
    if(color.size() != 3)
    {
      std::cout << __func__
                << "(): found face color with size != 3. Disabling "
                   "face color for all vertices."
                << std::endl;
      mvr.face_color_coords.clear();
      break;
    }
  }

  // check if vertex color is in range 0...255
  bool is_vertex_color_0_255 = false;
  if(! mvr.vertex_color_coords.empty())
  {
    for(auto &color : mvr.vertex_color_coords)
    {
      if(color[0] > 1 || color[1] > 1 || color[2] > 1)
      {
        is_vertex_color_0_255 = true;
        std::cout << __func__
                  << "(): vertex color will be converted from [0,255] "
                     "to [0.0,1.0]."
                  << std::endl;
        break;
      }
    }
  }

  // check if face color is in range 0...255
  bool is_face_color_0_255 = false;
  if(! mvr.face_color_coords.empty())
  {
    for(auto &color : mvr.face_color_coords)
    {
      if(color[0] > 1 || color[1] > 1 || color[2] > 1)
      {
        is_face_color_0_255 = true;
        std::cout << __func__
                  << "(): face color will be converted from [0,255] "
                     "to [0.0,1.0]."
                  << std::endl;
        break;
      }
    }
  }

  // material must be defined for all faces or none
  if(! mvr.face_material.empty())
  {
    if(mvr.face_material.size() != mvr.faces_indices.size())
    {
      std::cout << __func__
                << "(): found some faces with a material and some "
                   "without. Disabling material for all faces."
                << std::endl;
      mvr.face_material.clear();
    }
  }

  // check if all faces have a valid material
  if(! mvr.face_material.empty())
  {
    for(auto &mtl_id : mvr.face_material)
    {
      if(mtl_id == static_cast< index_type >(-1))
      {
        std::cout << __func__
                  << "(): found some faces with an invalid material. "
                     "Disabling material for all faces."
                  << std::endl;
        mvr.face_material.clear();
        break;
      }
    }
  }

  // store materials in a property map
  if(! mvr.materials.empty())
  {
    // load texture images for all materials
    for(size_t i = 0; i < mvr.materials.size(); i++)
    {
      auto &material = mvr.materials[i];

      for(const auto &texture_filename : { material.ambient_texture_filename,
                                           material.diffuse_texture_filename,
                                           material.specular_texture_filename,
                                           material.emissive_texture_filename,
                                           material.transparency_texture_filename,
                                           material.normal_map_filename,
                                           material.metallic_map_filename,
                                           material.roughness_map_filename })
      {
        if(texture_filename.empty())
          continue;

        try
        {
          material.images[texture_filename] =
              std::shared_ptr< FEVV::Types::Image >(
                  new FEVV::Types::Image(texture_filename.c_str()));
        }
        catch(const cimg_library::CImgIOException &/*e*/)
        {
          material.images[texture_filename] =
              std::shared_ptr< FEVV::Types::Image >(
                  make_no_texture_image()
                  );
        }
      }
    }

    // create a property map to store the material
    auto mtl_pm = make_property_map(FEVV::mesh_materials, g);

    // populate material property map
    for(size_t i = 0; i < mvr.materials.size(); i++)
      put(mtl_pm, i, mvr.materials[i]);

    // add property map to bag
    put_property_map(FEVV::mesh_materials, g, pmaps, mtl_pm);

    std::cout << __func__
              << "(): property map for materials created."
              << std::endl;
  }

  // import data into mesh datastructure

  bool use_vertex_normal = false;
  bool use_vertex_color = false;
  bool use_vertex_texture_coord = false;
  bool use_face_color = false;
  bool use_face_material = false;
  bool use_face_normal = false; // per-face vertices normals
  //unused  bool use_line_color = false;
  bool use_corner_texture_coord = params.m_use_corner_texcoord;

  if(! mvr.texture_coords.empty())
  {
    // choose between vertex-texcoord and corner-texcoord
    if(use_corner_texture_coord)
    {
      // corner-texture must be used
      // nothing to do
    }
    else if(mvr.texture_coords.size() == mvr.points_coords.size())
    {
      // texture coordinates are given by vertex
      use_vertex_texture_coord = true;
    }
    else if(mvr.texture_face_indices.size() == mvr.faces_indices.size())
    {
      // texture coordinates are given by corner
      use_corner_texture_coord = true;
    }

    // create property map to store texcoord
    if(use_vertex_texture_coord)
    {
      // create property map to store vertex texture-coord
      auto vtexcoord_pm = make_property_map(FEVV::vertex_texcoord, g);
      put_property_map(FEVV::vertex_texcoord, g, pmaps, vtexcoord_pm);

      std::cout << __func__
                << "(): property map for vertex-texture-coordinate created."
                << std::endl;
    }
    else if(use_corner_texture_coord)
    {
      // create property map to store corner texture-coord
      auto htexcoord_pm = make_property_map(FEVV::halfedge_texcoord, g);
      put_property_map(FEVV::halfedge_texcoord, g, pmaps, htexcoord_pm);

      std::cout << __func__
                << "(): property map for halfedge-texture-coordinate "
                   "created."
                << std::endl;
    }
  }
  else
  {
    // no texture coord provided in mesh vertex representation
    use_corner_texture_coord = false;
  }
  

  if(! mvr.normals_coords.empty())
  {
    // per-vertex normal or per-face-vertex normal ?
    // if there are as many normals as vertices, then it's per-vertex normal

    if(mvr.normals_coords.size() == mvr.points_coords.size())
    {
      use_vertex_normal = true;

      // create property map to store vertex normal
      auto vn_pm = make_property_map(FEVV::vertex_normal, g);
      put_property_map(FEVV::vertex_normal, g, pmaps, vn_pm);

      std::cout << __func__
                << "(): property map for vertex-normal created."
                << std::endl;
    }
    else if(! mvr.normal_face_indices[0].empty())
    {
      // per-face-vertex normal
      use_face_normal = true;

      // create property map to store face-vertex normal
      auto fn_pm = make_property_map(FEVV::face_normal, g);
      put_property_map(FEVV::face_normal, g, pmaps, fn_pm);

      std::cout << __func__
                << "(): property map for face-normal created."
                << std::endl;
    }
  }

  if(! mvr.vertex_color_coords.empty())
  {
    use_vertex_color = true;

    // create property map to store vertex color
    auto vc_pm = make_property_map(FEVV::vertex_color, g);
    put_property_map(FEVV::vertex_color, g, pmaps, vc_pm);

    std::cout << __func__
              << "(): property map for vertex-color created."
              << std::endl;
  }

  if(! mvr.face_color_coords.empty())
  {
    use_face_color = true;

    // create property map to store face color
    auto fc_pm = make_property_map(FEVV::face_color, g);
    put_property_map(FEVV::face_color, g, pmaps, fc_pm);

    std::cout << __func__
              << "(): property map for face-color created."
              << std::endl;
  }

  if(! mvr.face_material.empty())
  {
    use_face_material = true;

    // create property map to store face material
    auto fm_pm = make_property_map(FEVV::face_material, g);
    put_property_map(FEVV::face_material, g, pmaps, fm_pm);

    std::cout << __func__
              << "(): property map for face-material created."
              << std::endl;
  }

  /////////////////////////////// POPULATING MESH
  //////////////////////////////////////


  /////////////////////////////// VERTICES ADDITION
  //////////////////////////////////////

  std::vector< vertex_descriptor > vertices;
  vertices.reserve(mvr.points_coords.size());

  typedef typename std::vector< std::vector< coordP_type > >::const_iterator
      point_it;
  auto point_pm = get(boost::vertex_point, g);
  int vertex_nbr = 0;
  // loop over vertices/points
  for(point_it it_p = mvr.points_coords.begin();
      it_p != mvr.points_coords.end();
      ++it_p)
  {
    vertex_descriptor current_vertex = add_vertex(g);
    vertex_nbr++;

    put(point_pm, current_vertex, Point((*it_p)[0], (*it_p)[1], (*it_p)[2]));

    vertices.push_back(current_vertex);
      // note : must keep the vertex in a vector to handle
      // vertex indices for faces

    // Vertex normal/color/texture-coords addition when present
    if(use_vertex_normal || use_vertex_color || use_vertex_texture_coord)
    {
      if(use_vertex_normal)
      {
        // DBG std::cout << "store normal by vertex #" << vertexNbr <<
        // std::endl;
        std::vector< coordN_type > &coords = mvr.normals_coords[vertex_nbr - 1];

        if(coords.size() != 3)
          throw std::runtime_error("Support of normal coordinates of size != 3 "
                                   "not yet implemented!");
        // store normal in property map
        Vector vn(coords[0], coords[1], coords[2]);
        auto vn_pm = get_property_map(FEVV::vertex_normal, g, pmaps);
        put(vn_pm, current_vertex, vn);
      }
      if(use_vertex_color)
      {
        // DBG std::cout << "store color for vertex #" << vertexNbr <<
        // std::endl;
        std::vector< coordC_type > &coords =
            mvr.vertex_color_coords[vertex_nbr - 1];
        // at this point we are sure that coords.size() == 3 because of
        // previous sanity checks

        // fix colors in [0.0, 1.0] if provided in [0, 255]
        if(is_vertex_color_0_255)
        {
          // conversion copied from Mepp1
          // (src/mepp/Polyhedron/Polyhedron_OFF_CGALImporter.h)
          int rgb[3];
          rgb[0] = (int)floor(coords[0] + 0.5);
          rgb[1] = (int)floor(coords[1] + 0.5);
          rgb[2] = (int)floor(coords[2] + 0.5);
          coords[0] = (coordC_type)((float)rgb[0] / 255.0);
          coords[1] = (coordC_type)((float)rgb[1] / 255.0);
          coords[2] = (coordC_type)((float)rgb[2] / 255.0);
        }

        // store color in property map
        Vector vc(coords[0], coords[1], coords[2]);
        auto vc_pm = get_property_map(FEVV::vertex_color, g, pmaps);
        put(vc_pm, current_vertex, vc);
      }
      if(use_vertex_texture_coord)
      {
        // DBG std::cout << "store tex coord for vertex #" << vertexNbr <<
        // std::endl;
        std::vector< coordT_type > &coords = mvr.texture_coords[vertex_nbr - 1];

        if(coords.size() != 2)
          throw std::runtime_error("Support of texture coordinates of size != "
                                   "2 not yet implemented!");

        // store texture coordinates in property map
        Vector tex_coords(coords[0], coords[1], 0);
        auto vt_pm = get_property_map(FEVV::vertex_texcoord, g, pmaps);
        put(vt_pm, current_vertex, tex_coords);
      }
    }

    // save current vertex descriptor for later use
    if(params.m_vd_target)
    {
      params.m_vd_target->push_back(current_vertex);
    }
  }

  /////////////////////////////// FACES ADDITION
  //////////////////////////////////////

  dup_v_nbr = 0;
  typedef typename std::vector< std::vector< index_type > >::const_iterator
      face_it;
  typedef typename std::vector< index_type >::const_iterator  index_it;
  size_t face_id = 0;
  // loop over the faces
  for(face_it it_f = mvr.faces_indices.begin();
      it_f != mvr.faces_indices.end();
      ++it_f, ++face_id)
  {
    if(it_f->size() < 2)
    {
      throw std::runtime_error(std::string(__func__) +
                               "(): find a face with less than two "
                               "vertices : not conform.");
    }
    else if(it_f->size() == 2)
    {
      throw std::runtime_error(std::string(__func__) +
                               "(): found a face with two vertices, "
                               "isolated edge not supported.");
    }
    else
    {
      //--- create the face

      std::vector< vertex_descriptor > face_vertices;

      // loop over face's vertices
      for(index_it it = it_f->begin(); it != it_f->end(); ++it)
      {
        // add current vertex to vertices list
        face_vertices.push_back(vertices[(*it)]);
      }

      // add face to mesh

      // DBG  static unsigned int counter = 0;
      // DBG  counter ++;
      // DBG  std::cout << "call CGAL::Euler::add_face(face_vertices, g) #" <<
      // counter << std::endl;

      face_descriptor current_face = CGAL::Euler::add_face(face_vertices, g);

      if(current_face == GraphTraits::null_face())
      {
        // helper lambda function to duplicate a vertex and its attributes
        auto duplicate_v = [ &face_vertices,
                             &g,
                             &dup_v_nbr,
                             &pmaps,
                             &use_vertex_normal,
                             &use_vertex_color,
                             &use_vertex_texture_coord ](size_t i) -> void
        {
          vertex_descriptor vold = face_vertices[i];

          // create new vertex
          vertex_descriptor vnew = add_vertex(g);
          face_vertices[i] = vnew;
          dup_v_nbr++;

          // set new vertex geometry
          auto point_pm = get(boost::vertex_point, g);
          put(point_pm, vnew, get(point_pm, vold));

          // set new vertex attributes
          if(use_vertex_normal)
          {
            auto vn_pm = get_property_map(FEVV::vertex_normal, g, pmaps);
            put(vn_pm, vnew, get(vn_pm, vold));
          }
          if(use_vertex_color)
          {
            auto vc_pm = get_property_map(FEVV::vertex_color, g, pmaps);
            put(vc_pm, vnew, get(vc_pm, vold));
          }
          if(use_vertex_texture_coord)
          {
            auto vt_pm = get_property_map(FEVV::vertex_texcoord, g, pmaps);
            put(vt_pm, vnew, get(vt_pm, vold));
          }
        };

        // the code in this section addresses the two
        // first reasons why Euler::add_face(vr, g)
        // returns null_face():
        // - one face's vertex is not on a border
        // - one face's halfedge is not on a border
        // (see Euler_operations.h line 593 to 606)

        // DBG  std::cout << "current_face == GraphTraits::null_face()" <<
        // std::endl;

        // try to add face by duplicating vertices

        // duplicate faulty vertices
        bool face_vertex_duplicated = false;
        size_t n = face_vertices.size();
        for(size_t i = 0, ii = 1; i < n; ++i, ++ii, ii %= n)
        // the loop above and the 'faulty vertex' conditions below
        // are inspired by Euler_operations.h line 593 to 606
        {
          std::pair< halfedge_descriptor, bool > he =
              halfedge(face_vertices[i], face_vertices[ii], g);

          if((! CGAL::internal::is_isolated(face_vertices[i], g) &&
              ! CGAL::is_border(face_vertices[i], g)) ||
             (he.second && ! CGAL::is_border(he.first, g)))
          {
            // the vertex is not isolated and not on a border,
            // or the halfedge from this vertex to the next one
            // is not on a border ; means that the mesh is non manifold ;
            // the vertex must be duplicated to work with
            // non-manifold datastructures

            // DBG  std::cout << "DUPLICATE ONE VERTEX" << std::endl;

            duplicate_v(i);
            face_vertex_duplicated = true;
          }
        }

        // duplicate all face vertices if necessary

        if(! face_vertex_duplicated)
        {
          // no faulty vertex found, let's duplicate all face vertices

          // DBG  std::cout << "DUPLICATE ALL FACE-VERTICES" << std::endl;

          for(size_t i = 0; i < n; ++i)
          {
            duplicate_v(i);
          }
        }

        // add face with duplicated vertices
        current_face = CGAL::Euler::add_face(face_vertices, g);

        if(current_face == GraphTraits::null_face())
        {
          // despite all efforts the face could not be created
          std::cout << "Warning: failed to create one face." << std::endl;
          // go to next face
          continue;
        }
      }
      // the face was successfully created

      //--- fill property maps

      vertex_descriptor prev_vertex = face_vertices.back();
      vertex_descriptor current_vertex;

      // loop over face's vertices
      for(index_type face_vertex_id = 0;
          face_vertex_id < face_vertices.size();
          ++face_vertex_id)
      {
        current_vertex = face_vertices[face_vertex_id];

        // corner (aka halfedge) texture
        if(use_corner_texture_coord)
        {
          // DBG std::cout << "store tex coord for face #" << faceId << " and
          // vertex #" << faceVertexId << std::endl;

          // retrieve texture coords
          Vector texcoords;
          if(! mvr.texture_face_indices[face_id].empty())
          {
            index_type texture_id =
                mvr.texture_face_indices[face_id][face_vertex_id];
            std::vector< coordT_type > &coords = mvr.texture_coords[texture_id];
            if(coords.size() != 2)
              throw std::runtime_error(std::string(__func__) +
                                       "(): texture coordinates of "
                                       "size != 2 not yet supported!");
            texcoords = Vector(coords[0], coords[1], 0);
          }
          else
            texcoords = Vector(0, 0, 0);

          // retrieve halfedge
          // note: the most obvious solution
          //    h = halfedge(prevVertex, currentVertex, g).first
          // doesn't work with AIF in case of a non-manifold halfedge
          // because the (prevVertex, currentVertex) couple is not
          // sufficient to identify the right halfedge (the face is missing)
          halfedge_descriptor h = halfedge(current_face, g);
          while(
              !(source(h, g) == prev_vertex && target(h, g) == current_vertex))
            h = next(h, g);

          // fill in property map
          auto htexcoord_pm =
              get_property_map(FEVV::halfedge_texcoord, g, pmaps);
          put(htexcoord_pm, h, texcoords);

          // DBG std::cout << __FILE__ << ":" << __LINE__ << " he=" << he << "
          // (*pm)[he][0]=" << (*pm)[he][0] << "   (*pm)[he][1]=" << (*pm)[he][1]
          // << "   coords[0]=" << coords[0] << "   coords[1]=" << coords[1] << "
          // texCoords[0]=" << texCoords[0] << "   texCoords[1]=" << texCoords[1]
          // << std::endl;
        }

        // face normal
        // compute face normal incrementally using per-vertex normals
        // (face normal = mean of face's vertices normals)
        if(use_face_normal)
        {
          // DBG std::cout << "compute normal for face #" << faceId << " (vertex
          // #" << faceVertexId << ")" << std::endl;

          // retrieve current vertex normal
          index_type normal_id =
              mvr.normal_face_indices[face_id][face_vertex_id];
          std::vector< coordN_type > &normal = mvr.normals_coords[normal_id];
          if(normal.size() != 3)
            throw std::runtime_error(std::string(__func__) +
                                     "(): normals of size != 3 not yet "
                                     "supported!");
          Vector vnormal = Vector(normal[0], normal[1], normal[2]);

          // retrieve face normal (previous increment)
          auto fn_pm = get_property_map(FEVV::face_normal, g, pmaps);
          Vector fnormal = get(fn_pm, current_face);

          // compute face normal (mean of the vertices normals)
          fnormal = fnormal + vnormal;
          if(face_vertex_id == face_vertices.size() - 1)
            fnormal =
                fnormal / (static_cast< unsigned long >(face_vertex_id) + 1);

          // update face-normal map
          put(fn_pm, current_face, fnormal);
        }

        prev_vertex = current_vertex;
      }

      // set face color
      if(use_face_color)
      {
        std::vector< coordC_type > &color = mvr.face_color_coords[face_id];
        // at this point we are sure that color.size() == 3 because of
        // previous sanity checks

        // fix colors in [0.0, 1.0] if provided in [0, 255]
        if(is_face_color_0_255)
        {
          // conversion copied from Mepp1
          // (src/mepp/Polyhedron/Polyhedron_OFF_CGALImporter.h)
          int rgb[3];
          rgb[0] = (int)floor(color[0] + 0.5);
          rgb[1] = (int)floor(color[1] + 0.5);
          rgb[2] = (int)floor(color[2] + 0.5);
          color[0] = (coordC_type)((float)rgb[0] / 255.0);
          color[1] = (coordC_type)((float)rgb[1] / 255.0);
          color[2] = (coordC_type)((float)rgb[2] / 255.0);
        }

        // store color in property map
        Vector fc(color[0], color[1], color[2]);
        auto fc_pm = get_property_map(FEVV::face_color, g, pmaps);
        put(fc_pm, current_face, fc);
      }

      // set face material
      if(use_face_material)
      {
        size_t mtl_id = mvr.face_material[face_id];

        // store material id in property map
        auto fm_pm = get_property_map(FEVV::face_material, g, pmaps);
        put(fm_pm, current_face, mtl_id);
      }

      // save current face descriptor for later use
      if(params.m_fd_target)
      {
        params.m_fd_target->push_back(current_face);
      }
    }
  }
}


/**
 * \brief  Build the mesh from the given vector representation.
 *
 * \param  g           the output mesh
 * \param  pmaps       the output property maps bag of the mesh
 * \param  dup_v_nbr   number of duplicated vertices of the mesh
 * \param  mvr         the vector representation structure (input)
 * \param  params      the named parameters
 * 
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename HalfedgeGraph,
          typename coordP_type,
          typename coordN_type,
          typename coordT_type,
          typename coordC_type,
          typename index_type,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
mesh_from_vector_representation(
    const HalfedgeGraph                                 &g,
    const FEVV::PMapsContainer                          &pmaps,
    unsigned int                                        &dup_v_nbr,
    FEVV::Types::MVR < coordP_type,
                       coordN_type,
                       coordT_type,
                       coordC_type,
                       index_type >                     &mvr,
    MeshFromVectorReprParameters< HalfedgeGraph > const &params =
        MeshFromVectorReprParameters< HalfedgeGraph >())
{
  GeometryTraits gt(g);
  mesh_from_vector_representation(g, pmaps, dup_v_nbr, mvr, params, gt);
}


} // namespace Filters
} // namespace FEVV
