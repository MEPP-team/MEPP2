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
#include "FEVV/Wrappings/Graph_traits_extension.h"
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"

// DBG  #undef NDEBUG // disable some unwanted assertions in Debug mode


#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/internal/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h> // for add_face(vr,g)


#include "FEVV/Tools/IO/ObjFileReader.h"
#include "FEVV/Tools/IO/OffFileReader.h"
#include "FEVV/Tools/IO/PlyFileReader.h"
#ifdef FEVV_USE_VTK
#include "FEVV/Tools/IO/VtkFileReader.h"
#endif
#ifdef FEVV_USE_FBX
#include "FEVV/Tools/IO/FbxFileReader.h"
#endif
#include "FEVV/Tools/IO/MshFileReader.h"

#include "FEVV/Types/Material.h"

#include <chrono> // for time measurement

namespace FEVV {
namespace Filters {


/**
 *
 * \brief  Measure time since starting time and reset starting time.
 *
 * \param  Starting time
 * \return Elapsed time in seconds since starting time.
 */
inline double
get_time_and_reset(
    std::chrono::time_point< std::chrono::steady_clock > &time_start)
{
    auto time_now = std::chrono::steady_clock::now();
    std::chrono::duration< double > duration = time_now - time_start;
    time_start = time_now;

    return duration.count();
}


/**
 * \brief  Load mesh from file.
 *
 * \param  filename    name of the input mesh file
 * \param  g           mesh object to store the mesh being read
 * \param  pmaps       property maps bag to store the properties
 *                     of the mesh being read
 * \param  gt          the geometry traits to use
 *
 * \sa      the simplified variant that use the default geometry traits
 *          of the mesh.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
read_mesh(const std::string &filename,
          HalfedgeGraph &g,
          PMapsContainer &pmaps,
          const GeometryTraits &gt, bool only_pts/*=false*/)
{
  typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::edge_descriptor edge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_iterator face_iterator;
  typedef typename GeometryTraits::Point Point;
  typedef typename GeometryTraits::Vector Vector;

  // TODO-elo-make templatize this 4 types ? Extract them from the mesh type ?
  typedef double coordP_type; // position coordinate type
  typedef double coordN_type; // normal coordinate type
  typedef float coordC_type;  // color coordinate type
  typedef float coordT_type;  // texture coordinate type
  typedef size_t index_type;


  if(!boost::filesystem::exists(filename))
    throw std::invalid_argument(
        "read_mesh() -> input file path refer to no existing file.");
  else if(!FEVV::FileUtils::has_extension(filename))
    throw std::invalid_argument(
        "read_mesh() -> input file path find without any extension.");

  std::cout << "Generic reading of \""
            << FEVV::FileUtils::get_file_full_name(filename) << "\""
            << std::endl;

  std::vector< std::string > valid_extensions = {
      ".obj", ".off", ".coff", ".ply", ".msh"};
  std::vector< std::string > valid_vtk_extensions = {".vtk", ".vtp", ".vtu"};
  if(!(FEVV::FileUtils::has_extension(filename, valid_extensions)
#ifdef FEVV_USE_VTK
       || FEVV::FileUtils::has_extension(filename, valid_vtk_extensions)
#endif
#ifdef FEVV_USE_FBX
       || FEVV::FileUtils::has_extension(filename, ".fbx")
#endif
           ))
  {
    throw std::invalid_argument(
        "read_mesh() -> input file extension can't be read (yet).");
  }

  ///////////////////////////////////////////////////////////////////////////

  std::vector< std::vector< coordP_type > > points_coords;
  std::vector< std::vector< coordN_type > > normals_coords;
  std::vector< std::vector< coordT_type > > texture_coords;
  std::vector< std::vector< coordC_type > > vertex_color_coords;
  std::vector< std::vector< coordC_type > > face_color_coords;
  std::vector< std::vector< index_type > > lines_indices, faces_indices;
  std::vector< std::vector< index_type > > texture_face_indices;
  // index of the texture coordinates in texture_coords container
  std::vector< std::vector< index_type > > normal_face_indices;
  std::vector< std::vector< coordC_type > > points_colors, faces_colors,
      lines_colors;
  std::vector< std::vector< std::vector< double > > > field_attributes;
  std::vector< std::string > field_names;
  std::vector< FEVV::Types::Material > materials;
  std::vector< index_type > face_material;

  bool obj_file = false;

  auto time = std::chrono::steady_clock::now();

  if(FEVV::FileUtils::has_extension(filename, ".obj"))
  {
    IO::read_obj_file(filename,
                      points_coords,
                      normals_coords,
                      texture_coords,
                      vertex_color_coords,
                      faces_indices,
                      texture_face_indices,
                      normal_face_indices,
                      materials,
                      face_material);
    obj_file = true;
  }
  else if(FEVV::FileUtils::has_extension(filename, ".off") ||
          FEVV::FileUtils::has_extension(filename, ".coff"))
  {
    IO::read_off_file(filename,
                      points_coords,
                      normals_coords,
                      texture_coords,
                      vertex_color_coords,
                      faces_indices,
                      face_color_coords);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
    IO::read_ply_file(filename,
                      points_coords,
                      normals_coords,
                      texture_coords,
                      vertex_color_coords,
                      faces_indices,
                      texture_face_indices);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".msh"))
  {
    IO::read_gmsh_file(filename,
                       points_coords,
                       normals_coords,
                       vertex_color_coords,
                       lines_indices,
                       lines_colors,
                       faces_indices,
                       face_color_coords,
                       field_attributes,
                       field_names);
  }
#ifdef FEVV_USE_VTK
  else if(FEVV::FileUtils::has_extension(filename, ".vtk") ||
          FEVV::FileUtils::has_extension(filename, ".vtp") ||
          FEVV::FileUtils::has_extension(filename, ".vtu"))
  {
    IO::read_vtk_or_vtp_or_vtu_file(filename,
                                    points_coords,
                                    normals_coords,
                                    vertex_color_coords,
                                    lines_indices,
                                    lines_colors,
                                    faces_indices,
                                    face_color_coords,
                                    field_attributes,
                                    field_names);
  }
#endif
#ifdef FEVV_USE_FBX
  else if(FEVV::FileUtils::has_extension(filename, ".fbx"))
  {
    IO::read_fbx_file(filename,
                      points_coords,
                      normals_coords,
                      texture_coords,
                      vertex_color_coords,
                      faces_indices,
                      texture_face_indices,
                      normal_face_indices,
                      materials,
                      face_material);
  }
#endif

  auto parsing_duration = get_time_and_reset(time);

  if (only_pts)
  {
    //normals_coords.clear(); // we can also keep that...
    texture_coords.clear();
    //vertex_color_coords.clear();
    face_color_coords.clear();
    lines_indices.clear(); faces_indices.clear();
    texture_face_indices.clear();
    normal_face_indices.clear();
    points_colors.clear(); faces_colors.clear(); lines_colors.clear();
    field_attributes.clear();
    field_names.clear();
    materials.clear();
    face_material.clear();
  }

  ///////////////////////// CREATE NEEDED PROPERTY MAPS
  //////////////////////////////

  // vertex-color must be defined for all vertices or none
  if(vertex_color_coords.size() != points_coords.size())
    vertex_color_coords.clear();
  for(auto &color : vertex_color_coords)
  {
    if(color.size() != 3)
    {
      std::cout << "read_mesh(): found vertex color with size != 3. Disabling "
                   "vertex color for all vertices."
                << std::endl;
      vertex_color_coords.clear();
      break;
    }
  }

  // face-color must be defined for all faces or none
  if(face_color_coords.size() != faces_indices.size())
    face_color_coords.clear();
  for(auto &color : face_color_coords)
  {
    if(color.size() != 3)
    {
      std::cout << "read_mesh(): found face color with size != 3. Disabling "
                   "face color for all vertices."
                << std::endl;
      face_color_coords.clear();
      break;
    }
  }

  // check if vertex color is in range 0...255
  bool is_vertex_color_0_255 = false;
  if(!vertex_color_coords.empty())
  {
    for(auto &color : vertex_color_coords)
    {
      if(color[0] > 1 || color[1] > 1 || color[2] > 1)
      {
        is_vertex_color_0_255 = true;
        std::cout << "read_mesh(): vertex color will be converted from [0,255] "
                     "to [0.0,1.0]."
                  << std::endl;
        break;
      }
    }
  }

  // check if face color is in range 0...255
  bool is_face_color_0_255 = false;
  if(!face_color_coords.empty())
  {
    for(auto &color : face_color_coords)
    {
      if(color[0] > 1 || color[1] > 1 || color[2] > 1)
      {
        is_face_color_0_255 = true;
        std::cout << "read_mesh(): face color will be converted from [0,255] "
                     "to [0.0,1.0]."
                  << std::endl;
        break;
      }
    }
  }

  // material must be defined for all faces or none
  if(!face_material.empty())
  {
    if(face_material.size() != faces_indices.size())
    {
      std::cout << "read_mesh(): found some faces with a material and some "
                   "without. Disabling material for all faces."
                << std::endl;
      face_material.clear();
    }
  }

  // check if all faces have a valid material
  if(!face_material.empty())
  {
    for(auto &mtl_id : face_material)
    {
      if(mtl_id == -1)
      {
        std::cout << "read_mesh(): found some faces with an invalid material. "
                     "Disabling material for all faces."
                  << std::endl;
        face_material.clear();
        break;
      }
    }
  }

  // store materials in a property map
  if(!materials.empty())
  {
    // create a property map to store the material
    auto mtl_pm = make_property_map(FEVV::mesh_materials, g);

    // populate material property map
    for(size_t i = 0; i < materials.size(); i++)
      put(mtl_pm, i, materials[i]);

    // add property map to bag
    put_property_map(FEVV::mesh_materials, g, pmaps, mtl_pm);

    std::cout << "read_mesh(): property map for materials created."
              << std::endl;
  }

  // import data into mesh datastructure

  bool use_vertex_normal = false;
  bool use_vertex_color = false;
  bool use_vertex_texture_coord = false;
  bool use_corner_texture_coord = false;
  bool use_face_color = false;
  bool use_face_material = false;
  bool use_face_normal = false; // per-face vertices normals
  bool use_line_color = false;

  if(!texture_coords.empty())
  {
    // choose between vertex-texcoord and corner-texcoord
    if(obj_file)
    {
      // support corner-texture only with OBJ file
      use_corner_texture_coord = true;
    }
    else if(texture_coords.size() == points_coords.size())
    {
      // texture coordinates are given by vertex
      use_vertex_texture_coord = true;
    }
    else if(texture_face_indices.size() == faces_indices.size())
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

      std::cout
          << "read_mesh(): property map for vertex-texture-coordinate created."
          << std::endl;
    }
    else if(use_corner_texture_coord)
    {
      // create property map to store corner texture-coord
      auto htexcoord_pm = make_property_map(FEVV::halfedge_texcoord, g);
      put_property_map(FEVV::halfedge_texcoord, g, pmaps, htexcoord_pm);

      std::cout << "read_mesh(): property map for halfedge-texture-coordinate "
                   "created."
                << std::endl;
    }
  }

  if(!normals_coords.empty())
  {
    // per-vertex normal or per-face-vertex normal ?
    // if there are as many normals as vertices, then it's per-vertex normal

    if(normals_coords.size() == points_coords.size())
    {
      use_vertex_normal = true;

      // create property map to store vertex normal
      auto vn_pm = make_property_map(FEVV::vertex_normal, g);
      put_property_map(FEVV::vertex_normal, g, pmaps, vn_pm);

      std::cout << "read_mesh(): property map for vertex-normal created."
                << std::endl;
    }
    else if(!normal_face_indices[0].empty())
    {
      // per-face-vertex normal
      use_face_normal = true;

      // create property map to store face-vertex normal
      auto fn_pm = make_property_map(FEVV::face_normal, g);
      put_property_map(FEVV::face_normal, g, pmaps, fn_pm);

      std::cout << "read_mesh(): property map for face-normal created."
                << std::endl;
    }
  }

  if(!vertex_color_coords.empty())
  {
    use_vertex_color = true;

    // create property map to store vertex color
    auto vc_pm = make_property_map(FEVV::vertex_color, g);
    put_property_map(FEVV::vertex_color, g, pmaps, vc_pm);

    std::cout << "read_mesh(): property map for vertex-color created."
              << std::endl;
  }

  if(!face_color_coords.empty())
  {
    use_face_color = true;

    // create property map to store face color
    auto fc_pm = make_property_map(FEVV::face_color, g);
    put_property_map(FEVV::face_color, g, pmaps, fc_pm);

    std::cout << "read_mesh(): property map for face-color created."
              << std::endl;
  }

  if(!face_material.empty())
  {
    use_face_material = true;

    // create property map to store face material
    auto fm_pm = make_property_map(FEVV::face_material, g);
    put_property_map(FEVV::face_material, g, pmaps, fm_pm);

    std::cout << "read_mesh(): property map for face-material created."
              << std::endl;
  }

  /////////////////////////////// POPULATING MESH
  //////////////////////////////////////


  /////////////////////////////// VERTICES ADDITION
  //////////////////////////////////////

  std::vector< vertex_descriptor > vertices;
  vertices.reserve(points_coords.size());

  typedef std::vector< std::vector< coordP_type > >::const_iterator point_it;
  auto point_pm = get(boost::vertex_point, g);
  int vertex_nbr = 0;
  // loop over vertices/points
  for(point_it it_p = points_coords.begin(); it_p != points_coords.end();
      ++it_p)
  {
    vertex_descriptor current_vertex = add_vertex(g);
    vertex_nbr++;

    put(point_pm, current_vertex, Point((*it_p)[0], (*it_p)[1], (*it_p)[2]));

    vertices.push_back(
        current_vertex); // note : must keep the vertex in a vector to handle
                         // vertex indices for faces

    // Vertex normal/color/texture-coords addition when present
    if(use_vertex_normal || use_vertex_color || use_vertex_texture_coord)
    {
      if(use_vertex_normal)
      {
        // DBG std::cout << "store normal by vertex #" << vertexNbr <<
        // std::endl;
        std::vector< coordN_type > &coords = normals_coords[vertex_nbr - 1];

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
            vertex_color_coords[vertex_nbr - 1];
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
        std::vector< coordT_type > &coords = texture_coords[vertex_nbr - 1];

        if(coords.size() != 2)
          throw std::runtime_error("Support of texture coordinates of size != "
                                   "2 not yet implemented!");

        // store texture coordinates in property map
        Vector tex_coords(coords[0], coords[1], 0);
        auto vt_pm = get_property_map(FEVV::vertex_texcoord, g, pmaps);
        put(vt_pm, current_vertex, tex_coords);
      }
    }
  }

  /////////////////////////////// FACES ADDITION
  //////////////////////////////////////

  unsigned int duplicated_vertices_nbr = 0;
  typedef std::vector< std::vector< index_type > >::const_iterator face_it;
  typedef std::vector< index_type >::const_iterator index_it;
  size_t face_id = 0;
  // loop over the faces
  for(face_it it_f = faces_indices.begin(); it_f != faces_indices.end();
      ++it_f, ++face_id)
  {
    if(it_f->size() < 2)
    {
      throw std::runtime_error("read_mesh(): find a face with less than two "
                               "vertices : not conform.");
    }
    else if(it_f->size() == 2)
    {
      throw std::runtime_error("read_mesh(): found a face with two vertices, "
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
        size_t n = face_vertices.size();
        for(size_t i = 0, ii = 1; i < n; ++i, ++ii, ii %= n)
        // the loop above and the 'faulty vertex' conditions below
        // are inspired by Euler_operations.h line 593 to 606
        {
          std::pair< halfedge_descriptor, bool > he =
              halfedge(face_vertices[i], face_vertices[ii], g);

          if((!CGAL::internal::is_isolated(face_vertices[i], g) &&
              !CGAL::is_border(face_vertices[i], g)) ||
             (he.second && !CGAL::is_border(he.first, g)))
          {
            // the vertex is not isolated and not on a border,
            // or the halfedge from this vertex to the next one
            // is not on a border ; means that the mesh is non manifold ;
            // the vertex must be duplicated to work with
            // non-manifold datastructures

            // DBG  std::cout << "DUPLICATE A VERTEX" << std::endl;

            vertex_descriptor vold = face_vertices[i];

            // create new vertex
            vertex_descriptor vnew = add_vertex(g);
            face_vertices[i] = vnew;
            duplicated_vertices_nbr++;

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
          }
        }

        // add face with duplicated vertices
        current_face = CGAL::Euler::add_face(face_vertices, g);
      }
      assert(current_face != GraphTraits::null_face());

      //--- fill property maps

      vertex_descriptor prev_vertex = face_vertices.back();
      vertex_descriptor current_vertex;

      // loop over face's vertices
      for(index_type face_vertex_id = 0; face_vertex_id < face_vertices.size();
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
          if(!texture_face_indices[face_id].empty())
          {
            index_type texture_id =
                texture_face_indices[face_id][face_vertex_id];
            std::vector< coordT_type > &coords = texture_coords[texture_id];
            if(coords.size() != 2)
              throw std::runtime_error("read_mesh(): texture coordinates of "
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
          index_type normal_id = normal_face_indices[face_id][face_vertex_id];
          std::vector< coordN_type > &normal = normals_coords[normal_id];
          if(normal.size() != 3)
            throw std::runtime_error(
                "read_mesh(): normals of size != 3 not yet supported!");
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
        std::vector< coordC_type > &color = face_color_coords[face_id];
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
        size_t mtl_id = face_material[face_id];

        // store material id in property map
        auto fm_pm = get_property_map(FEVV::face_material, g, pmaps);
        put(fm_pm, current_face, mtl_id);
      }
    }
  }

  auto populating_duration = get_time_and_reset(time);

  // display some information

  if(duplicated_vertices_nbr > 0)
  {
    std::cout << "read_mesh(): " << duplicated_vertices_nbr
              << " vertices duplicated (non-manifold mesh)." << std::endl;
  }
  std::cout << "read_mesh(): the mesh has " << size_of_vertices(g)
            << " vertices, " << size_of_edges(g) << " edges, "
            << size_of_faces(g) << " faces." << std::endl;
  std::cout << "Generic reading of \""
            << FEVV::FileUtils::get_file_full_name(filename) << "\" done"
            << " in " << parsing_duration + populating_duration << "s"
            << " (parsing " << parsing_duration << "s"
            << ", populating DS " << populating_duration << "s)."
            << std::endl;
}


/**
 * \brief  Load mesh from file.
 *
 * \param  filename    name of the input mesh file
 * \param  g           mesh object to store the mesh being read
 * \param  pmaps       property maps bag to store the properties
 *                     of the mesh being read
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
read_mesh(const std::string &filename, HalfedgeGraph &g, PMapsContainer &pmaps, bool only_pts=false)
{
  GeometryTraits gt(g);
  read_mesh< HalfedgeGraph, GeometryTraits >(filename, g, pmaps, gt, only_pts);
}

} // namespace Filters
} // namespace FEVV
