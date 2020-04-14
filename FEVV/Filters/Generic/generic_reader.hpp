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

#include "FEVV/Types/Mesh_vector_representation.h"
#include "FEVV/Filters/Generic/mesh_from_vector_representation.hpp"

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
  // TODO-elo-make templatize this 4 types ? Extract them from the mesh type ?
  typedef double coordP_type; // position coordinate type
  typedef double coordN_type; // normal coordinate type
  typedef float coordC_type;  // color coordinate type
  typedef float coordT_type;  // texture coordinate type
  typedef size_t index_type;

  ///////////////////////////////////////////////////////////////////////////
  // check that input file exists and format is supported

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
  // parse input file

  FEVV::Types::MVR< coordP_type,
                    coordN_type,
                    coordT_type,
                    coordC_type,
                    index_type > mvr;

  bool obj_file = false;

  auto time = std::chrono::steady_clock::now();

  if(FEVV::FileUtils::has_extension(filename, ".obj"))
  {
    IO::read_obj_file(filename,
                      mvr.points_coords,
                      mvr.normals_coords,
                      mvr.texture_coords,
                      mvr.vertex_color_coords,
                      mvr.faces_indices,
                      mvr.texture_face_indices,
                      mvr.normal_face_indices,
                      mvr.materials,
                      mvr.face_material);
    obj_file = true;
  }
  else if(FEVV::FileUtils::has_extension(filename, ".off") ||
          FEVV::FileUtils::has_extension(filename, ".coff"))
  {
    IO::read_off_file(filename,
                      mvr.points_coords,
                      mvr.normals_coords,
                      mvr.texture_coords,
                      mvr.vertex_color_coords,
                      mvr.faces_indices,
                      mvr.face_color_coords);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
    IO::read_ply_file(filename,
                      mvr.points_coords,
                      mvr.normals_coords,
                      mvr.texture_coords,
                      mvr.vertex_color_coords,
                      mvr.faces_indices,
                      mvr.texture_face_indices);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".msh"))
  {
    IO::read_gmsh_file(filename,
                       mvr.points_coords,
                       mvr.normals_coords,
                       mvr.vertex_color_coords,
                       mvr.lines_indices,
                       mvr.lines_colors,
                       mvr.faces_indices,
                       mvr.face_color_coords,
                       mvr.field_attributes,
                       mvr.field_names);
  }
#ifdef FEVV_USE_VTK
  else if(FEVV::FileUtils::has_extension(filename, ".vtk") ||
          FEVV::FileUtils::has_extension(filename, ".vtp") ||
          FEVV::FileUtils::has_extension(filename, ".vtu"))
  {
    IO::read_vtk_or_vtp_or_vtu_file(filename,
                                    mvr.points_coords,
                                    mvr.normals_coords,
                                    mvr.vertex_color_coords,
                                    mvr.lines_indices,
                                    mvr.lines_colors,
                                    mvr.faces_indices,
                                    mvr.face_color_coords,
                                    mvr.field_attributes,
                                    mvr.field_names);
  }
#endif
#ifdef FEVV_USE_FBX
  else if(FEVV::FileUtils::has_extension(filename, ".fbx"))
  {
    IO::read_fbx_file(filename,
                      mvr.points_coords,
                      mvr.normals_coords,
                      mvr.texture_coords,
                      mvr.vertex_color_coords,
                      mvr.faces_indices,
                      mvr.texture_face_indices,
                      mvr.normal_face_indices,
                      mvr.materials,
                      mvr.face_material);
  }
#endif

  auto parsing_duration = get_time_and_reset(time);

  ///////////////////////////////////////////////////////////////////////////
  // populate mesh object

  if (only_pts)
  {
    //normals_coords.clear(); // we can also keep that...
    mvr.texture_coords.clear();
    //vertex_color_coords.clear();
    mvr.face_color_coords.clear();
    mvr.lines_indices.clear();
    mvr.faces_indices.clear();
    mvr.texture_face_indices.clear();
    mvr.normal_face_indices.clear();
    mvr.points_colors.clear();
    mvr.faces_colors.clear();
    mvr.lines_colors.clear();
    mvr.field_attributes.clear();
    mvr.field_names.clear();
    mvr.materials.clear();
    mvr.face_material.clear();
  }

  unsigned int duplicated_vertices_nbr;
  mesh_from_vector_representation(g, pmaps, duplicated_vertices_nbr,
                                  mvr, obj_file, gt);
    // corner texture are supported with OBJ file only

  auto populating_duration = get_time_and_reset(time);

  ///////////////////////////////////////////////////////////////////////////
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
