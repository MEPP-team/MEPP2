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
#include "FEVV/Wrappings/Graph_traits_extension.h" // for size_of_edges()
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"

#include "FEVV/Types/Mesh_vector_representation.h"
#include "FEVV/Filters/Generic/mesh_to_vector_representation.hpp"

#include "FEVV/Tools/IO/ObjFileWriter.h"
#include "FEVV/Tools/IO/OffFileWriter.h"
#ifdef FEVV_USE_VTK
#include "FEVV/Tools/IO/VtkFileWriter.h"
#endif
#include "FEVV/Tools/IO/MshFileWriter.h"


namespace FEVV {
namespace Filters {


/**
 * \brief  Write mesh to file.
 *
 * \param  filename    name of the output mesh file
 * \param  g           mesh to write to file
 * \param  pmaps       property maps bag of the mesh
 * \param  gt          the geometry traits to use
 *
 * \sa     the simplified variant that use the default geometry traits
 *         of the mesh.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits >
void
write_mesh(const std::string &filename,
           HalfedgeGraph &g,
           PMapsContainer &pmaps,
           const GeometryTraits &gt)
{
  // TODO-elo  templatize this 4 types ? Extract them from the mesh type ?
  typedef double coordP_type; // position coordinate type
  typedef double coordN_type; // normal coordinate type
  typedef float coordC_type;  // color coordinate type
  typedef float coordT_type;  // texture coordinate type
  typedef size_t index_type;


  if(!FEVV::FileUtils::has_extension(filename))
    throw std::invalid_argument(
        "write_mesh(): output file path has no extension.");

  std::cout << "Generic writing of \""
            << FEVV::FileUtils::get_file_full_name(filename) << "\""
            << std::endl;

  std::vector< std::string > valid_extensions = {".obj", ".off", ".coff",
                                                 /*".ply",*/ ".msh"};
  std::vector< std::string > valid_vtk_extensions = {".vtk", ".vtp", ".vtu"};
  if(!(FEVV::FileUtils::has_extension(filename, valid_extensions)
#ifdef FEVV_USE_VTK
       || FEVV::FileUtils::has_extension(filename, valid_vtk_extensions)
#endif
           ))
  {
    throw std::invalid_argument(
        "write_mesh(): output file format not yet supported.");
  }

  ///////////////////////////////////////////////////////////////////////////

  // build the vector representation of the mesh

  FEVV::Types::MVR< coordP_type,
                    coordN_type,
                    coordT_type,
                    coordC_type,
                    index_type > mvr;

  mesh_to_vector_representation(g, pmaps, mvr, gt);

  // WRITE TO FILE

  if(FEVV::FileUtils::has_extension(filename, ".obj"))
  {
    IO::write_obj_file(filename,
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
  else if(FEVV::FileUtils::has_extension(filename, ".off") ||
          FEVV::FileUtils::has_extension(filename, ".coff"))
  {
    IO::write_off_file(filename,
                       mvr.points_coords,
                       mvr.normals_coords,
                       mvr.texture_coords,
                       mvr.vertex_color_coords,
                       mvr.faces_indices,
                       mvr.face_color_coords,
                       static_cast< unsigned long >(size_of_edges(g)));
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
    throw std::runtime_error("write_mesh(): ply writer not yet implemented");
  }
  else if(FEVV::FileUtils::has_extension(filename, ".msh"))
  {
    IO::write_gmsh_file(filename,
                        mvr.points_coords,
                        mvr.normals_coords,
                        mvr.points_colors,
                        mvr.lines_indices,
                        mvr.lines_colors,
                        mvr.faces_indices,
                        mvr.faces_colors,
                        mvr.field_attributes,
                        mvr.field_names);
  }
#ifdef FEVV_USE_VTK
  else if(FEVV::FileUtils::has_extension(filename, ".vtk") ||
          FEVV::FileUtils::has_extension(filename, ".vtp") ||
          FEVV::FileUtils::has_extension(filename, ".vtu"))
  {
    IO::write_vtk_or_vtp_or_vtu_file(filename,
                                     mvr.points_coords,
                                     mvr.normals_coords,
                                     mvr.points_colors,
                                     mvr.lines_indices,
                                     mvr.lines_colors,
                                     mvr.faces_indices,
                                     mvr.faces_colors,
                                     mvr.field_attributes,
                                     mvr.field_names);
  }
#endif

  std::cout << "Generic writing of \""
            << FEVV::FileUtils::get_file_full_name(filename) << "\" Done."
            << std::endl;
}


/**
 * \brief  Write mesh to file.
 *
 * \param  filename    name of the output mesh file
 * \param  g           mesh to write to file
 * \param  pmaps       property maps bag of the mesh
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template< typename HalfedgeGraph,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
write_mesh(const std::string &filename, HalfedgeGraph &g, PMapsContainer &pmaps)
{
  GeometryTraits gt(g);
  write_mesh< HalfedgeGraph, GeometryTraits >(filename, g, pmaps, gt);
}

} // namespace Filters
} // namespace FEVV
