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

#include <iostream>
#include <vector>
#include <cstdio>
#include <string>
#include <exception>

#include <boost/filesystem.hpp>

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "FEVV/Tools/IO/StringUtilities.hpp"
#include "FEVV/Types/Material.h"

namespace FEVV {
namespace IO {

using namespace StrUtils;
using namespace FileUtils;

/**
 * Write materials to MTL file.
 */
template< typename MaterialType >
void
write_mtl_file(const std::string &mtl_file_name,
               std::vector< MaterialType > &materials)
{
  // create MTL file
  std::ofstream mtl_file(mtl_file_name, std::ios::out | std::ios::trunc);
  if(!mtl_file.is_open())
    throw std::runtime_error("write_MTL_file(): failed to create MTL file.");

  std::string parent_directory = get_parent_directory(mtl_file_name);
  if(parent_directory.empty())
    parent_directory = ".";
  const std::string texture_filepath_base = parent_directory + '/';

  // populate MTL file
  mtl_file << "# MTL File produced by MEPP2" << std::endl;
  for(auto &material : materials)
  {
    mtl_file << "newmtl " << material.name << std::endl;

    // Writing material coefficients
    mtl_file << "    Kd " << material.diffuse_red_component << ' '
             << material.diffuse_green_component << ' '
             << material.diffuse_blue_component << std::endl;

    if(material.type == Types::MaterialType::MATERIAL_TYPE_PBR)
    {
      mtl_file << "    Pm " << material.metallic_factor << std::endl;
      mtl_file << "    Pr " << material.roughness_factor << std::endl;
    }
    else
    {
      mtl_file << "    Ka " << material.ambient_red_component << ' '
               << material.ambient_green_component << ' '
               << material.ambient_blue_component << std::endl;
      mtl_file << "    Ks " << material.specular_red_component << ' '
               << material.specular_green_component << ' '
               << material.specular_blue_component << std::endl;
      mtl_file << "    Ke " << material.emissive_red_component << ' '
               << material.emissive_green_component << ' '
               << material.emissive_blue_component << std::endl;
      mtl_file << "    d " << material.transparency << std::endl;
    }

    // Diffuse/albedo map
    if(!material.diffuse_texture_filename.empty())
    {
      const std::string texture_filename =
          get_file_full_name(material.diffuse_texture_filename);
      mtl_file << "    map_Kd " << texture_filename << std::endl;

      copy_file(material.diffuse_texture_filename,
                texture_filepath_base + texture_filename);
    }

    // Ambient/ambient occlusion map
    if(!material.ambient_texture_filename.empty())
    {
      const std::string texture_filename =
          get_file_full_name(material.ambient_texture_filename);
      mtl_file << "    map_Ka " << texture_filename << std::endl;

      copy_file(material.ambient_texture_filename,
                texture_filepath_base + texture_filename);
    }

    // Normal/bump map
    if(!material.normal_map_filename.empty())
    {
      const std::string texture_filename =
          get_file_full_name(material.normal_map_filename);
      mtl_file << "    "
               << (material.type == Types::MaterialType::MATERIAL_TYPE_PBR
                       ? "norm"
                       : "bump")
               << ' ';
      mtl_file << texture_filename << std::endl;

      copy_file(material.normal_map_filename,
                texture_filepath_base + texture_filename);
    }

    if(material.type == Types::MaterialType::MATERIAL_TYPE_PBR)
    {
      // Metallic map
      if(!material.metallic_map_filename.empty())
      {
        const std::string texture_filename =
            get_file_full_name(material.metallic_map_filename);
        mtl_file << "    map_Pm " << texture_filename << std::endl;

        copy_file(material.metallic_map_filename,
                  texture_filepath_base + texture_filename);
      }

      // Roughness map
      if(!material.roughness_map_filename.empty())
      {
        const std::string texture_filename =
            get_file_full_name(material.roughness_map_filename);
        mtl_file << "    map_Pr " << texture_filename << std::endl;

        copy_file(material.roughness_map_filename,
                  texture_filepath_base + texture_filename);
      }
    }
    else
    {
      // Specular map
      if(!material.specular_texture_filename.empty())
      {
        const std::string texture_filename =
            get_file_full_name(material.specular_texture_filename);
        mtl_file << "    map_Ks " << texture_filename << std::endl;

        copy_file(material.specular_texture_filename,
                  texture_filepath_base + texture_filename);
      }

      // Transparency map
      if(!material.transparency_texture_filename.empty())
      {
        const std::string texture_filename =
            get_file_full_name(material.transparency_texture_filename);
        mtl_file << "    map_d " << texture_filename << std::endl;

        copy_file(material.transparency_texture_filename,
                  texture_filepath_base + texture_filename);
      }
    }

    mtl_file << std::endl;
  }

  mtl_file.close();
}

/**
 * Write mesh data to OBJ file.
 */
template< typename CoordType,
          typename CoordNType,
          typename CoordTType,
          typename CoordCType,
          typename IndexType,
          typename MaterialType >
void
write_obj_file(const std::string &file_path,
               std::vector< std::vector< CoordType > > &points_coords,
               std::vector< std::vector< CoordNType > > &normals_coords,
               std::vector< std::vector< CoordTType > > &texture_coords,
               std::vector< std::vector< CoordCType > > &vertex_color_coords,
               std::vector< std::vector< IndexType > > &face_indices,
               std::vector< std::vector< IndexType > > &texture_face_indices,
               std::vector< std::vector< IndexType > > &normal_face_indices,
               std::vector< MaterialType > &materials,  // list of materials
               std::vector< IndexType > &face_material) // material of each face
{
  size_t nb_points = points_coords.size(), nb_normals = normals_coords.size(),
         nb_tex_coords = texture_coords.size(),
         nb_point_color = vertex_color_coords.size(),
         nb_faces = face_indices.size(),
         nb_texture_ind = texture_face_indices.size(),
         nb_normal_ind = normal_face_indices.size(),
         nb_face_material = face_material.size(),
         nb_materials = materials.size();

  bool use_normal = (nb_normals != 0), use_texture = (nb_tex_coords != 0),
       use_color = (nb_point_color != 0),
       use_face_material = (nb_face_material != 0),
       use_materials = (nb_materials != 0);

  if((nb_point_color != nb_points) && (nb_point_color != 0))
    throw std::invalid_argument(
        "Writer::write_obj_file -> invalid number of vertex colors.");

  if((nb_texture_ind != nb_faces) && (nb_texture_ind != 0))
    throw std::invalid_argument(
        "Writer::write_obj_file -> invalid number of texture indices.");

  if((nb_normal_ind != nb_faces) && (nb_normal_ind != 0))
    throw std::invalid_argument(
        "Writer::write_obj_file -> invalid number of normal indices.");

  if((nb_face_material != nb_faces) && (nb_face_material != 0))
    throw std::invalid_argument(
        "Writer::write_obj_file -> invalid number of face materials.");

  // WRITE MATERIALS TO MTL FILE

  std::string mtl_file_name;

  if(use_materials)
  {
    std::string parent_directory = get_parent_directory(file_path);
    if(parent_directory.empty())
      parent_directory = ".";
    mtl_file_name = get_file_name(file_path) + ".mtl";
    std::string mtl_file_path = parent_directory + "/" + mtl_file_name;
    write_mtl_file(mtl_file_path, materials);
  }

  // WRITE MESH TO OBJ FILE

  std::ofstream file(file_path, std::ios::out | std::ios::trunc);

  if(file.is_open())
  {
    file << "# OBJ File produced by MEPP2" << std::endl
         << "# The M2Disco Framework" << std::endl
         << "#" << std::endl;
    file << "# Mesh with :" << std::endl;
    file << "# \t" << nb_points << " vertex positions" << std::endl;
    file << "# \t" << nb_faces << " faces" << std::endl;
    file << "# \t" << nb_tex_coords << " UV coordinates" << std::endl;
    file << "# \t" << nb_normals << " vertex normals" << std::endl << std::endl;

    // WRITE MTL FILENAME
    if(use_materials)
    {
      // store MTL file name in OBJ file
      file << "mtllib " << mtl_file_name << std::endl;
    }

    // VERTICES WRITING
    typename std::vector< std::vector< CoordCType > >::const_iterator it_c =
        vertex_color_coords.cbegin();
    for(typename std::vector< std::vector< CoordType > >::const_iterator it_p =
            points_coords.cbegin();
        it_p != points_coords.cend();
        ++it_p)
    {
      file << "v ";
      for(typename std::vector< CoordType >::const_iterator it_pc =
              it_p->cbegin();
          it_pc != it_p->cend();
          ++it_pc)
        file << *it_pc << " ";

      if(use_color)
      {
        for(typename std::vector< CoordCType >::const_iterator it_cc =
                it_c->cbegin();
            it_cc != it_c->cend();
            ++it_cc)
          file << *it_cc << " ";
        ++it_c;
      }

      file << std::endl;
    }

    file << std::endl;

    // NORMALS WRITING
    for(typename std::vector< std::vector< CoordNType > >::const_iterator it_n =
            normals_coords.cbegin();
        it_n != normals_coords.cend();
        ++it_n)
    {
      file << "vn ";
      for(typename std::vector< CoordNType >::const_iterator it_nc =
              it_n->cbegin();
          it_nc != it_n->cend();
          ++it_nc)
        file << *it_nc << " ";
      file << std::endl;
    }

    // TEXTURE COORDS WRITING
    for(typename std::vector< std::vector< CoordTType > >::const_iterator it_t =
            texture_coords.cbegin();
        it_t != texture_coords.cend();
        ++it_t)
    {
      file << "vt ";
      for(typename std::vector< CoordTType >::const_iterator it_tc =
              it_t->cbegin();
          it_tc != it_t->end();
          ++it_tc)
        file << *it_tc << " ";
      file << std::endl;
    }

    // FACES WRITING
    IndexType current_material_id = -1;
    std::string current_material_name;
    for(unsigned int no_face = 0; no_face < nb_faces; ++no_face)
    {
      // write a 'usemtl' line when material changes
      if(use_face_material && face_material[no_face] != current_material_id)
      {
        current_material_id = face_material[no_face];
        file << "usemtl " << materials[current_material_id].name << std::endl;
      }

      // write face line
      file << "f ";
      for(unsigned int no_vertex = 0; no_vertex < face_indices[no_face].size();
          ++no_vertex)
      {
        file << (face_indices[no_face][no_vertex] +
                 1); // +1 because obj format starts the indices to 1 (not 0)

        if(use_texture || use_normal)
          file << "/";

        if(use_texture)
          file << (texture_face_indices[no_face][no_vertex] + 1);

        if(use_normal)
          file << "/" << (normal_face_indices[no_face][no_vertex] + 1);

        file << " ";
      }
      file << std::endl;
    }

    file.close();
  }
  else
  {
    throw std::runtime_error(
        "Writer::write_obj_file -> output file failed to open.");
  }
}

} // namespace IO
} // namespace FEVV

