#pragma once

/*
 * OBJ ref:
 * - http://www.martinreddy.net/gfx/3d/OBJ.spec
 * - https://en.wikipedia.org/wiki/Wavefront_.obj_file
 * - http://paulbourke.net/dataformats/obj
 * - http://paulbourke.net/dataformats/mtl
 */

#include <iostream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <string>
#include <exception>

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "FEVV/Tools/IO/StringUtilities.hpp"
#include "FEVV/Types/Material.h"


namespace FEVV {
namespace IO {


using namespace StrUtils;
using namespace FileUtils;

/**
 * Print information about the materials.
 */
template< typename MaterialType >
void
dbg_display_materials(const std::vector< MaterialType > &materials)
{
  std::cout << "materials:" << std::endl;
  for(auto mat : materials)
  {
    std::cout << "  " << mat.name << std::endl;
    std::cout << "\tdiffuse/albedo: " << mat.diffuse_texture_filename
              << std::endl;
    std::cout << "\tspecular: " << mat.specular_texture_filename << std::endl;
    std::cout << "\ttransparency: " << mat.transparency_texture_filename
              << std::endl;
    std::cout << "\tbump/normal: " << mat.normal_map_filename << std::endl;
    std::cout << "\tmetallic: " << mat.metallic_map_filename << std::endl;
    std::cout << "\troughness: " << mat.roughness_map_filename << std::endl;
    std::cout << "\temissive: " << mat.emissive_texture_filename << std::endl;
    std::cout << "\tambient/ambient occlusion: " << mat.ambient_texture_filename
              << std::endl;
  }
}


/**
 * Read the provided MTL file and populate the materials datastructure.
 */
template< typename MaterialType >
void
read_mtl_file(const std::string &mtl_file_name,
              std::vector< MaterialType > &materials)
{
  // DBG std::cout << "mtl_file_name = " << mtlFileName << std::endl;

  std::ifstream file(mtl_file_name);
  if(!file.is_open())
  {
    std::cout << "Failed to open MTL file '" << mtl_file_name << "'"
              << std::endl;
    return;
  }

  const std::string parent_directory = get_parent_directory(mtl_file_name);

  // read materials
  Types::Material new_material;

  std::string word;

  while(!file.eof())
  {
    file >> word;
    // DBG std::cout << "word=" << word << std::endl;

    // RM: import all material's features (coefficients & textures names)
    if(word == "newmtl")
    {
      if(!new_material.name.empty())
      {
        materials.push_back(new_material);
        new_material = Types::Material();
      }

      file >> new_material.name;
    }
    else if(word[0] == 'm') // Texture [map_*]
    {
      // we suppose there is only one 'map_Kd' per material
      std::string texture_filename;
      file >> texture_filename;
      texture_filename.insert(0, parent_directory + "/");

      if(word[4] == 'K') // Standard textures [map_K*]
      {
        if(word[5] ==
           'a') // Ambient texture [map_Ka]; used for PBR as ambient occlusion
        {
          new_material.ambient_texture_filename = texture_filename;
        }
        else if(word[5] ==
                'd') // Diffuse texture [map_Kd]; used for PBR as albedo
        {
          new_material.diffuse_texture_filename = texture_filename;
        }
        else if(word[5] == 's') // Specular texture [map_Ks]
        {
          new_material.specular_texture_filename = texture_filename;
        }
        else if(word[5] == 'e') // Emissive texture [map_Ke]
        {
          new_material.emissive_texture_filename = texture_filename;
        }
      }
      else if(word[4] == 'P') // PBR texture [map_P*]
      {
        if(word[5] == 'm') // Metallic map [map_Pm]
        {
          new_material.metallic_map_filename = texture_filename;
        }
        else if(word[5] == 'r') // Roughness map [map_Pr]
        {
          new_material.roughness_map_filename = texture_filename;
        }

        new_material.type = Types::MaterialType::MATERIAL_TYPE_PBR;
      }
      else if(word[4] == 'd') // Transparency texture [map_d]
      {
        new_material.transparency_texture_filename = texture_filename;
      }
      else if(word[4] == 'b') // Bump or normal map [map_bump]
      {
        new_material.normal_map_filename = texture_filename;

        // RM: if normal map filename contains "normal" or "nrm", consider it a
        // normal map instead of a bump map
        const std::string lower_tex_filename =
            boost::algorithm::to_lower_copy(texture_filename);
        if(lower_tex_filename.find("normal") != std::string::npos ||
           lower_tex_filename.find("nrm") != std::string::npos)
          new_material.has_normal_map = true;
      }
    }
    else if(word[0] == 'K')
    {
      std::string token_r, token_g, token_b;
      file >> token_r >> token_g >> token_b;

      if(word[1] == 'a') // Ambient color [Ka]
      {
        convert(token_r, new_material.ambient_red_component);
        convert(token_g, new_material.ambient_green_component);
        convert(token_b, new_material.ambient_blue_component);
      }
      else if(word[1] == 'd') // Diffuse color [Kd]; used for PBR as albedo
      {
        convert(token_r, new_material.diffuse_red_component);
        convert(token_g, new_material.diffuse_green_component);
        convert(token_b, new_material.diffuse_blue_component);
      }
      else if(word[1] == 's') // Specular color [Ks]
      {
        convert(token_r, new_material.specular_red_component);
        convert(token_g, new_material.specular_green_component);
        convert(token_b, new_material.specular_blue_component);
      }
      else if(word[1] == 'e') // Emissive color [Ke]
      {
        convert(token_r, new_material.emissive_red_component);
        convert(token_g, new_material.emissive_green_component);
        convert(token_b, new_material.emissive_blue_component);
      }
    }
    else if(word[0] == 'P') // PBR factor [P*]
    {
      std::string factor;
      file >> factor;

      if(word[1] == 'm') // Metallic factor [Pm]
      {
        convert(factor, new_material.metallic_factor);
      }
      else if(word[1] == 'r') // Roughness factor [Pr]
      {
        convert(factor, new_material.roughness_factor);
      }

      new_material.type = Types::MaterialType::MATERIAL_TYPE_PBR;
    }
    else if(word[0] == 'T')
    {
      std::string value;
      file >> value;

      if(word[1] == 'r') // Transparency (reversed, 1 - value) [Tr]
      {
        convert(value, new_material.transparency);
        new_material.transparency = 1.0 - new_material.transparency;
      }
      /*else if( word[1] == 'i' ) // Transmission filter [Ti]
      {

      }*/
    }
    else if(word[0] == 'b') // Bump or normal map [bump]
    {
      std::string texture_filename;
      file >> texture_filename;
      texture_filename.insert(0, parent_directory + "/");

      new_material.normal_map_filename = texture_filename;

      // RM: if normal map filename contains "normal" or "nrm", consider it a
      // normal map instead of a bump map
      const std::string lower_tex_filename =
          boost::algorithm::to_lower_copy(texture_filename);
      if(lower_tex_filename.find("normal") != std::string::npos ||
         lower_tex_filename.find("nrm") != std::string::npos)
        new_material.has_normal_map = true;
    }
    else if(word == "norm") // Normal map [norm]
    {
      std::string texture_filename;
      file >> texture_filename;
      texture_filename.insert(0, parent_directory + "/");

      new_material.normal_map_filename = texture_filename;
      new_material.has_normal_map = true;
    }
    else if(word[0] == 'd') // Transparency [d]
    {
      std::string value;
      file >> value;

      convert(value, new_material.transparency);
    }
    else if(word[0] == 'i') // Illumination [i]
    {
      // empty canvas, not yet supported

      std::string value;
      file >> value;

      switch(std::stoul(value))
      {
      case 0: // Color ON, ambient OFF
        break;

      case 1: // Color ON, ambient ON
        break;

      case 2: // Highlight ON
        break;

      case 3: // Reflection ON, raytrace ON
        break;

      case 4: // Transparency: glass ON; reflection: raytrace ON
        break;

      case 5: // Reflection: Fresnel ON, raytrace ON
        break;

      case 6: // Transparency: refraction ON; reflection: Fresnel OFF, raytrace
              // ON
        break;

      case 7: // Transparency: refraction ON; reflection: Fresnel ON, raytrace
              // ON
        break;

      case 8: // Reflection ON, raytrace OFF
        break;

      case 9: // Transparency: glass ON; reflection: raytrace OFF
        break;

      case 10: // Casts shadows onto invisible surfaces
        break;

      default:
        break;
      }
    }

    std::getline(file, word); // Recovering the rest of the line, to avoid the
                              // same line being taken twice
  }

  materials.push_back(new_material);
}

/**
 * Read OBJ file and store mesh data in an intermediate vector representation.
 */
template< typename CoordType,
          typename CoordNType,
          typename CoordTType,
          typename CoordCType,
          typename IndexType,
          typename MaterialType >
void
read_obj_file(const std::string &file_path,
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
  unsigned int cpt_l = count_file_lines(file_path);
  points_coords.clear();
  points_coords.reserve(
      (size_t)(cpt_l * 0.4f)); // TO DO find exactly the nb of points
  normals_coords.clear();
  normals_coords.reserve((size_t)(cpt_l * 0.4f));
  face_indices.clear();
  face_indices.reserve(
      (size_t)(cpt_l * 0.6f)); // TO DO find exactly the nb of faces
  texture_face_indices.clear();

  bool usemtl_found = false;
  IndexType current_material_id = -1;

  std::ifstream file(file_path);

  if(file.is_open())
  {
    std::string line;
    CoordType coord;
    CoordNType n_coord;
    CoordTType t_coord;
    CoordCType c_coord;
    int64_t f_ind, t_ind, n_ind;

    while(safe_getline(file, line))
    {
      std::vector< std::string > tokens = split(line, "\t ");

      size_t tokens_size = tokens.size();
      if(!tokens.empty())
      {
        if(is_equal(tokens[0], "o")) // MESH NAME
        {
          // TODO : find a place to keep the mesh name
        }
        if(is_equal(tokens[0], "mtllib"))
        {
          // append mtl file name to obj file path
          std::string str_parent_directory = get_parent_directory(file_path);
          if(str_parent_directory.empty())
            str_parent_directory = ".";
          std::string mtl_file_name = str_parent_directory + "/" + tokens[1];

          // remove ending '\r' if present
          // TODO-elo-note:  is it still necessary because
          //                now we use safe_getline ?
          if(mtl_file_name.back() == '\r')
            mtl_file_name.pop_back();

          // search texture file name in mtl file
          read_mtl_file(mtl_file_name, materials);
          /*//TODO-elo-dbg*/ dbg_display_materials(materials);
        }
        else if(is_equal(tokens[0], "vn")) // NORMAL POINT
        {
          std::vector< CoordNType > normal;

          for(unsigned int i = 1; i < 4; ++i)
          {
            convert(tokens[i], n_coord);
            normal.push_back(n_coord);
          }
          normals_coords.push_back(normal);
        }
        else if(is_equal(tokens[0], "vt")) // TEXTURE POINT
        {
          std::vector< CoordTType > tex_coord;

          for(unsigned int i = 1; i < 3; ++i)
          {
            convert(tokens[i], t_coord);
            tex_coord.push_back(t_coord);
          }
          texture_coords.push_back(tex_coord);
        }
        else if(is_equal(tokens[0], "v")) // POINT
        {
          std::vector< CoordType > point;
          std::vector< CoordCType > color;

          if(tokens_size == 5 || tokens_size == 8)
          { // homogeneous point
            throw std::runtime_error(
                "Reader::read_obj_file -> homogeneous point find in file : not "
                "implemented yet.");
          }
          else if(tokens_size == 4 || tokens_size == 7)
          { // cartesian point
            for(unsigned int i = 1; i < 4; ++i)
            {
              convert(tokens[i], coord);
              point.push_back(coord);
            }
          }
          else
          {
            throw std::runtime_error("Reader::read_obj_file -> wrong number of "
                                     "coordinates for a point.");
          }
          if(tokens_size == 7 || tokens_size == 8)
          { // colored point
            for(size_t i = tokens_size - 3; i < tokens_size; ++i)
            {
              convert(tokens[i], c_coord);
              color.push_back(c_coord);
            }
          }
          points_coords.push_back(point);
          vertex_color_coords.push_back(color);
        }
        else if(is_equal(tokens[0], "f")) // FACE
        {
          std::vector< IndexType > face_points, face_textures, face_normals;
          bool has_normal = false, has_texture = false,
               is_first_face_token = true;
          for(size_t i = 1; i < tokens_size; ++i)
          {
            std::vector< std::string > face_tokens = split(
                tokens[i],
                "/",
                true); // the last "true" parameter stands for keep empty tokens

            // Determined whether the face has normals and/or texture
            //  point/texture/normal
            if(is_first_face_token)
            {
              if(face_tokens.size() == 2)
              {
                has_texture = true;
              }
              else if(face_tokens.size() == 3)
              {
                has_normal = true;
                if(!face_tokens[1].empty())
                  has_texture = true;
              }
              is_first_face_token = false;
            }

            convert(face_tokens[0], f_ind);
            face_points.push_back(static_cast< IndexType >(
                f_ind < 0
                    ? f_ind + points_coords.size()
                    : f_ind - 1)); // the indices needs to start to 0, which is
                                   // not the case in obj format (starts to 1)

            if(has_texture) // TODO : maybe have to check if all the face tokens
                            // are formed the same way (eg : not 0//12 1/1 2/3/9)
                            // as the standard explicitly set it
            {
              convert(face_tokens[1], t_ind);
              face_textures.push_back(static_cast< IndexType >(
                  t_ind < 0 ? t_ind + texture_coords.size() : t_ind - 1));
            }
            if(has_normal) // TODO : maybe have to check if all the face tokens
                           // are formed the same way (eg : not 0//12 1/1 2/3/9)
                           // as the standard explicitly set it
            {
              convert(face_tokens[2], n_ind);
              face_normals.push_back(static_cast< IndexType >(
                  n_ind < 0 ? n_ind + normals_coords.size() : n_ind - 1));
            }
          }

          face_indices.push_back(face_points);
          texture_face_indices.push_back(face_textures);
          normal_face_indices.push_back(face_normals);
          if(usemtl_found)
            face_material.push_back(current_material_id);
        }
        else if(is_equal(tokens[0], "usemtl")) // MATERIAL
        {
          usemtl_found = true;
          std::string mtl_name = tokens[1];
          current_material_id = -1;
          // std::cout << "---------> mtl_name: " << mtl_name << std::endl;
          for(IndexType i = 0; i < materials.size(); i++)
          {
            // std::cout << "-------------------> materials[i].name: " <<
            // materials[i].name << std::endl;
            if(materials[i].name == mtl_name)
            {
              current_material_id = i;
              break;
            }
          }
        }
      }
    }
    file.close();
  }
  else
  {
    throw std::runtime_error(
        "Reader::read_obj_file -> input file failed to open.");
  }
}

} // namespace IO
} // namespace FEVV

