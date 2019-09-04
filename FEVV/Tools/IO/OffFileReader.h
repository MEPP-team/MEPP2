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
#include <sstream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <string>
#include <stdexcept> // for runtime_error

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "FEVV/Tools/IO/StringUtilities.hpp"

namespace FEVV {
namespace IO {

using namespace StrUtils;
using namespace FileUtils;


/**
 * Read _OFF file and store mesh data in an intermediate vector representation.
 */
template< typename CoordType,
          typename CoordNType,
          typename CoordTType,
          typename CoordCType,
          typename IndexType >
void
read_off_file(std::string file_path,
              std::vector< std::vector< CoordType > > &points_coords,
              std::vector< std::vector< CoordNType > > &normals_coords,
              std::vector< std::vector< CoordTType > > &texture_coords,
              std::vector< std::vector< CoordCType > > &vertex_color_coords,
              std::vector< std::vector< IndexType > > &face_indices,
              std::vector< std::vector< CoordCType > > &face_color_coords)
{
  points_coords.clear();
  normals_coords.clear();
  texture_coords.clear();
  face_indices.clear();
  vertex_color_coords.clear();
  face_color_coords.clear();

  std::ifstream file(file_path);

  if(! file.is_open())
  {
    throw std::runtime_error(
        "Reader::read_off_file -> failed to open input file.");
  }

  std::string line_str, word;
  std::istringstream line_ss;

  bool vertices_have_color = false;
  bool vertices_have_normals = false;
  bool vertices_have_texture = false;

  unsigned int vertex_dim = 3u,
               normal_dim = 3u,
               texture_dim = 2u,
               color_dim = 3u;
  unsigned int face_degree;
  unsigned long nb_read_vertices = 0ul,
                nb_total_vertices = 0ul;
  unsigned long nb_total_edges = 0ul;
  unsigned long nb_read_faces = 0ul,
                nb_total_faces = 0ul;

  CoordType coord;
  CoordNType n_coord;
  CoordTType t_coord;
  CoordCType c_coord;
  IndexType f_ind;

  // Find the OFF flag
  getline_skip_comment(file, line_str, line_ss);
  line_ss >> word;
  if((word.size() < 3) || (word.substr(word.size() - 3) != "OFF"))
  {
    throw std::runtime_error("Reader::read_off_file -> file format is "
                             "unknown (no [ST][C][N]OFF header). Other "
                             "OFF formats are not implemented yet.");
  }

  vertices_have_normals = (word.size() >= 4) &&
                          (word.substr(word.size() - 4) == "NOFF");
                          // end with NOFF
  vertices_have_color = (word == "COFF") || (word == "CNOFF") ||
                        (word == "STCOFF") || (word == "STCNOFF");
  vertices_have_texture = (word.substr(0, 2) == "ST"); // begin with ST

  // Find the element numbers
  getline_skip_comment(file, line_str, line_ss);
  line_ss >> nb_total_vertices >> nb_total_faces >> nb_total_edges;

  points_coords.reserve(nb_total_vertices);
  face_indices.reserve(nb_total_faces);
  if(vertices_have_normals)
    normals_coords.reserve(nb_total_vertices);
  if(vertices_have_color)
    vertex_color_coords.reserve(nb_total_vertices);
  if(vertices_have_texture)
    texture_coords.reserve(nb_total_vertices);

  // read vertices
  std::vector< CoordType > point;
  std::vector< CoordNType > normal;
  std::vector< CoordCType > color;
  std::vector< CoordTType > tex_coord;
  while(nb_read_vertices < nb_total_vertices)
  {
    if(getline_skip_comment(file, line_str, line_ss))
    {
      // get vertex coordinates
      point.clear();
      for(unsigned int i = 0; i < vertex_dim; ++i)
      {
        line_ss >> coord;
        point.push_back(coord);
      }
      points_coords.push_back(point);

      // get vertex normal
      if(vertices_have_normals)
      {
        normal.clear();
        for(unsigned int i = 0; i < normal_dim; ++i)
        {
          line_ss >> n_coord;
          normal.push_back(n_coord);
        }
        normals_coords.push_back(normal);
      }

      // get vertex color
      if(vertices_have_color)
      {
        color.clear();
        for(unsigned int i = 0; i < color_dim; ++i)
        {
          line_ss >> c_coord;
          color.push_back(c_coord);
        }
        vertex_color_coords.push_back(color);
      }

      // get vertex texture coordinates
      if(vertices_have_texture)
      {
        tex_coord.clear();
        for(unsigned int i = 0; i < texture_dim; ++i)
        {
          line_ss >> t_coord;
          tex_coord.push_back(t_coord);
        }
        texture_coords.push_back(tex_coord);
      }

      nb_read_vertices++;
      }
    else
    {
      throw std::runtime_error(
          "Reader::read_off_file -> failed to read vertex line.");
    }
  }

  // read faces
  std::vector< IndexType > face;
  while(nb_read_faces < nb_total_faces)
  {
    if(getline_skip_comment(file, line_str, line_ss))
    {
      // read face vertices
      face.clear();
      line_ss >> face_degree;
      for(unsigned int i = 0; i < face_degree; ++i)
      {
        line_ss >> f_ind;
        face.push_back(f_ind);
      }
      face_indices.push_back(face);

      // read face color if present
      if(line_ss.good())
      {
        // note: .good() is true if previous read was ok and
        //       end of line is not yet reached, so it means
        //       that extra data are available ;
        //       beware : some files have lines that end with
        //       a whitespace, so .good() is not enough to be
        //       sure that a color is provided for this face

        color.clear();
        for(unsigned int i = 0; i < color_dim; ++i)
        {
          line_ss >> c_coord;
          color.push_back(c_coord);
        }

        // for the case where the line ends with a whitespace,
        // or the color information is ill-formed
        if(! line_ss.fail())
          face_color_coords.push_back(color);
      }

      nb_read_faces++;
    }
    else
    {
      throw std::runtime_error(
          "Reader::read_off_file -> failed to read face line.");
    }
  }

  file.close();
}

} // namespace IO
} // namespace FEVV

