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
#include <cassert>
#include <cstdio>
#include <string>
#include <exception>

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

  if(file.is_open())
  {
    std::string line;
    char comment_tag = '#';
    size_t tokens_size;
    std::vector< std::string > tokens;

    bool vertices_have_color = false, vertices_have_normals = false,
         vertices_have_texture = false;

    unsigned int vertex_dim = 3u, normal_dim = 3u, texture_dim = 2u;
    unsigned int face_degree;
    unsigned long nb_read_vertices = 0ul, nb_total_vertices = 0ul;
    unsigned long nb_total_edges = 0ul;
    unsigned long nb_read_faces = 0ul, nb_total_faces = 0ul;

    CoordType coord;
    CoordNType n_coord;
    CoordTType t_coord;
    CoordCType c_coord;
    IndexType f_ind;

    while(safe_getline(
        file,
        line)) // Find the OFF flag even if it is not on the very first line
    {
      tokens = split(line, "\t ");
      tokens_size = tokens.size();

      if(tokens_size == 0) // Nothing on the line
        continue;
      else if(tokens[0].front() == comment_tag) // Comment line
        continue;
      else if(tokens_size > 1) // Not just one OFF flag
        throw std::runtime_error(
            "Reader::read_off_file -> file format is unknown (no [ST][C][N]OFF "
            "header). Other OFF formats are not implemented yet.");
      else
      {
        vertices_have_normals = (tokens[0].compare("NOFF") == 0) ||
                                (tokens[0].compare("CNOFF") == 0) ||
                                (tokens[0].compare("STNOFF") == 0) ||
                                (tokens[0].compare("STCNOFF") == 0);
        vertices_have_color = (tokens[0].compare("COFF") == 0) ||
                              (tokens[0].compare("CNOFF") == 0) ||
                              (tokens[0].compare("STCOFF") == 0) ||
                              (tokens[0].compare("STCNOFF") == 0);
        vertices_have_texture = (tokens[0].compare("STOFF") == 0) ||
                                (tokens[0].compare("STNOFF") == 0) ||
                                (tokens[0].compare("STCOFF") == 0) ||
                                (tokens[0].compare("STCNOFF") == 0);

        if(!vertices_have_normals && !vertices_have_color &&
           !vertices_have_texture && tokens[0].compare("OFF") != 0)
          throw std::runtime_error("Reader::read_off_file -> file format is "
                                   "unknown (no [ST][C][N]OFF header). Other "
                                   "OFF formats are not implemented yet.");
        break;
      }
    }

    while(safe_getline(file, line)) // Find the element numbers
    {
      tokens = split(line, "\t ");
      tokens_size = tokens.size();

      if(tokens_size == 0) // Nothing on the line
        continue;
      else if(tokens[0].front() == comment_tag) // Comment line
        continue;
      else if(tokens_size == 1) // Vertex dimension detetected
        convert(tokens[0], vertex_dim);
      else if(tokens_size == 3)
      { // Right number of tokens found
        convert(tokens[0], nb_total_vertices);
        convert(tokens[1], nb_total_faces);
        convert(tokens[2], nb_total_edges);
        break;
      }
      else // Not the right number of tokens found
        throw std::runtime_error("Reader::read_off_file -> wrong tokens number "
                                 "for vertices/faces/edges count.");
    }

    points_coords.reserve(nb_total_vertices);
    face_indices.reserve(nb_total_faces);
    if(vertices_have_normals)
      normals_coords.reserve(nb_total_vertices);
    if(vertices_have_color)
      vertex_color_coords.reserve(nb_total_vertices);
    if(vertices_have_texture)
      texture_coords.reserve(nb_total_vertices);

    while(nb_read_vertices < nb_total_vertices)
    {
      if((safe_getline(file, line)))
      {
        tokens = split(line, "\t ");
        tokens_size = tokens.size();
        unsigned int i = 0u;

        if(tokens_size == 0) // Nothing on the line
          continue;
        else if(tokens[0].front() == comment_tag) // Comment line
          continue;
        else if(tokens_size >=
                vertex_dim + (vertices_have_normals ? normal_dim : 0) +
                    (vertices_have_texture
                         ? texture_dim
                         : 0)) // Color is not taking in account as it can
                               // appear in different formats
        {
          std::vector< CoordType > point;
          for(; i < vertex_dim; ++i)
          {
            convert(tokens[i], coord);
            point.push_back(coord);
          }
          points_coords.push_back(point);

          if(tokens_size > vertex_dim) // If there is more than only point coord
          {
            if(vertices_have_normals)
            {
              std::vector< CoordNType > normal;
              for(; i < vertex_dim + normal_dim; ++i)
              {
                convert(tokens[i], n_coord);
                normal.push_back(n_coord);
              }
              normals_coords.push_back(normal);
            }

            if(vertices_have_color)
            {
              std::vector< CoordCType > color;
              for(;
                  i < tokens_size - ((vertices_have_texture) ? texture_dim : 0);
                  ++i)
              {
                convert(tokens[i], c_coord);
                color.push_back(c_coord);
              }

              // if (color[0] == color.size() - 1) // the first element is
              // optionally the size of the color vector: we do not take it into
              // account for the time being 	color.erase(color.begin());

              vertex_color_coords.push_back(color);
            }

            if(vertices_have_texture)
            {
              std::vector< CoordTType > tex_coord;
              for(; i < tokens_size; ++i)
              {
                convert(tokens[i], t_coord);
                tex_coord.push_back(t_coord);
              }
              texture_coords.push_back(tex_coord);
            }
          }
          nb_read_vertices++;
        }
        else
          throw std::runtime_error(
              "Reader::read_off_file -> vertex line do not fit to the format.");
      }
      else
        throw std::runtime_error(
            "Reader::read_off_file -> wrong number of vertices specified or "
            "read error occured while parsing vertices.");
    }

    while(nb_read_faces < nb_total_faces)
    {
      if((safe_getline(file, line)))
      {
        tokens = split(line, "\t ");
        tokens_size = tokens.size();
        unsigned int i = 1u;

        if(tokens_size == 0) // Nothing on the line
          continue;
        else if(tokens[0].front() == comment_tag) // Comment line
          continue;
        else
        {
          std::vector< IndexType > face;
          convert(tokens[0], face_degree);

          if(tokens_size >= face_degree + 1)
          {
            for(; i < face_degree + 1; ++i)
            {
              convert(tokens[i], f_ind);
              face.push_back(f_ind);
            }
            face_indices.push_back(face);

            if(tokens_size > face_degree + 1)
            {
              std::vector< CoordCType > color;
              for(; i < tokens_size; ++i)
              {
                convert(tokens[i], c_coord);
                color.push_back(c_coord);
              }

              face_color_coords.push_back(color);
            }
            nb_read_faces++;
          }
          else
            throw std::runtime_error(
                "Reader::read_off_file -> face line do not fit to the format.");
        }
      }
      else
        throw std::runtime_error(
            "Reader::read_off_file -> wrong number of faces specified or read "
            "error occured while parsing faces.");
    }
    file.close();
  }
  else
  {
    throw std::runtime_error(
        "Reader::read_off_file -> input file failed to open.");
  }
}

} // namespace IO
} // namespace FEVV

