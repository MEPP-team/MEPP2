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


enum ply_property { VERTEX_PROPERTY, EDGE_PROPERTY, FACE_PROPERTY };

/**
 * Read PLY file and store mesh data in an intermediate vector representation.
 */
template< typename CoordType,
          typename CoordNType,
          typename CoordTType,
          typename CoordCType,
          typename IndexType >
void
read_ply_file(std::string file_path,
              std::vector< std::vector< CoordType > > &points_coords,
              std::vector< std::vector< CoordNType > > &normals_coords,
              std::vector< std::vector< CoordTType > > &texture_coords,
              std::vector< std::vector< CoordCType > > &vertex_color_coords,
              std::vector< std::vector< IndexType > > &face_indices,
              std::vector< std::vector< IndexType > > &texture_face_indices)
{
  points_coords.clear();
  normals_coords.clear();
  texture_coords.clear();
  vertex_color_coords.clear();
  face_indices.clear();
  texture_face_indices.clear();

  std::ifstream file(file_path);

  if(file.is_open())
  {
    std::string line;
    size_t tokens_size;
    std::vector< std::string > tokens;

    ply_property prop = VERTEX_PROPERTY;
    bool vertices_have_color = false, vertices_have_alpha = false,
         vertices_have_normal = false, faces_have_texture = false;

    unsigned int vertex_dim = 3u, normal_dim = 3u, color_dim = 3u,
                 alpha_dim = 3u, tex_dim = 2u, offset, face_degree;

    unsigned long nb_read_vertices = 0ul, nb_total_vertices = 0ul,
                  nb_read_edges = 0ul, nb_total_edges = 0ul,
                  nb_read_faces = 0ul, nb_total_faces = 0ul,
                  nb_face_texture_point = 0u;

    CoordType coord;
    CoordNType n_coord;
    CoordTType t_coord;
    CoordCType c_coord;
    IndexType f_ind;

    // HEADER READ
    while(safe_getline(file, line))
    {
      tokens = split(line, "\t ");
      tokens_size = tokens.size();

      if(tokens_size == 0) // Nothing on the line
        continue;
      else
      {
        if(tokens[0].compare("element") == 0)
        {
          if(tokens_size > 2)
          {
            if(tokens[1].compare("vertex") == 0)
            {
              convert(tokens[2], nb_total_vertices);
              prop = VERTEX_PROPERTY;
            }
            else if(tokens[1].compare("edge") == 0)
            {
              convert(tokens[2], nb_total_edges);
              prop = EDGE_PROPERTY;
            }
            else if(tokens[1].compare("face") == 0)
            {
              convert(tokens[2], nb_total_faces);
              prop = FACE_PROPERTY;
            }
          }
          else
            throw std::runtime_error("Reader::read_ply_file -> wrong tokens "
                                     "number for header line.");
        }
        else if(tokens[0].compare("property") == 0)
        {
          if(tokens_size > 2)
          {
            switch(prop) // TODO : we can do a better properties handling (and
                         // add material properties handling)
            {
            case VERTEX_PROPERTY:
              if(tokens[2].compare("nx") == 0) // eg : property float nx
                vertices_have_normal = true;
              if(tokens[2].compare("red") == 0) // eg : property uchar red
                vertices_have_color = true;
              if(tokens[2].compare("alpha") == 0) // eg : property uchar alpha
                vertices_have_alpha = true;
              break;
            case EDGE_PROPERTY:
              // No edge properties handle yet
              break;
            case FACE_PROPERTY:
              if(tokens_size > 4) // eg : property uint8 float texcoord
              {
                if(tokens[4].compare("texcoord") == 0)
                  faces_have_texture = true;
              }
              break;
            }
          }
          else
            throw std::runtime_error("Reader::read_ply_file -> wrong tokens "
                                     "number for header line.");
        }
        else if((tokens[0].compare("end_header") == 0))
          break;
      }
    }

    // VERTICES READ
    while(nb_read_vertices < nb_total_vertices)
    {
      if((safe_getline(file, line)))
      {
        tokens = split(line, "\t ");
        tokens_size = tokens.size();
        unsigned int i = 0u;

        if(tokens_size == 0) // Nothing on the line (the comment line case is
                             // not detected, as the comment tag is a full
                             // string, it would have decrease the performance)
          continue;
        else if(tokens_size == vertex_dim +
                                   (vertices_have_normal ? normal_dim : 0) +
                                   (vertices_have_color ? color_dim : 0) +
                                   (vertices_have_alpha ? alpha_dim : 0))
        {
          std::vector< CoordType > point;
          for(; i < vertex_dim; ++i)
          {
            convert(tokens[i], coord);
            point.push_back(coord);
          }
          points_coords.push_back(point);

          if(vertices_have_normal)
          {
            std::vector< CoordNType > normal;
            offset = vertex_dim + normal_dim;
            for(; i < offset; ++i)
            {
              convert(tokens[i], n_coord);
              normal.push_back(n_coord);
            }
            normals_coords.push_back(normal);
          }

          if(vertices_have_color)
          {
            std::vector< CoordCType > color;
            offset = vertex_dim + (vertices_have_normal ? normal_dim : 0) +
                     color_dim +
                     (vertices_have_alpha ? alpha_dim
                                          : 0); // corrected by Vincent V.
            for(; i < offset; ++i)
            {
              convert(tokens[i], c_coord);
              color.push_back(c_coord);
            }
            vertex_color_coords.push_back(color);
          }

          nb_read_vertices++;
        }
        else
        {
          // DBG std::cout << "tokensSize=" << tokensSize << " vertexDim=" <<
          // vertexDim << " verticesHaveNormal=" << verticesHaveNormal << "
          // normalDim=" << normalDim << " verticesHaveColor=" <<
          // verticesHaveColor << " colorDim=" << colorDim << "
          // verticesHaveAlpha=" << verticesHaveAlpha << " alphaDim=" << alphaDim
          // << std::endl;
          // TODO-elo  vertex texture coordinates not yet supported in this PLY
          // reader
          // TODO-elo  supported vertex line format:
          // TODO-elo     vt_coordinates_x_y_z  [vt_normal_x_y_z]
          // [vt_color_x_y_z [vt_color_alpha_a]]
          throw std::runtime_error(
              std::string("Reader::read_ply_file -> vertex line '") + line +
              std::string("'do not fit to the format."));
        }
      }
      else
        throw std::runtime_error(
            "Reader::read_ply_file -> wrong number of vertices specified or "
            "read error occured while parsing vertices.");
    }

    // FACES READ
    while(nb_read_faces < nb_total_faces)
    {
      if((safe_getline(file, line)))
      {
        tokens = split(line, "\t ");
        tokens_size = tokens.size();
        unsigned int i = 1u;

        if(tokens_size == 0) // Nothing on the line (the comment line case is
                             // not detected, as the comment tag is a full
                             // string, it would have decrease the performance)
          continue;
        else
        {
          std::vector< IndexType > face;
          convert(tokens[0], face_degree);
          offset = face_degree + 1;

          if(tokens_size >= offset)
          {
            for(; i < offset; ++i)
            {
              convert(tokens[i], f_ind);
              face.push_back(f_ind);
            }
            face_indices.push_back(face);

            if(faces_have_texture)
            {
              i++; // The number of tex coord is not taken into account for now
                   // (a valid number of tex coord is faceDegree*texDim)

              std::vector< IndexType > face_textures;
              for(unsigned int j = 0; j < face_degree; ++i)
              {
                std::vector< CoordTType > tex_coord;
                for(; i < tex_dim; ++i)
                {
                  convert(tokens[i], t_coord);
                  tex_coord.push_back(t_coord);
                }
                texture_coords.push_back(tex_coord);
                face_textures.push_back(nb_face_texture_point++);
              }
              texture_face_indices.push_back(face_textures);
            }
            nb_read_faces++;
          }
          else
            throw std::runtime_error(
                "Reader::read_ply_file -> face line do not fit to the format.");
        }
      }
      else
        throw std::runtime_error(
            "Reader::read_ply_file -> wrong number of faces specified or read "
            "error occured while parsing faces.");
    }

    // EDGES READ
    while(nb_read_edges < nb_total_edges)
    {
      if((safe_getline(file, line)))
      {
        tokens = split(line, "\t ");
        tokens_size = tokens.size();
        unsigned int i = 0u;

        if(tokens_size == 0) // Nothing on the line (the comment line case is
                             // not detected, as the comment tag is a full
                             // string, it would have decrease the performance)
          continue;
        else
        {
          std::vector< IndexType > edge;

          if(tokens_size >= 2)
          {
            for(; i < 2; ++i)
            {
              convert(tokens[i], f_ind);
              edge.push_back(f_ind);
            }
            face_indices.push_back(edge); // add edge as a two vertices face

            nb_read_edges++;
          }
          else
            throw std::runtime_error(
                "Reader::read_ply_file -> edge line do not fit to the format.");
        }
      }
      else
        throw std::runtime_error(
            "Reader::read_ply_file -> wrong number of edges specified or read "
            "error occured while parsing edges.");
    }

    file.close();
  }
  else
  {
    throw std::runtime_error(
        "Reader::read_ply_file -> input file failed to open.");
  }
}

} // namespace IO
} // namespace FEVV

