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
#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "FEVV/Tools/IO/FaceIndicesUtilities.h"

#include "FEVV/Tools/IO/ObjFileReader.h"
#include "FEVV/Tools/IO/OffFileReader.h"
#include "FEVV/Tools/IO/PlyFileReader.h"
#include "FEVV/Tools/IO/MshFileReader.h"

#include <vector>
#include <array>
#include <string>

#ifdef FEVV_USE_VTK
#include "FEVV/Tools/IO/VtkFileReader.h"
#endif

#include "FEVV/Types/Material.h"

namespace FEVV {
namespace DataStructures {
namespace AIF {

inline
AIFMeshReader::ptr_output
AIFMeshReader::read(const std::string &filePath)
{
  using namespace FileUtils;
  using namespace StrUtils;

  if(!boost::filesystem::exists(filePath))
    throw std::invalid_argument(
        "Reader::read -> input file path refer to no existing file.");
  else if(!has_extension(filePath))
    throw std::invalid_argument(
        "Reader::read -> input file path find without any extension.");

  // std::cout << "AIFReading of \"" << get_file_full_name(filePath) << "\"" <<
  // std::endl;

  std::vector< std::string > validExtensions = {
      ".obj", ".off", ".coff", ".ply", ".msh"};
  std::vector< std::string > validVtkExtensions = {".vtk", ".vtp", ".vtu"};
  if(!(has_extension(filePath, validExtensions)
#ifdef FEVV_USE_VTK
       || has_extension(filePath, validVtkExtensions)
#endif
           ))
  {
    throw std::invalid_argument(
        "Reader::read -> input file extension can't be read (yet).");
  }

  std::vector< std::vector< coord_type > > points_coords;
  std::vector< std::vector< coordN_type > > normals_coords;
  std::vector< std::vector< coordT_type > > texture_coords;
  // if points_coords.size() == texture_coords.size(), then
  // texture coordinates are given by vertex
  // else texture coordinates are given by corner
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
  // std::vector< Material > materials;
  std::string texture_file_name; // TODO-elo-rm  when this file really supports
                                 // multitexture
  std::vector< FEVV::Types::Material > materials;
  std::vector< index_type > face_material;

  if(has_extension(filePath, ".obj"))
  {
    IO::read_obj_file(filePath,
                      points_coords,
                      normals_coords,
                      texture_coords,
                      vertex_color_coords,
                      faces_indices,
                      texture_face_indices,
                      normal_face_indices,
                      materials,
                      face_material);
    if(!materials.empty())
      texture_file_name = materials[0].diffuse_texture_filename;
  }
  else if(has_extension(filePath, ".off") || has_extension(filePath, ".coff"))
  {
    IO::read_off_file(filePath,
                      points_coords,
                      normals_coords,
                      texture_coords,
                      vertex_color_coords,
                      faces_indices,
                      face_color_coords);
  }
  else if(has_extension(filePath, ".ply"))
  {
    IO::read_ply_file(filePath,
                      points_coords,
                      normals_coords,
                      texture_coords,
                      vertex_color_coords,
                      faces_indices,
                      texture_face_indices);
  }
  else if(has_extension(filePath, ".msh"))
  {
    IO::read_gmsh_file(filePath,
                       points_coords,
                       normals_coords,
                       vertex_color_coords,
                       lines_indices,
                       lines_colors,
                       faces_indices,
                       faces_colors,
                       field_attributes,
                       field_names);
  }
#ifdef FEVV_USE_VTK
  else if(has_extension(filePath, ".vtk") || has_extension(filePath, ".vtp") ||
          has_extension(filePath, ".vtu"))
  {
    IO::read_vtk_or_vtp_or_vtu_file(filePath,
                                    points_coords,
                                    normals_coords,
                                    vertex_color_coords,
                                    lines_indices,
                                    lines_colors,
                                    faces_indices,
                                    faces_colors,
                                    field_attributes,
                                    field_names);
  }
#endif
  /////////////////////////////// MESH CREATION
  //////////////////////////////////////

  typedef AIFTopologyHelpers helpers;
  typedef AIFPropertiesHelpers PropHelpers;
  ptr_output outputMesh = AIFMesh::New();

  ///////////////////////// CREATE NEEDED PROPERTY MAPS
  //////////////////////////////


  // vertex-color must be defined for all vertices or none
  if(vertex_color_coords.size() != points_coords.size())
    vertex_color_coords.clear();
  for(auto &color : vertex_color_coords)
  {
    if(color.size() < 3 || color.size() > 4)
    {
      std::cout << "AIFMeshReader::read: found vertex color with size neither "
                   "3 nor 4. Disabling vertex color for all vertices."
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
    if(color.size() < 3 || color.size() > 4)
    {
      std::cout << "read_mesh(): found face color with size neither 3 nor 4. "
                   "Disabling face color for all vertices."
                << std::endl;
      face_color_coords.clear();
      break;
    }
  }

  // store texture filename
  if(!texture_file_name.empty())
  {
    // create a property map to store the texture name
    auto tfn_pm =
        outputMesh->AddAssocPropertyMap< helpers::ptr_mesh, std::string >(
            "m:texturefilename");
    if(!outputMesh->isAssocPropertyMap("m:texturefilename"))
      throw std::runtime_error(
          "Failed to create texture-filename property map.");

    (*tfn_pm)[outputMesh.get()] = texture_file_name;

    std::cout << "AIF property map for texture-filename created." << std::endl;
  }

  if (!field_names.empty())
  {
    auto it = field_names.begin();
    auto ite = field_names.end();
    for (; it != ite; ++it)
    {
      if (it->find("POINT_DATA_") != std::string::npos )
      { // property map associated to vertices

        // create property map to store vertex current data field
        outputMesh->AddPropertyMap< AIFVertex::ptr, std::vector<double> >(
          "v:datafield:"+ it->substr(11));
        if (!outputMesh->isPropertyMap< AIFVertex::ptr >("v:datafield:" + it->substr(11)))
          throw std::runtime_error(
            "Failed to create a vertex data field property map.");
      }
      else if (it->find("ELEMENT_DATA_") != std::string::npos)
      { // property map associated to faces

        // create property map to store face current data field
        outputMesh->AddPropertyMap< AIFFace::ptr, std::vector<double> >(
          "f:datafield:" + it->substr(13));
        if (!outputMesh->isPropertyMap< AIFFace::ptr >("f:datafield:" + it->substr(13)))
          throw std::runtime_error(
            "Failed to create a face data field property map.");
      }
      else if (it->find("CELL_DATA_") != std::string::npos)
      { // property map associated to faces

        // create property map to store face current data field
        outputMesh->AddPropertyMap< AIFFace::ptr, std::vector<double> >(
          "f:datafield:" + it->substr(10));
        if (!outputMesh->isPropertyMap< AIFFace::ptr >("f:datafield:" + it->substr(10)))
          throw std::runtime_error(
            "Failed to create a face data field property map.");
      }
    }
  }

  // import data into AIF datastructure

  bool useVertexNormal = false;
  bool useVertexColor = false;
  bool useVertexTextureCoord = false;
  bool useCornerTextureCoord = false;
  bool useFaceColor = false;
  bool useFaceNormal = false; // per-face vertices normals
  bool useLineColor = false;

  if(!texture_coords.empty())
  {
    if(texture_coords.size() == points_coords.size())
    {
      // texture coordinates are given by vertex
      useVertexTextureCoord = true;

      // create property map to store vertex texture-coord
      outputMesh->AddPropertyMap< AIFVertex::ptr, AIFMesh::PointUV >(
          "v:texcoord");
      if(!outputMesh->isPropertyMap< AIFVertex::ptr >("v:texcoord"))
        throw std::runtime_error(
            "Failed to create vertex-texture-coordinate property map.");

      std::cout << "AIF property map for vertex-texture-coordinate created."
                << std::endl;
    }
    else if(texture_face_indices.size() == faces_indices.size())
    {
      // texture coordinates are given by corner
      useCornerTextureCoord = true;

      // create property map to store corner texture-coord
      outputMesh->AddAssocPropertyMap< helpers::halfedge_descriptor,
                                       AIFMesh::PointUV >("h:texcoord");
      if(!outputMesh->isAssocPropertyMap("h:texcoord"))
        throw std::runtime_error(
            "Failed to create halfedge-texture-coordinate property map.");

      std::cout << "AIF property map for halfedge-texture-coordinate created."
                << std::endl;
    }
  }

  if(!normals_coords.empty())
  {
    // per-vertex normal or per-face-vertex normal ?
    // if there are as many normals as vertices, then it's per-vertex normal

    if(normals_coords.size() == points_coords.size())
    {
      useVertexNormal = true;

      // create property map to store vertex normal
      outputMesh->AddPropertyMap< AIFVertex::ptr, AIFMesh::Normal >("v:normal");
      if(!outputMesh->isPropertyMap< AIFVertex::ptr >("v:normal"))
        throw std::runtime_error(
            "Failed to create vertex-normal property map.");

      std::cout << "AIF property map for vertex-normal created." << std::endl;
    }
    else if(!normal_face_indices[0].empty())
    {
      // per-face-vertex normal
      useFaceNormal = true;

      // create property map to store face-vertex normal
      outputMesh->AddPropertyMap< AIFFace::ptr, AIFMesh::Normal >("f:normal");
      if(!outputMesh->isPropertyMap< AIFFace::ptr >("f:normal"))
        throw std::runtime_error("Failed to create face-normal property map.");

      std::cout << "AIF property map for face-normal created." << std::endl;
    }
  }

  if(!vertex_color_coords.empty())
  {
    useVertexColor = true;

    // create property map to store vertex color
    outputMesh->AddPropertyMap< AIFVertex::ptr, AIFMesh::Vector >("v:color");
    if(!outputMesh->isPropertyMap< AIFVertex::ptr >("v:color"))
      throw std::runtime_error("Failed to create vertex-color property map.");

    std::cout << "AIF property map for vertex-color created." << std::endl;
  }

  if(!face_color_coords.empty())
  {
    useFaceColor = true;

    // create property map to store face color
    outputMesh->AddPropertyMap< AIFFace::ptr, AIFMesh::Vector >("f:color");
    if(!outputMesh->isPropertyMap< AIFFace::ptr >("f:color"))
      throw std::runtime_error("Failed to create face-color property map.");

    std::cout << "AIF property map for face-color created." << std::endl;
  }

  /////////////////////////////// VERTICES ADDITION
  //////////////////////////////////////
  // bool useLineColor=(lines.size() == lines_c.size()) ;
  std::vector< helpers::vertex_type::ptr > vertices;
  vertices.reserve(points_coords.size());

  typedef std::vector< std::vector< coord_type > >::const_iterator point_it;
  int vertexNbr = 0;
  for(point_it itP = points_coords.begin(); itP != points_coords.end(); ++itP)
  {
    helpers::vertex_descriptor currentVertex = helpers::add_vertex(outputMesh);
    vertexNbr++;

    PropHelpers::set_point(
        outputMesh, currentVertex, (*itP)[0], (*itP)[1], (*itP)[2]);

    vertices.push_back(
        currentVertex); // Warning : must keep the vertex in a vector to handle
                        // vertex indices for faces

    // Vertex normal/color/texture-coords addition when present
    if(useVertexNormal || useVertexColor || useVertexTextureCoord)
    {
      if(useVertexNormal)
      {
        // DBG std::cout << "store normal by vertex #" << vertexNbr <<
        // std::endl;
        std::vector< coordN_type > &coords = normals_coords[vertexNbr - 1];

        if(coords.size() < 3 || coords.size() > 4)
          throw std::runtime_error("Support of normal coordinates of size "
                                   "neither 3 nor 4 not supported!");
        // store normal in property map
        AIFMesh::Normal vn(coords[0], coords[1], coords[2]);
        outputMesh->SetProperty< AIFVertex::ptr, AIFMesh::Normal >(
            "v:normal", currentVertex->GetIndex(), vn);
      }
      if(useVertexColor)
      {
        // DBG std::cout << "store color for vertex #" << vertexNbr <<
        // std::endl;
        std::vector< coordC_type > &coords = vertex_color_coords[vertexNbr - 1];
        // at this point we are sure that coords.size() == 3 because of
        // previous sanity checks

        // store color in property map
        AIFMesh::Vector vc(coords[0], coords[1], coords[2]);
        outputMesh->SetProperty< AIFVertex::ptr, AIFMesh::Vector >(
            "v:color", currentVertex->GetIndex(), vc);
      }
      if(useVertexTextureCoord)
      {
        // DBG std::cout << "store tex coord for vertex #" << vertexNbr <<
        // std::endl;
        std::vector< coordT_type > &coords = texture_coords[vertexNbr - 1];

        if(coords.size() != 2)
          throw std::runtime_error("Support of texture coordinates of size != "
                                   "2 not yet implemented!");

        // store texture coordinates in property map
        AIFMesh::PointUV texCoords = {coords[0], coords[1]};
        outputMesh->SetProperty< AIFVertex::ptr, AIFMesh::PointUV >(
            "v:texcoord", currentVertex->GetIndex(), texCoords);
      }
    }
    if (!field_names.empty())
    {
      auto it = field_names.begin();
      auto ite = field_names.end();
      auto it_d = field_attributes.begin();
      for (; it != ite; ++it, ++it_d)
      {
        if (it->find("POINT_DATA_") != std::string::npos)
        { // property map associated to vertices

          outputMesh->SetProperty< AIFVertex::ptr, std::vector<double> >(
            ("v:datafield:" + it->substr(11)), currentVertex->GetIndex(), (*it_d)[vertexNbr - 1]);
        }
      }
    }
  }

  /////////////////////////////// LINES ADDITION
  //////////////////////////////////////
  typedef std::vector< std::vector< index_type > >::const_iterator line_it;
  for(line_it itL = lines_indices.begin(); itL != lines_indices.end(); ++itL)
  {
    if((*itL).size() >
       3) // case where a line is used to describe a volume element (we do not
          // allows polyline loading, only segment loading+a particular use of
          // lines to descrime volume cells)
      continue;

    // find the corresponding vertices in the vertices std::vector (care about
    // the potential decrementation)
    helpers::vertex_descriptor firstVertex = vertices[(*itL)[0]];
    helpers::vertex_descriptor secondVertex = vertices[(*itL)[1]];

    // Dangling edge
    // create the edge in the mesh
    helpers::edge_type::ptr danglingEdge;

    if((danglingEdge = helpers::common_edge(firstVertex, secondVertex)) ==
       helpers::null_edge())
    { // we avoid multiple similar edges
      danglingEdge = helpers::add_edge(outputMesh);

      // link the edge and the vertices
      helpers::link_vertex_and_edge(
          firstVertex, danglingEdge, helpers::vertex_pos::FIRST);
      helpers::link_vertex_and_edge(
          secondVertex, danglingEdge, helpers::vertex_pos::SECOND);
    }
  }
  /////////////////////////////// FACES ADDITION
  //////////////////////////////////////
  typedef std::vector< std::vector< index_type > >::/*const_*/ iterator face_it;
  typedef std::vector< index_type >::const_iterator index_it;
  size_t faceId = 0;
  for(face_it itF = faces_indices.begin(); itF != faces_indices.end();
      ++itF, ++faceId)
  {
    if(itF->size() < 2)
    {
      throw std::runtime_error("Reader::read -> find a face with less than two "
                               "vertices : not complying case.");
    }
    else if(itF->size() == 2)
    {
      // find the corresponding vertices in the vertices std::vector (care about
      // the potential decrementation)
      helpers::vertex_descriptor firstVertex = vertices[(*itF)[0]];
      helpers::vertex_descriptor secondVertex = vertices[(*itF)[1]];

      // Dangling edge
      // create the edge in the mesh
      helpers::edge_type::ptr danglingEdge;
      if((danglingEdge = helpers::common_edge(firstVertex, secondVertex)) ==
         helpers::null_edge())
      { // we avoid multiple similar edges
        danglingEdge = helpers::add_edge(outputMesh);

        // link the edge and the vertices
        helpers::link_vertex_and_edge(
            firstVertex, danglingEdge, helpers::vertex_pos::FIRST);
        helpers::link_vertex_and_edge(
            secondVertex, danglingEdge, helpers::vertex_pos::SECOND);
      }
    }
    else
    {
      // create the face in the mesh
      helpers::face_descriptor currentFace = helpers::add_face(outputMesh);
      // helpers::face_descriptor currentFace = helpers::add_face(*outputMesh);
      /////////////////////////////////////////////////////////////////////////
      if(!FEVV::FaceIndicesUtils::are_all_face_indices_unique< index_type >(
             *itF))
      {
        std::vector< index_type > out_face_indices;
        std::vector< std::vector< index_type > > out_lines_indices;
        FEVV::FaceIndicesUtils::face_indices_to_face_and_lines_indices<
            index_type >(*itF, out_face_indices, out_lines_indices);
        *itF = out_face_indices;
        FEVV::FaceIndicesUtils::lines_indices_to_segments_indices(
            out_lines_indices);
        /////////////////////////////// LINES ADDITION
        //////////////////////////////////////
        typedef std::vector< std::vector< index_type > >::const_iterator
            line_it;
        for(line_it itL = out_lines_indices.begin();
            itL != out_lines_indices.end();
            ++itL)
        {
          // find the corresponding vertices in the vertices std::vector (care
          // about the potential decrementation)
          helpers::vertex_descriptor firstVertex = vertices[(*itL)[0]];
          helpers::vertex_descriptor secondVertex = vertices[(*itL)[1]];

          // Dangling edge
          // create the edge in the mesh
          helpers::edge_type::ptr danglingEdge;
          if((danglingEdge = helpers::common_edge(firstVertex, secondVertex)) ==
             helpers::null_edge())
          { // we avoid multiple similar edges
            danglingEdge = helpers::add_edge(outputMesh);

            // link the edge and the vertices
            helpers::link_vertex_and_edge(
                firstVertex, danglingEdge, helpers::vertex_pos::FIRST);
            helpers::link_vertex_and_edge(
                secondVertex, danglingEdge, helpers::vertex_pos::SECOND);
          }
        }
      }
      /////////////////////////////////////////////////////////////////////////
      // get the first vertex of the face to successively constructs the face's
      // edges (we can use back() here as vertices is explicitly an std::vector)
      helpers::vertex_descriptor prevVertex = vertices[itF->back()];
      helpers::vertex_descriptor currentVertex;
      helpers::edge_descriptor currentEdge;

      // loop over face's vertices
      index_type faceVertexId = 0;
      for(index_it it = itF->begin(); it != itF->end(); ++it, ++faceVertexId)
      {
        currentVertex = vertices[(*it)];

        // check if the 2 vertices ain't already link by an edge
        helpers::edge_descriptor currentEdge;
        if((currentEdge = helpers::common_edge(prevVertex, currentVertex)) ==
           helpers::null_edge())
        {
          // create the edge in the mesh
          currentEdge = helpers::add_edge(outputMesh);
          // currentEdge = helpers::add_edge(*outputMesh);

          // link the edge and the vertices
          helpers::link_vertex_and_edge(
              prevVertex, currentEdge, helpers::vertex_pos::FIRST);
          helpers::link_vertex_and_edge(
              currentVertex, currentEdge, helpers::vertex_pos::SECOND);
        }
        bool face_has_that_incident_edge =
            helpers::are_incident(currentFace, currentEdge);
        bool edge_has_that_incident_face =
            helpers::are_incident(currentEdge, currentFace);
        if(!face_has_that_incident_edge ||
           !edge_has_that_incident_face) // Even if data in file is not normal
                                         // (e.g. double dangling edge), file
                                         // loading can be performed.
        {
          // link the edge and the face
          helpers::link_edge_and_face(currentEdge, currentFace);
        }

        // corner (aka halfedge) texture
        if(useCornerTextureCoord)
        {
#if 0
          switch (texture_face_indices[faceId].size())
          {
          case 0: throw std::runtime_error("corner texture: found a face (index " + FEVV::StrUtils::convert(faceId) + ") without any texture coord index");
          case 1: throw std::runtime_error("corner texture: found a face (index " + FEVV::StrUtils::convert(faceId) + ") associated with only one texture coord index");
          case 2: throw std::runtime_error("corner texture: found a face (index " + FEVV::StrUtils::convert(faceId) + ") associated with only two texture coord indices");
          }
#endif
          if(texture_face_indices[faceId].size() ==
             faces_indices[faceId].size())
          { // in some files, it happens that a part of the faces have no
            // texture corners, while others have
            // retrieve texture coords
            index_type textureId = texture_face_indices[faceId][faceVertexId];
            std::vector< coordT_type > &coords = texture_coords[textureId];
            if(coords.size() != 2)
              throw std::runtime_error("Support of texture coordinates of size "
                                       "!= 2 not yet implemented!");
            AIFMesh::PointUV texCoords = {coords[0], coords[1]};

            // retrieve halfedge
            helpers::halfedge_descriptor he =
                helpers::AIFHalfEdge(prevVertex, currentVertex, currentFace);

            // fill in property map
            auto pm =
                outputMesh->GetAssocPropertyMap< helpers::halfedge_descriptor,
                                                 AIFMesh::PointUV >(
                    "h:texcoord");
            (*pm)[he] = texCoords;

            // DBG std::cout << __FILE__ << ":" << __LINE__ << " he=" << he << "
            // (*pm)[he][0]=" << (*pm)[he][0] << "   (*pm)[he][1]=" <<
            // (*pm)[he][1] << "   coords[0]=" << coords[0] << "   coords[1]=" <<
            // coords[1] << "   texCoords[0]=" << texCoords[0] << "
            // texCoords[1]=" << texCoords[1] << std::endl;
          }
        }

        // face-vertex-normal
        if(useFaceNormal)
        {
          // DBG std::cout << "compute normal for face #" << faceId << " (vertex
          // #" << faceVertexId << ")" << std::endl;

          // retrieve current vertex normal
          index_type normalId = normal_face_indices[faceId][faceVertexId];
          std::vector< coordN_type > &normal = normals_coords[normalId];
          if(normal.size() != 3)
            throw std::runtime_error(
                "Support of normals of size != 3 not yet implemented!");
          AIFMesh::Normal vnormal = {normal[0], normal[1], normal[2]};

          // retrieve face normal
          auto pm = outputMesh->GetPropertyMap< AIFFace::ptr, AIFMesh::Normal >(
              "f:normal");
          AIFMesh::Normal fnormal = (*pm)[currentFace];

          // computation of face normal (mean of the vertices normals)
          fnormal = fnormal + vnormal;
          if(it == itF->end() - 1)
            fnormal =
                fnormal / (faceVertexId +
                           1); // \todo: why to compute face normal not based on
                               // its underlying plane and just use the vertex
                               // normals to ensure correct orientation?

          // update face-normal map
          (*pm)[currentFace] = fnormal;
        }

        prevVertex = currentVertex;
      }

      if(currentFace->GetDegree() == 2)
      { // case of n>=2 edges describing the same non-oriented edge
        helpers::remove_face(
            currentFace,
            outputMesh); // we just remove that face (we always must ensure that
                         // a face is valid <=> degree>=3)
        continue;
      }

      if(useFaceColor)
      {
        // DBG std::cout << "store color for face #" << faceId << std::endl;
        std::vector< coordC_type > &color = face_color_coords[faceId];
        // at this point we are sure that color.size() == 3 because of
        // previous sanity checks

        // store color in property map
        AIFMesh::Vector fc(color[0], color[1], color[2]);
        outputMesh->SetProperty< AIFFace::ptr, AIFMesh::Vector >(
            "f:color", currentFace->GetIndex(), fc);
      }
      if (!field_names.empty())
      {
        auto it = field_names.begin();
        auto ite = field_names.end();
        auto it_d = field_attributes.begin();
        for (; it != ite; ++it, ++it_d)
        {
          if (it->find("ELEMENT_DATA_") != std::string::npos)
          { // property map associated to faces

            outputMesh->SetProperty< AIFFace::ptr, std::vector<double> >(
              ("f:datafield:" + it->substr(13)), currentFace->GetIndex(), (*it_d)[faceId]);
          }
          else if (it->find("CELL_DATA_") != std::string::npos)
          { // property map associated to faces

            outputMesh->SetProperty< AIFFace::ptr, std::vector<double> >(
              ("f:datafield:" + it->substr(10)), currentFace->GetIndex(), (*it_d)[faceId]);
          }
        }
      }
    }
  }

  // std::cout << "AIFReading of \"" << get_file_full_name(filePath) << "\"
  // Done." << std::endl;

  return outputMesh;
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV
