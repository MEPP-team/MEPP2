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

#include "FEVV/Tools/IO/ObjFileWriter.h"
#include "FEVV/Tools/IO/OffFileWriter.h"
#ifdef FEVV_USE_VTK
#include "FEVV/Tools/IO/VtkFileWriter.h"
#endif
#include "FEVV/Tools/IO/MshFileWriter.h"

#include "FEVV/Types/Material.h"

namespace FEVV {
namespace DataStructures {
namespace AIF {

inline
void
AIFMeshWriter::write(const ptr_input inputMesh, const std::string &filePath)
{
  write(*inputMesh, filePath);
}

inline
void
AIFMeshWriter::write(/*const*/ input_type &inputMesh,
                     const std::string &filePath)
{
  using namespace FileUtils;

  if(!has_extension(filePath))
    throw std::invalid_argument(
        "Writer::write -> output file path find without any extension.");

  // std::cout << "AIFWriting of \"" << get_file_full_name(filePath) << "\"" <<
  // std::endl;

  std::vector< std::string > validExtensions = {
      ".obj", ".off", ".coff", "ply", ".msh"};
  std::vector< std::string > validVtkExtensions = {".vtk", ".vtp", ".vtu"};
  if(!(has_extension(filePath, validExtensions)
#ifdef FEVV_USE_VTK
       || has_extension(filePath, validVtkExtensions)
#endif
           ))
  {
    throw std::invalid_argument(
        "Writer::write -> output file extension can't be read (yet).");
  }

  ///////////////////////////////////////////////////////////////////////////
  typedef AIFTopologyHelpers helpers;
  typedef AIFPropertiesHelpers PropHelpers;
  typedef helpers::vertex_container::const_iterator vertex_iterator;
  typedef boost::iterator_range< vertex_iterator > vertex_range;
  typedef helpers::edge_container::const_iterator edge_iterator;
  typedef boost::iterator_range< edge_iterator > edge_range;
  typedef helpers::face_container::const_iterator face_iterator;
  typedef boost::iterator_range< face_iterator > face_range;
  typedef helpers::vertex_container_in_face::const_iterator
      vertex_face_iterator;
  typedef boost::iterator_range< vertex_face_iterator > vertex_face_range;

  std::vector< std::vector< coord_type > > points_coords;
  std::vector< std::vector< coordN_type > > normals_coords;
  std::vector< std::vector< coordT_type > > texture_coords;
  std::vector< std::vector< coordC_type > > vertex_color_coords;
  std::vector< std::vector< coordC_type > > face_color_coords;
  std::vector< std::vector< index_type > > lines_indices, faces_indices;
  std::vector< std::vector< index_type > > texture_face_indices;
  std::vector< std::vector< index_type > > normal_face_indices;
  std::vector< std::vector< coordC_type > > line_colors_coords;
  std::vector< index_type > face_material;
  std::vector< FEVV::Types::Material > materials;

  std::vector< std::vector< std::vector< double > > > field_attributes;
  std::vector< std::string > field_names;
  // which vertex attribute are available for export ?

  bool useNorm = false;
  bool useVertexColor = false;
  bool useFaceColor = false;
  bool useVertexTextureCoord = false;
  bool useVertexDatafield = false;
  bool useFaceDatafield = false;

#ifdef FEVV_USE_VTK
  bool useEdgeColor =
      false; // not a common usage, therfore we restric it to vtk
  if(inputMesh.isPropertyMap< AIFEdge::ptr >("e:color"))
    useEdgeColor = true;
#endif

  if(inputMesh.isPropertyMap< AIFVertex::ptr >("v:normal"))
    useNorm = true;
  if(inputMesh.isPropertyMap< AIFVertex::ptr >("v:color"))
    useVertexColor = true;
  if(inputMesh.isPropertyMap< AIFFace::ptr >("f:color"))
    useFaceColor = true;
  if(inputMesh.isPropertyMap< AIFVertex::ptr >("v:texcoord"))
    useVertexTextureCoord = true;
  if (inputMesh.isAPropertyMapStartingWithPrefix< AIFVertex::ptr >("v:datafield:"))
    useVertexDatafield = true;
  if (inputMesh.isAPropertyMapStartingWithPrefix< AIFFace::ptr >("f:datafield:"))
    useFaceDatafield = true;
  
  std::vector<std::string> dfv_names = inputMesh.GetPropertyMapNamesStartingWithPrefix< AIFVertex::ptr >("v:datafield:");
  std::vector<std::string> dff_names = inputMesh.GetPropertyMapNamesStartingWithPrefix< AIFFace::ptr >("f:datafield:");
  
  unsigned int pointDim = 3;
  long vertexIndex = 0; // to be sure to start to 0
  std::map< helpers::vertex_descriptor, long >
      indexMap; // for correct corresponding between vertex order in file and
                // vertices around face indices

  // VERTICES IN CONTAINER
  vertex_range vRange = helpers::vertices(inputMesh);

  for(vertex_iterator itV = vRange.begin(); itV != vRange.end(); ++itV)
  {
    PropHelpers::Point p = PropHelpers::get_point(inputMesh, *itV);
    std::vector< coord_type > point;

    for(unsigned int i = 0; i < p.size(); ++i)
      point.push_back(p[i]);

    points_coords.push_back(point);

    indexMap[*itV] = vertexIndex;
    ++vertexIndex;

    // deal with vertex attributes here

    if(useNorm)
    {
      // store vertex normal in standard container
      AIFMesh::Normal &vn =
          inputMesh.GetProperty< AIFVertex::ptr, AIFMesh::Normal >(
              "v:normal", (*itV)->GetIndex());
      std::vector< coordN_type > normal(vn.cbegin(), vn.cend());
      normals_coords.push_back(normal);
    }

    if(useVertexColor)
    {
      // store vertex color in standard container
      AIFMesh::Vector &vc =
          inputMesh.GetProperty< AIFVertex::ptr, AIFMesh::Vector >(
              "v:color", (*itV)->GetIndex());
      std::vector< coordC_type > color(vc.cbegin(), vc.cend());
      vertex_color_coords.push_back(color);
    }

    if(useVertexTextureCoord)
    {
      // store vertex texture-coordinates in standard container
      AIFMesh::PointUV &uv =
          inputMesh.GetProperty< AIFVertex::ptr, AIFMesh::PointUV >(
              "v:texcoord", (*itV)->GetIndex());
      std::vector< coordT_type > vt(uv.cbegin(), uv.cend());
      texture_coords.push_back(vt);
    }

    if (useVertexDatafield)
    {
      auto it_s = dfv_names.begin(), it_se = dfv_names.end();
      size_t i = 0;
      for (; it_s != it_se; ++it_s, ++i)
      {
        const std::vector< double >&vdata =
          inputMesh.GetProperty< AIFVertex::ptr, std::vector< double > >(
            *it_s, (*itV)->GetIndex());

        if (itV == vRange.begin())
        {
          field_names.push_back("POINT_DATA_" + it_s->substr(12));
          field_attributes.resize(field_attributes.size() + 1);
        }

        field_attributes[i].push_back(vdata);
      }
    }
  }

  // EDGES IN CONTAINER (needed to take into accound dangling or isolated edges)
  edge_range eRange = helpers::edges(inputMesh);

  for(edge_iterator itE = eRange.begin(); itE != eRange.end(); ++itE)
  {
    if(helpers::is_dangling_edge(*itE) || helpers::is_isolated_edge(*itE))
    {
      std::vector< index_type > line;

      line.push_back(indexMap[(*itE)->get_first_vertex()]);
      line.push_back(indexMap[(*itE)->get_second_vertex()]);

      lines_indices.push_back(line);
#ifdef FEVV_USE_VTK
      if(useEdgeColor)
      {
        // store vertex color in standard container
        AIFMesh::Vector &vc =
            inputMesh.GetProperty< AIFEdge::ptr, AIFMesh::Vector >(
                "e:color", (*itE)->GetIndex());
        std::vector< coordC_type > color(vc.cbegin(), vc.cend());
        line_colors_coords.push_back(color);
      }
#endif
    }
  }

  // FACES IN CONTAINER
  std::set< AIFVertex::ptr > already_processed; // use to export dangling polylines only once
  face_range fRange = helpers::faces(inputMesh);

  for(face_iterator itF = fRange.begin(); itF != fRange.end(); ++itF)
  {
    std::vector< index_type > face;
    std::vector< index_type > normal_indices;
    std::vector< index_type > texture_indices;

    vertex_face_range vfRange = helpers::incident_vertices(*itF);

    for(vertex_face_iterator itVF = vfRange.begin(); itVF != vfRange.end();
        ++itVF)
    {
	  if ( has_extension(filePath, ".off") )
      {
        auto nb_incident_dangling_edges = helpers::num_incident_dangling_edge(*itVF);
        if(nb_incident_dangling_edges>1)
          throw std::runtime_error("Writer::write -> multiple dangling edge integration into face indices is not available for .off");
        else if ( (nb_incident_dangling_edges == 1) && (already_processed.find(*itVF)==already_processed.end()))
        {
          already_processed.insert(*itVF);
          //throw std::runtime_error("Writer::write -> dangling edge integration into face indices for .off has not been implemented yet");
          // 1) insert the dangling polyline
          AIFVertex::ptr current_v = *itVF;
          AIFEdge::ptr current_e = helpers::get_incident_dangling_edge(current_v);
          face.push_back(indexMap[*itVF]);
          int nb_indices_before = static_cast<int>(face.size());
          do
          {
            current_v = helpers::opposite_vertex(current_e, current_v);
            face.push_back(indexMap[current_v]);
            ///////////////////////////////////////////////////////////////////
            nb_incident_dangling_edges = helpers::num_incident_dangling_edge(current_v);
            auto evRange = helpers::incident_edges(current_v);
            auto iterE = evRange.begin();
            for(; iterE != evRange.end(); ++iterE)
            {
              if (*iterE != current_e)
              {
                current_e = *iterE;
                break;
              }
            }
          } while (nb_incident_dangling_edges==2); // we only deal with dangling polylines without bisection
          int nb_indices_after = static_cast<int>(face.size());
          std::vector<index_type> additional_indices;
          for (int i = nb_indices_after-2; i >= nb_indices_before; --i)
          {
            additional_indices.push_back(face[i]);
          }
          face.insert(face.end(), additional_indices.begin(), additional_indices.end());
          // 2) remove it from lines_indices
          // not done since a) it is costly and b) lines_indices are not used by off writer 
        }
      }
      face.push_back(
          indexMap[*itVF]); // for correct corresponding between vertex order in
                            // file and vertices around face indices

      // face's vertices attributes (normal & texture)
      if(useNorm)
        normal_indices.push_back(indexMap[*itVF]);
      if(useVertexTextureCoord)
        texture_indices.push_back(indexMap[*itVF]);
    }

    faces_indices.push_back(face);
    if(useNorm)
      normal_face_indices.push_back(normal_indices);
    if(useVertexTextureCoord)
      texture_face_indices.push_back(texture_indices);

    // face color
    if(useFaceColor)
    {
      // retrieve face-color from property map
      AIFMesh::Vector &fc =
          inputMesh.GetProperty< AIFFace::ptr, AIFMesh::Vector >(
              "f:color", (*itF)->GetIndex());

      // store face-color in standard container
      std::vector< coordC_type > color(fc.cbegin(), fc.cend());
      face_color_coords.push_back(color);
    }

    if (useFaceDatafield)
    {
      auto it_s = dff_names.begin(), it_se = dff_names.end();
      size_t i = 0;
      for (; it_s != it_se; ++it_s, ++i)
      {
        const std::vector< double >&vdata =
          inputMesh.GetProperty< AIFFace::ptr, std::vector< double > >(
            *it_s, (*itF)->GetIndex());

        if (itF == fRange.begin())
        {
          field_names.push_back("CELL_DATA_" + it_s->substr(12));
          field_attributes.resize(field_attributes.size() + 1);
        }
        if (useVertexDatafield && (it_s == dff_names.begin()))
          i = dfv_names.size();

        field_attributes[i].push_back(vdata);
      }
    }
  }

  if(has_extension(filePath, ".obj"))
  {
    IO::write_obj_file(filePath,
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
  else if(has_extension(filePath, ".off") || has_extension(filePath, ".coff"))
  {
    IO::write_off_file(filePath,
                       points_coords,
                       normals_coords,
                       texture_coords,
                       vertex_color_coords,
                       faces_indices,
                       face_color_coords,
                       helpers::num_edges(inputMesh));
  }
  else if(has_extension(filePath, ".ply"))
  {
    throw std::runtime_error(
        "Writer::write -> ply writer has not been implemented yet");
  }
  else if(has_extension(filePath, ".msh"))
  {
    IO::write_gmsh_file(filePath,
                         points_coords,
                         normals_coords,
                         vertex_color_coords,
                         lines_indices,
                         line_colors_coords,
                         faces_indices,
                         face_color_coords,
                         field_attributes,
                         field_names);
  }
#ifdef FEVV_USE_VTK
  else if(has_extension(filePath, ".vtk") || has_extension(filePath, ".vtp") ||
          has_extension(filePath, ".vtu"))
  {
    IO::write_vtk_or_vtp_or_vtu_file(filePath,
                                     points_coords,
                                     normals_coords,
                                     vertex_color_coords,
                                     lines_indices,
                                     line_colors_coords,
                                     faces_indices,
                                     face_color_coords,
                                     field_attributes,
                                     field_names);
  }
#endif

  // std::cout << "AIFWriting of \"" << get_file_full_name(filePath) << "\"
  // Done." << std::endl;
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV
