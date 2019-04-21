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

// FBX SDK
#include <fbxsdk.h>

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "FEVV/Tools/IO/StringUtilities.hpp"
#include "FEVV/Types/Material.h"


/*
 * Doc : http://docs.autodesk.com/FBX/2014/ENU/FBX-SDK-Documentation/
 */

namespace FEVV {
namespace IO {


using namespace StrUtils;
using namespace FileUtils;


/**
 * Return a string-based representation based on the attribute type.
 */
inline FbxString
dbg_get_attr_type_name(FbxNodeAttribute::EType type)
{
  switch(type)
  {
  case FbxNodeAttribute::eUnknown:
    return "unidentified";
  case FbxNodeAttribute::eNull:
    return "null";
  case FbxNodeAttribute::eMarker:
    return "marker";
  case FbxNodeAttribute::eSkeleton:
    return "skeleton";
  case FbxNodeAttribute::eMesh:
    return "mesh";
  case FbxNodeAttribute::eNurbs:
    return "nurbs";
  case FbxNodeAttribute::ePatch:
    return "patch";
  case FbxNodeAttribute::eCamera:
    return "camera";
  case FbxNodeAttribute::eCameraStereo:
    return "stereo";
  case FbxNodeAttribute::eCameraSwitcher:
    return "camera switcher";
  case FbxNodeAttribute::eLight:
    return "light";
  case FbxNodeAttribute::eOpticalReference:
    return "optical reference";
  case FbxNodeAttribute::eOpticalMarker:
    return "marker";
  case FbxNodeAttribute::eNurbsCurve:
    return "nurbs curve";
  case FbxNodeAttribute::eTrimNurbsSurface:
    return "trim nurbs surface";
  case FbxNodeAttribute::eBoundary:
    return "boundary";
  case FbxNodeAttribute::eNurbsSurface:
    return "nurbs surface";
  case FbxNodeAttribute::eShape:
    return "shape";
  case FbxNodeAttribute::eLODGroup:
    return "lodgroup";
  case FbxNodeAttribute::eSubDiv:
    return "subdiv";
  default:
    return "unknown";
  }
}

/**
 * Print an attribute.
 */
inline void
dbg_print_attribute(FbxNodeAttribute *pAttribute)
{
  if(!pAttribute)
    return;

  FbxString typeName = dbg_get_attr_type_name(pAttribute->GetAttributeType());
  FbxString attrName = pAttribute->GetName();

  // Note: to retrieve the character array of a FbxString, use its Buffer()
  // method.
  std::cout << '\n'
            << typeName.Buffer() << ' ' << attrName.Buffer() << std::endl;
}

/**
 * Print a node, its attributes, and all its children recursively.
 */
inline void
dbg_print_node(FbxNode *pNode)
{
  const auto nodeName = pNode->GetName();
  FbxDouble3 translation = pNode->LclTranslation.Get();
  FbxDouble3 rotation = pNode->LclRotation.Get();
  FbxDouble3 scaling = pNode->LclScaling.Get();

  // Print the contents of the node.
  std::cout << nodeName << ' ' << translation[0] << ' ' << translation[1] << ' '
            << translation[2] << ' ' << rotation[0] << ' ' << rotation[1] << ' '
            << rotation[2] << ' ' << scaling[0] << ' ' << scaling[1] << ' '
            << scaling[2];

  // Print the node's attributes.
  for(int i = 0; i < pNode->GetNodeAttributeCount(); ++i)
    dbg_print_attribute(pNode->GetNodeAttributeByIndex(i));

  // Recursively print the children.
  for(int i = 0; i < pNode->GetChildCount(); ++i)
    dbg_print_node(pNode->GetChild(i));

  std::cout << std::endl;
}

inline void
dbg_display_infos(const FbxScene &scene)
{
  // Print the nodes of the scene and their attributes recursively.
  // Note that we are not printing the root node because it should
  // not contain any attributes.
  FbxNode *root_node = scene.GetRootNode();
  if(root_node)
  {
    for(int i = 0; i < root_node->GetChildCount(); ++i)
      dbg_print_node(root_node->GetChild(i));
  }
}


/**
 * \brief  Read the FBX file and load its data into given parameters.
 *
 * \param  filePath              path to the FBX file.
 * \param  points_coords         list of positions' data.
 * \param  normals_coords        list of normals' data.
 * \param  texture_coords        list of texture coordinates' data.
 * \param  vertex_color_coords   list of vertex colors' data.
 * \param  face_indices          list of face indices.
 * \param  texture_face_indices  list of texture coordinates' indices.
 * \param  normal_face_indices   list of normals' indices.
 * \param  materials             list of all the  materials the file contains.
 * \param  face_material         list of material index per face.
 */
template< typename coord_type,
          typename coordN_type,
          typename coordT_type,
          typename coordC_type,
          typename index_type,
          typename material_type >
void
read_fbx_file(const std::string &filePath,
              std::vector< std::vector< coord_type > > &points_coords,
              std::vector< std::vector< coordN_type > > &normals_coords,
              std::vector< std::vector< coordT_type > > &texture_coords,
              std::vector< std::vector< coordC_type > > &vertex_color_coords,
              std::vector< std::vector< index_type > > &face_indices,
              std::vector< std::vector< index_type > > &texture_face_indices,
              std::vector< std::vector< index_type > > &normal_face_indices,
              std::vector< material_type > &materials,  // list of materials
              std::vector< index_type > &face_material) // material of each face
{
  FbxManager *manager = FbxManager::Create();

  FbxIOSettings *ioSettings = FbxIOSettings::Create(manager, IOSROOT);
  manager->SetIOSettings(ioSettings);

  FbxImporter *importer = FbxImporter::Create(manager, "");

  if(!importer->Initialize(filePath.c_str(), -1, manager->GetIOSettings()))
    throw std::invalid_argument(
        {"Reader::read_fbx_file -> Error: couldn't load the given FBX file.\n",
         importer->GetStatus().GetErrorString()});

  FbxScene *scene =
      FbxScene::Create(manager, FileUtils::get_file_name(filePath).c_str());

  // Importing the contents of the file into the scene
  importer->Import(scene);
  importer->Destroy();

  // RM: may be needed if the mesh is not oriented correctly
  FbxAxisSystem::OpenGL.ConvertScene(scene);

  // dbg_display_infos(*scene);

  auto &settings = scene->GetGlobalSettings();

  // Recovering geometry
  for(int meshIndex = 0; meshIndex < scene->GetGeometryCount(); ++meshIndex)
  {
    const auto mesh = static_cast< FbxMesh * >(scene->GetGeometry(meshIndex));

    std::cout << "[FbxFileReader] --- Mesh no." << meshIndex
              << "'s poly count: " << mesh->GetPolygonCount() << std::endl;

    ////////////
    // Values //
    ////////////

    // Recovering positions
    std::size_t currentVertIndex = points_coords.size();
    points_coords.resize(currentVertIndex + mesh->GetControlPointsCount(),
                         std::vector< coord_type >(3));

    for(int vertIndex = 0; vertIndex < mesh->GetControlPointsCount();
        ++vertIndex, ++currentVertIndex)
    {
      const auto &vertPos = mesh->GetControlPointAt(vertIndex);

      points_coords[currentVertIndex][0] = vertPos[0]; // X position
      points_coords[currentVertIndex][1] = vertPos[1]; // Y position
      points_coords[currentVertIndex][2] = vertPos[2]; // Z position
    }

    // Recovering normals
    const auto &meshNormals = mesh->GetElementNormal();

    if(meshNormals)
    {
      std::size_t currentNormIndex = normals_coords.size();
      normals_coords.resize(currentNormIndex +
                                meshNormals->GetDirectArray().GetCount(),
                            std::vector< coordN_type >(3));

      for(int normIndex = 0;
          normIndex < meshNormals->GetDirectArray().GetCount();
          ++normIndex, ++currentNormIndex)
      {
        const auto &normal = meshNormals->GetDirectArray()[normIndex];

        normals_coords[currentNormIndex][0] = normal[0]; // X
        normals_coords[currentNormIndex][1] = normal[1]; // Y
        normals_coords[currentNormIndex][2] = normal[2]; // Z
      }
    }

    // Recovering texture coordinates (UVs)
    const auto &meshTexcoords = mesh->GetElementUV();

    if(meshTexcoords)
    {
      std::size_t currentTexIndex = texture_coords.size();
      texture_coords.resize(currentTexIndex +
                                meshTexcoords->GetDirectArray().GetCount(),
                            std::vector< coordT_type >(2));

      for(int texIndex = 0;
          texIndex < meshTexcoords->GetDirectArray().GetCount();
          ++texIndex, ++currentTexIndex)
      {
        const auto &texcoords = meshTexcoords->GetDirectArray()[texIndex];

        texture_coords[currentTexIndex][0] =  static_cast<coordT_type>(texcoords[0]); // U
        texture_coords[currentTexIndex][1] = static_cast<coordT_type>(texcoords[1]); // V
      }
    }

    /////////////
    // Indices //
    /////////////

    // Fetching positions' indices

    // We need to make submeshes contiguous in memory; FBX meshes' indices are
    // starting back from 0 for every mesh We then calculate the base position
    // index as a stride, to handle contiguous data in memory. Adding it to
    // indices makes them absolute For example, the first index of the second
    // mesh will be equal to the amount of entries for the first mesh
    const std::size_t basePosIndex =
        points_coords.size() - mesh->GetControlPointsCount();
    std::size_t currentPosPolyIndex = face_indices.size();
    face_indices.resize(currentPosPolyIndex + mesh->GetPolygonCount());

    for(int polyIndex = 0; polyIndex < mesh->GetPolygonCount();
        ++polyIndex, ++currentPosPolyIndex)
    {
      const auto polySize = mesh->GetPolygonSize(polyIndex);
      face_indices[currentPosPolyIndex].resize(polySize);

      for(int polyVertIndex = 0; polyVertIndex < polySize; ++polyVertIndex)
        face_indices[currentPosPolyIndex][polyVertIndex] =
            basePosIndex + mesh->GetPolygonVertex(polyIndex, polyVertIndex);
    }

    // Fetching texture coordinates' indices
    if(meshTexcoords)
    {
      if(meshTexcoords->GetMappingMode() ==
         FbxLayerElement::EMappingMode::eByControlPoint)
      {
        std::cout
            << "[FbxFileReader] Mapping mesh's texcoords by control point."
            << std::endl;

        texture_face_indices.resize(texture_face_indices.size() +
                                    mesh->GetPolygonCount());

        // Getting texcoords for each vertex, since the mapping mode of
        // texcoords element is by control point
        for(std::size_t faceIndex =
                face_indices.size() - mesh->GetPolygonCount();
            faceIndex < face_indices.size();
            ++faceIndex)
        {
          texture_face_indices[faceIndex].resize(
              face_indices[faceIndex].size());

          for(int vertIndex = 0; vertIndex < face_indices[faceIndex].size();
              ++vertIndex)
          {
            index_type texIndex{};

            if(meshTexcoords->GetReferenceMode() ==
               FbxGeometryElement::eDirect) // Texcoords index is the same as
                                            // vertex index
              texIndex = face_indices[faceIndex][vertIndex];
            else if(meshTexcoords->GetReferenceMode() ==
                    FbxGeometryElement::eIndexToDirect) // Get texcoords index
                                                        // by polygon vertex
              texIndex = meshTexcoords->GetIndexArray().GetAt(
                  static_cast< int >(face_indices[faceIndex][vertIndex]));

            texture_face_indices[faceIndex][vertIndex] = texIndex;
          }
        }
      }
      else if(meshTexcoords->GetMappingMode() ==
              FbxLayerElement::EMappingMode::eByPolygonVertex)
      {
        std::cout
            << "[FbxFileReader] Mapping mesh's texcoords by face vertices."
            << std::endl;

        const std::size_t baseTexIndex =
            texture_coords.size() - meshTexcoords->GetDirectArray().GetCount();
        std::size_t globalTexIndex = texture_face_indices.size();
        texture_face_indices.resize(globalTexIndex + mesh->GetPolygonCount());
        int currentTexPolyIndex = 0;

        // Getting texcoords of each polygon, since the mapping mode of
        // texcoords element is by polygon-vertex
        for(int polyIndex = 0; polyIndex < mesh->GetPolygonCount();
            ++polyIndex, ++globalTexIndex)
        {
          texture_face_indices[globalTexIndex].resize(
              mesh->GetPolygonSize(polyIndex));

          // Retrieve each vertex of current polygon
          for(int polyVertIndex = 0;
              polyVertIndex < mesh->GetPolygonSize(polyIndex);
              ++polyVertIndex, ++currentTexPolyIndex)
          {
            index_type texIndex{};

            if(meshTexcoords->GetReferenceMode() == FbxGeometryElement::eDirect)
              texIndex = currentTexPolyIndex;
            else if(meshTexcoords->GetReferenceMode() ==
                    FbxGeometryElement::eIndexToDirect)
              texIndex =
                  meshTexcoords->GetIndexArray().GetAt(currentTexPolyIndex);

            texture_face_indices[globalTexIndex][polyVertIndex] =
                baseTexIndex + texIndex;
          }
        }
      }
      else
      {
        std::cerr
            << "[FbxFileReader] Couldn't handle mesh's texcoords' mapping mode."
            << std::endl;
      }
    }

    if(meshNormals)
    {
      if(meshNormals->GetMappingMode() ==
         FbxLayerElement::EMappingMode::eByControlPoint)
      {
        std::cout << "[FbxFileReader] Mapping mesh's normals by control point."
                  << std::endl;

        normal_face_indices.resize(normal_face_indices.size() +
                                   mesh->GetPolygonCount());

        // Getting normals for each vertex, since the mapping mode of normal
        // element is by control point
        for(std::size_t faceIndex =
                face_indices.size() - mesh->GetPolygonCount();
            faceIndex < face_indices.size();
            ++faceIndex)
        {
          normal_face_indices[faceIndex].resize(face_indices[faceIndex].size());

          for(int vertIndex = 0; vertIndex < face_indices[faceIndex].size();
              ++vertIndex)
          {
            index_type normIndex{};

            if(meshNormals->GetReferenceMode() ==
               FbxGeometryElement::eDirect) // Normal index is the same as
                                            // vertex index
              normIndex = face_indices[faceIndex][vertIndex];
            else if(meshNormals->GetReferenceMode() ==
                    FbxGeometryElement::eIndexToDirect) // Get normal index by
                                                        // polygon vertex
              normIndex = meshNormals->GetIndexArray().GetAt(
                  static_cast< int >(face_indices[faceIndex][vertIndex]));

            normal_face_indices[faceIndex][vertIndex] = normIndex;
          }
        }
      }
      else if(meshNormals->GetMappingMode() ==
              FbxLayerElement::EMappingMode::eByPolygonVertex)
      {
        std::cout << "[FbxFileReader] Mapping mesh's normals by face vertices."
                  << std::endl;

        const std::size_t baseNormIndex =
            normals_coords.size() - meshNormals->GetDirectArray().GetCount();
        std::size_t globalNormIndex = normal_face_indices.size();
        normal_face_indices.resize(globalNormIndex + mesh->GetPolygonCount());
        int currentNormPolyIndex = 0;

        // Getting normals of each polygon, since the mapping mode of normal
        // element is by polygon-vertex
        for(int polyIndex = 0; polyIndex < mesh->GetPolygonCount();
            ++polyIndex, ++globalNormIndex)
        {
          normal_face_indices[globalNormIndex].resize(
              mesh->GetPolygonSize(polyIndex));

          // Retrieve each vertex of current polygon
          for(int polyVertIndex = 0;
              polyVertIndex < mesh->GetPolygonSize(polyIndex);
              ++polyVertIndex, ++currentNormPolyIndex)
          {
            index_type normIndex{};

            if(meshNormals->GetReferenceMode() == FbxGeometryElement::eDirect)
              normIndex = currentNormPolyIndex;
            else if(meshNormals->GetReferenceMode() ==
                    FbxGeometryElement::eIndexToDirect)
              normIndex =
                  meshNormals->GetIndexArray().GetAt(currentNormPolyIndex);

            normal_face_indices[globalNormIndex][polyVertIndex] =
                baseNormIndex + normIndex;
          }
        }
      }
      else
      {
        std::cerr
            << "[FbxFileReader] Couldn't handle mesh's normals' mapping mode."
            << std::endl;
      }
    }

    const auto &meshMaterial = mesh->GetElementMaterial();
    if(meshMaterial)
    {
      if(meshMaterial->GetMappingMode() ==
         FbxLayerElement::EMappingMode::eByPolygon)
      {
        face_material.resize(mesh->GetPolygonCount());

        for(int polyIndex = 0; polyIndex < mesh->GetPolygonCount(); ++polyIndex)
          face_material[polyIndex] = meshMaterial->GetIndexArray()[polyIndex];
      }
      else if(meshMaterial->GetMappingMode() ==
              FbxLayerElement::EMappingMode::eAllSame)
      {
        const std::size_t prevFaceMatSize = face_material.size();

        face_material.resize(prevFaceMatSize + mesh->GetPolygonCount());
        // REALLY bad hack here: without this std::min, that works for meshes
        // with the same amount of materials and meshes Clamping here
        // temporarily prevents segfaulting, but obviously isn't a long-term
        // solution
        std::fill(face_material.begin() + prevFaceMatSize,
                  face_material.end(),
                  std::min(meshIndex, scene->GetMaterialCount() - 1));
      }
    }
  }

  // Recovering materials
  materials.resize(scene->GetMaterialCount());

  for(int matIndex = 0; matIndex < scene->GetMaterialCount(); ++matIndex)
  {
    const FbxSurfaceMaterial *material = scene->GetMaterial(matIndex);

    materials[matIndex].name = material->GetName();

    // Recovering properties from the material

    const FbxPropertyT< FbxDouble3 > &ambient =
        material->FindProperty(FbxSurfaceMaterial::sAmbient);
    if(ambient.IsValid())
    {
      materials[matIndex].ambient_red_component = ambient.Get()[0];
      materials[matIndex].ambient_green_component = ambient.Get()[1];
      materials[matIndex].ambient_blue_component = ambient.Get()[2];
    }

    const FbxPropertyT< FbxDouble3 > &diffuse =
        material->FindProperty(FbxSurfaceMaterial::sDiffuse);
    if(diffuse.IsValid())
    {
      materials[matIndex].diffuse_red_component = diffuse.Get()[0];
      materials[matIndex].diffuse_green_component = diffuse.Get()[1];
      materials[matIndex].diffuse_blue_component = diffuse.Get()[2];
    }

    const FbxPropertyT< FbxDouble3 > &specular =
        material->FindProperty(FbxSurfaceMaterial::sSpecular);
    if(specular.IsValid())
    {
      materials[matIndex].specular_red_component = specular.Get()[0];
      materials[matIndex].specular_green_component = specular.Get()[1];
      materials[matIndex].specular_blue_component = specular.Get()[2];
    }

    const FbxPropertyT< FbxDouble3 > &emissive =
        material->FindProperty(FbxSurfaceMaterial::sEmissive);
    if(emissive.IsValid())
    {
      materials[matIndex].emissive_red_component = emissive.Get()[0];
      materials[matIndex].emissive_green_component = emissive.Get()[1];
      materials[matIndex].emissive_blue_component = emissive.Get()[2];
    }

    const FbxPropertyT< FbxDouble > &transparency =
        material->FindProperty(FbxSurfaceMaterial::sTransparencyFactor);
    if(transparency.IsValid())
      materials[matIndex].transparency = transparency.Get();

    // Recovering textures associated to the material
    const auto parent_directory = get_parent_directory(filePath);

    const auto ambientTexture = static_cast< FbxFileTexture * >(
        ambient.GetSrcObject(FbxCriteria::ObjectType(FbxFileTexture::ClassId)));
    if(ambientTexture)
    {
      materials[matIndex].ambient_texture_filename =
          ambientTexture->GetRelativeFileName();
      materials[matIndex].ambient_texture_filename.insert(
          0, parent_directory + "/");
    }

    const auto diffuseTexture = static_cast< FbxFileTexture * >(
        diffuse.GetSrcObject(FbxCriteria::ObjectType(FbxFileTexture::ClassId)));
    if(diffuseTexture)
    {
      materials[matIndex].diffuse_texture_filename =
          diffuseTexture->GetRelativeFileName();
      materials[matIndex].diffuse_texture_filename.insert(
          0, parent_directory + "/");
    }

    const auto specularTexture =
        static_cast< FbxFileTexture * >(specular.GetSrcObject(
            FbxCriteria::ObjectType(FbxFileTexture::ClassId)));
    if(specularTexture)
    {
      materials[matIndex].specular_texture_filename =
          specularTexture->GetRelativeFileName();
      materials[matIndex].specular_texture_filename.insert(
          0, parent_directory + "/");
    }

    const auto emissiveTexture =
        static_cast< FbxFileTexture * >(emissive.GetSrcObject(
            FbxCriteria::ObjectType(FbxFileTexture::ClassId)));
    if(emissiveTexture)
    {
      materials[matIndex].emissive_texture_filename =
          emissiveTexture->GetRelativeFileName();
      materials[matIndex].emissive_texture_filename.insert(
          0, parent_directory + "/");
    }

    const auto normalMapProp =
        material->FindProperty(FbxSurfaceMaterial::sNormalMap);
    if(normalMapProp.IsValid())
    {
      const auto normalMap =
          static_cast< FbxFileTexture * >(normalMapProp.GetSrcObject(
              FbxCriteria::ObjectType(FbxFileTexture::ClassId)));
      if(normalMap)
      {
        materials[matIndex].normal_map_filename =
            normalMap->GetRelativeFileName();
        materials[matIndex].normal_map_filename.insert(0,
                                                       parent_directory + "/");

        materials[matIndex].has_normal_map = true;
      }
    }
  }

  manager->Destroy(); // Finally destroying the FbxManager now that we've
                      // imported our scene
}

} // namespace IO
} // namespace FEVV

