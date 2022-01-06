// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
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
#include "FEVV/Wrappings/Geometry_traits.h"

#include <iostream>
#include <vector>
#include <map>
#include <utility>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4146 26812 26451)
#endif

#include "coarsemeshdecoder_draco_nowarning.h"

#if defined _MSC_VER
#pragma warning(pop)
#endif

namespace FEVV {
namespace Filters {
template<
    typename HalfedgeGraph,
    typename PointMap,
    typename Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector,
    typename Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point,
    typename vertex_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
    typename halfedge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor >
class CoarseMeshDecoder
{
public:

  /// CoarseMeshDecoder
  /// \param[in] g an empty graph
  /// \param[in] pm an empty point map associated with g
  CoarseMeshDecoder(HalfedgeGraph &g, /// empty graph
	  PointMap &pm /// associated point map
  ) : _g(g), _pm(pm) {}

  ~CoarseMeshDecoder() {}

  /// decode draco coarse mesh (vertex position and topology only)
  void decodeCoarseMesh(draco::DecoderBuffer &buffer)
  {
    draco::Decoder decoder;
    decoder.SetSkipAttributeTransform(
        draco::GeometryAttribute::POSITION); // do not touch attributes
                                             // quantization (already done)
    std::unique_ptr< draco::PointCloud > pc;

    draco::Mesh *mesh = nullptr;
	// if our mesh is a triangular mesh, decode it
    auto type_statusor = decoder.GetEncodedGeometryType(&buffer);
    if(!type_statusor.ok())
    {
      std::cerr << "decodeCoarseMesh: failed to retrieve geometry " << std::endl;
      return;
    }
    const draco::EncodedGeometryType geom_type = type_statusor.value();
    if(geom_type == draco::TRIANGULAR_MESH)
    {
      auto statusor = decoder.DecodeMeshFromBuffer(&buffer);
      if(!statusor.ok())
      {
        std::cerr << "decodeCoarseMesh: failed to retrieve coarse mesh " << std::endl;
        return;
      }
      std::unique_ptr< draco::Mesh > in_mesh = std::move(statusor).value();
      if(in_mesh)
      {
        mesh = in_mesh.get();
        pc = std::move(in_mesh);
      }
    }

    if(pc == nullptr){
      std::cerr << "decodeCoarseMesh: Failed to decode the input file" << std::endl;
    }
	else
      getMeshFromBinaryFile(*mesh);
  }
 
  bool getMeshFromBinaryFile(draco::Mesh &mesh)
  {
    const draco::PointCloud *p_point_cloud =
        &(static_cast< const draco::PointCloud & >(mesh));
    const draco::Mesh *p_mesh = &mesh;
    // create a map that permits retriving a vertex from its index
    std::map< int, vertex_descriptor > map_index_to_vertex;

    // Encode points position into a buffer and faces into an other buffer
    const draco::PointAttribute *const att =
        p_point_cloud->GetNamedAttribute(draco::GeometryAttribute::POSITION);
    if(!GetMeppMeshVertex(p_point_cloud, att, &map_index_to_vertex))
      return false;
    if(!GetMeppMeshFace(p_mesh, att, map_index_to_vertex))
      return false;
    return true;
  }

  bool
  GetMeppMeshVertex(const draco::PointCloud * /*p_point_cloud*/,
                    const draco::PointAttribute *const att,
                    std::map< int, vertex_descriptor > *map_index_to_vertex)
  {
    if(att == nullptr || att->size() == 0)
      return false; // Position attribute must be valid.
    std::array< float, 3 > value;
    for(draco::AttributeValueIndex i(0);
        i < static_cast< uint32_t >(att->size());
        ++i)
    {
      if(!att->ConvertValue< float, 3 >(i, &value[0]))
        return false;
      // insert new vertex in the mesh
      vertex_descriptor v = add_vertex(_g);
      map_index_to_vertex->insert(
          std::pair< int, vertex_descriptor >(i.value(), v));
      // get vertex position
      Point vertex_pos = Point(value[0], value[1], value[2]);
      // insert vertex position in map pm
      put(_pm, v, vertex_pos);
    }
    return true;
  }

  bool
  GetMeppMeshFace(const draco::Mesh *p_mesh,
                  const draco::PointAttribute *const att,
                  const std::map< int, vertex_descriptor > &map_index_to_vertex)
  {
    for(draco::FaceIndex i(0); i < p_mesh->num_faces(); ++i)
    {
      std::vector< int > face;
	  face.reserve(3);
      for(int j = 0; j < 3; ++j)
      {
        if(!GetFaceCorner(i, j, p_mesh, face, att))
          return false;
      }

      // update mesh connectivity
      setMeshConnectivityFromFace(face, map_index_to_vertex);
    }
    return true;
  }

  bool GetFaceCorner(draco::FaceIndex face_id,
                     int local_corner_id,
                     const draco::Mesh *p_mesh,
                     std::vector< int >& face,
                     const draco::PointAttribute *const att)
  {
    const draco::PointIndex vert_index = p_mesh->face(face_id)[local_corner_id];
    // get vertices that makes the face
    face.push_back(att->mapped_index(vert_index).value());
    return true;
  }


  void setMeshConnectivityFromFace(
      const std::vector< int >& face,
      const std::map< int, vertex_descriptor > &map_index_to_vertex)
  {
    // make a vertex range from the face vertices
    std::vector< vertex_descriptor > vertexRange;
	vertexRange.reserve(3);
    for(int i = 0; i < 3; i++)
    {
      vertex_descriptor v = map_index_to_vertex.find(face[i])->second;
      vertexRange.push_back(v);
    }

    // add face
    CGAL::Euler::add_face< HalfedgeGraph, std::vector< vertex_descriptor > >(
        vertexRange, _g);
  }

private:
  HalfedgeGraph &_g; // HalfedgeGraph: at the beginning, will be empty. the 
                     // function decodeCoarseMesh will
                     // add the vertices and corresponding connectivity from a 
                     // draco buffer
  PointMap &_pm; // Empty point map at the beginning. It will be filled thanks
                 // to the same draco buffer as _g
};


} // namespace Filters
} // namespace FEVV
