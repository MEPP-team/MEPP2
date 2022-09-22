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
#include <fstream>
#include <iomanip>
#include <list>
#include <vector>
#include <utility>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4146 26812 26451)
#endif

#include "Binary_batch_encoder_draco_nowarning.h"

#if defined _MSC_VER
#pragma warning(pop)
#endif

namespace FEVV {
namespace Filters {

/** Binary_batch_encoder class.
 *  \brief Encode binary data with draco.
 *         Also contains static methods to
 *         generate csv files.
 */
template<
    typename HalfedgeGraph,
    typename PointMap,
    typename Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector,
    typename Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point,
    typename vertex_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
    typename Geometry = FEVV::Geometry_traits< HalfedgeGraph > >
class Binary_batch_encoder
{
private:
  draco::Mesh g_draco; /// draco mesh to fill with right data
                       /// before its encoding
public:
  /// \brief Add column titles to cvs file (to call before
  /// the addition of measurement lines).
  static void init_measure_file(const std::string& csv_file_path)
  {
    std::ofstream csv_file;
    csv_file.open(csv_file_path);
    csv_file
        << "NUM_BATCH,SIZE_BITMASK,SIZE_CONNECTIVITY,SIZE_RESIDUALS,BITMASK_"
           "ENTROPY,CONNECTIVITY_ENTROPY,RESIDUALS_ENTROPY,DISTORTION_RMSE,"
           "DISTORTION_HAUSDORFF,NBVALUESBITMASK,NBVALUESCONNECTIVITY,NBVALUESOTHERINFO,"
           "SIZEOTHERINFO,ENTROPYOTHERINFO,SIZE_HEADER,NUM_VERTICES"
        << "\n";

    csv_file.close();
  }
  /// \brief Add a line of measurements to cvs file.
  static void add_line_measure_file(const std::string &csv_file_path,
                          int num_batch,
                          std::pair< int, double > &data_bitmask, // returned by Binary_batch_encoder::encode_bitmask (gives the the actual size of the bitmask and the entropy)
                          std::pair< int, double > &data_connectivity,// returned by Binary_batch_encoder::encode_bitmask
                          std::pair< int, double > &data_residuals,// returned by Binary_batch_encoder::encode_residuals
                          std::pair< int, double > &data_other_info,// returned by Binary_batch_encoder::encode_bitmask
                          double dist_RMSE, // Distortion value
                          double dist_hausdorff, // Distortion value
                          size_t num_values_bitmask, 
                          size_t num_values_connectivity,
                          size_t num_values_other_info,
                          size_t size_header, // size of encoded header
                          size_t num_vertices)
  {
    std::ofstream csv_file;

    csv_file.open(csv_file_path, std::ofstream::app);
    csv_file << num_batch << "," << data_bitmask.first << ","
             << data_connectivity.first << "," << data_residuals.first << ","
             << data_bitmask.second << "," << data_connectivity.second << ","
             << data_residuals.second << "," << dist_RMSE << "," << dist_hausdorff << ","
             << num_values_bitmask << "," << num_values_connectivity << ","
             << num_values_other_info << "," << data_other_info.first << ","
             << data_other_info.second << "," << size_header << ","
             << num_vertices << "\n";

    csv_file.close();
  }
 
  /// \brief Encodes a non-textured quantized mesh with draco single-rate
  ///        mesh encoder.
  size_t quantize_encode_coarse_mesh(HalfedgeGraph &g, /// Mesh to encode
                                     PointMap &pm, /// Vertex positions property map
                                     const Geometry &gt, /// Geometry trait object
                                     draco::EncoderBuffer &buffer /// Output encoded buffer (using draco)
								  ) 
  {
    auto iterator_pair = vertices(g);
    auto faces_pair = faces(g);
    auto num_faces_ = std::distance(faces_pair.first, faces_pair.second);
    auto num_positions_ =
        std::distance(iterator_pair.first, iterator_pair.second);
    size_t original_size = buffer.size();
    if(num_faces_ == 0)
    {
      return 0;
    }

    g_draco.SetNumFaces(num_faces_);
    draco::PointCloud *out_point_cloud_ =
        static_cast< draco::PointCloud * >(&g_draco);
    out_point_cloud_->set_num_points(
        3 * static_cast<unsigned int>(num_faces_)); 

    draco::GeometryAttribute va;
    va.Init(draco::GeometryAttribute::POSITION,
            nullptr,
            3,
            draco::DT_INT32, // \todo Why to not use uint32 instead?
            false,
            sizeof(int32_t) * 3,
            0);
    int pos_att_id_ = out_point_cloud_->AddAttribute(va, false, static_cast<unsigned int>(num_positions_));
	
    // Add vertex position coordinates
    std::map< vertex_descriptor, int > vertex_pos;

    int pos_num = 0;
    auto vi = iterator_pair.first;
    std::vector< uint32_t > vertex_data;
	vertex_data.reserve(num_positions_ * 3);
    for( ; vi != iterator_pair.second; ++vi, ++pos_num)
    {
      vertex_pos.insert(std::pair< vertex_descriptor, int >(*vi, pos_num));
      uint32_t val[3];
      const Point& pos = get(pm, *vi);
      val[0] = static_cast< uint32_t >(gt.get_x(pos));
      val[1] = static_cast< uint32_t >(gt.get_y(pos));
      val[2] = static_cast< uint32_t >(gt.get_z(pos));
      vertex_data.push_back(val[0]);
      vertex_data.push_back(val[1]);
      vertex_data.push_back(val[2]);
      out_point_cloud_->attribute(pos_att_id_)
          ->SetAttributeValue(draco::AttributeValueIndex(pos_num), val);
    }

    // Add connectivity data
    int pos_face = 0;
    auto iterator_pair_f = faces(g);
    auto fi = iterator_pair_f.first;
    for( ; fi != iterator_pair_f.second; ++fi, ++pos_face)
    {
      std::vector< vertex_descriptor > vertices_of_face;
	  vertices_of_face.reserve(3);
      boost::iterator_range<
          CGAL::Vertex_around_face_iterator< HalfedgeGraph > >
          iterator_range_vf = CGAL::vertices_around_face(halfedge(*fi, g), g);
      for(auto vf : iterator_range_vf)
        vertices_of_face.push_back(vf);
      for(int i = 0; i < 3; ++i)
      {
        const draco::PointIndex vert_id(3 * pos_face + i);
        vertex_descriptor vertex = vertices_of_face[i];
        int vertex_id = (vertex_pos.find(vertex))->second;
        out_point_cloud_->attribute(pos_att_id_)
            ->SetPointMapEntry(vert_id, draco::AttributeValueIndex(vertex_id));
      }
    }
    draco::Mesh::Face face;
    for(draco::FaceIndex i(0); i < static_cast<uint32_t>(num_faces_); ++i)
    {
      for(int c = 0; c < 3; ++c)
        face[c] = 3 * i.value() + c;
      g_draco.SetFace(i, face);
    }

    draco::Encoder encoder;
    /// desactivate draco quantization for attributes
    encoder.SetAttributeQuantization(draco::GeometryAttribute::POSITION, -1);
    const draco::Status status = encoder.EncodeMeshToBuffer(g_draco, &buffer);

    return (buffer.size() - original_size) * 8;
  }
  /// \brief Encodes a bit mask with draco's non-adaptive RAns coder 
  ///        (RAnsBitEncoder).
  std::pair< int, double > encode_bitmask(const std::list< bool > &bitmask,
                                         draco::EncoderBuffer &buffer)
  {
    draco::RAnsBitEncoder encoder;
    int size_buffer = static_cast< int >(buffer.size());
    draco::EncodeVarint((int)bitmask.size(), &buffer);
    int nb_true = 0;

    encoder.StartEncoding();
	auto it = bitmask.begin(), it_e = bitmask.end();
    for( ; it != it_e; ++it)
    {
      encoder.EncodeBit(*it);
      if(*it == true)
      {
        nb_true += 1;
      }
    }
    encoder.EndEncoding(&buffer);
	
    double bitmask_entropy = draco::ComputeBinaryShannonEntropy(
        static_cast< uint32_t >(bitmask.size()), nb_true);

    return std::make_pair(static_cast< int >((buffer.size()) - size_buffer) * 8,
                          bitmask_entropy);
  }
  /// \brief  Encodes a set of residuals with draco's entropy coder 
  ///         (SymbolBitEncoder)
  std::pair< int, double >
  encode_residuals(const std::list< std::vector< Vector > > &residuals,
                  draco::EncoderBuffer &buffer,
                  int bit_quantization)
  {
    int size_buffer = static_cast< int >(buffer.size());
    draco::EncodeVarint((int)residuals.size(), &buffer);
    if((int)residuals.size() > 0)
    {
      size_t nb_residuals = residuals.front().size(); // either one or two residuals per vertex split

      int i = 0;

      int32_t *p_in = new int32_t[residuals.size() * nb_residuals * 3];
      uint32_t *p_out = new uint32_t[residuals.size() * nb_residuals * 3];

      auto residuals_it = residuals.begin(), residuals_it_e = residuals.end();

      for( ; residuals_it != residuals_it_e; ++residuals_it, ++i)
      {
        for(size_t j = 0; j < nb_residuals; j++)
        {

          int32_t x = static_cast< int32_t >(((*residuals_it)[j])[0]);
          int32_t y = static_cast< int32_t >(((*residuals_it)[j])[1]);
          int32_t z = static_cast< int32_t >(((*residuals_it)[j])[2]);

          p_in[nb_residuals * 3 * i + 3 * j + 0] = x;
          p_in[nb_residuals * 3 * i + 3 * j + 1] = y;
          p_in[nb_residuals * 3 * i + 3 * j + 2] = z;
        }
      }
      draco::ConvertSignedIntsToSymbols(
          p_in, static_cast< int >(residuals.size() * nb_residuals * 3), p_out);
      // entropy computation
      uint32_t max = 0;
      for(size_t i = 0; i < residuals.size() * nb_residuals * 3; i++)
      {
        if(p_out[i] > max)
        {
          max = p_out[i];
        }
      }
      auto stream_2_entropy_bits = draco::ComputeShannonEntropy(
          p_out,
          static_cast< int >(residuals.size() * nb_residuals * 3),
          max,
          nullptr);


      // delete pointers in
      delete[] p_in;


      // encode the prediction
      draco::SymbolBitEncoder symbolBitEncoder;
      symbolBitEncoder.StartEncoding();
      for(size_t count = 0; count < residuals.size() * nb_residuals * 3; count++)
      {
        //auto test = p_out[count];
        symbolBitEncoder.EncodeLeastSignificantBits32(bit_quantization + 1,
                                                      p_out[count]);
      }

      symbolBitEncoder.EndEncoding(&buffer);

      size_t num_values = residuals.size() * nb_residuals * 3;
      int num_values_int = (int)num_values;
      draco::EncoderBuffer buffer_max_comp;

      draco::Options *opt = new draco::Options();
      draco::SetSymbolEncodingCompressionLevel(opt, static_cast< int >(10));
      draco::EncodeSymbols(p_out, num_values_int, 1, opt, &buffer_max_comp);
      int size_original_method =
          (static_cast< int >(buffer.size()) - size_buffer) * 8;

      delete opt;
      delete[] p_out;
      return std::make_pair(size_original_method, stream_2_entropy_bits / static_cast<double>(residuals.size() * nb_residuals * 3));
    }
    return std::make_pair(0, 0);
  }
};
} // namespace Filters
} // namespace FEVV
