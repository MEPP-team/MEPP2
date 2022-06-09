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

#include "FEVV/Wrappings/Geometry_traits.h"
#include <boost/graph/graph_traits.hpp>
#include "FEVV/Wrappings/properties.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Parameters.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/Edge_length_metric.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/Volume_preserving.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/QEM_3D.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Decompression/Batch_decompressor.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Raw_positions.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Delta_predictor.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Helpers/Header_handler.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Decompression/Binary_batch_decoder.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Uniform_quantization.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Uniform_dequantization.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Butterfly.h"

//#define TIMING
#ifdef TIMING
#include <chrono>
#endif

namespace FEVV {
namespace Filters {

  /**
   * \brief Takes a buffer as an input, will decode the compression settings, 
   *    the coarse mesh and refine the mesh until the end of buffer is reached.
   *    This function implements the Algorithm 2 (Progressive decoder).
   *
   * @param[out] g Empty halfedge graph that is used to reconstruct the
   *                decoded mesh topology.
   * @param[out] pm Empty pointmap associated with g that is used to
   *                reconstruct the decoded mesh geometry (vertex positions).
   * @param[out] v_cm Empty vertex color map associated with g, used to debug
   *                topological issues during the decoding.
   * @param[in] gt The geometry trait object. Used to automatically find
   *            the GeometryTraits type.
   * @param[in] buffer The draco buffer to decode. Not a const param because
   *            of call to non-const draco methods such as Decode.
   * @param[in] dequantize True to dequantize the vertex positions at
   *            the end of the decompression, else vertex positions remain
   *            quantized, used to test (without dequantization).
   */
template< typename HalfedgeGraph,
    typename PointMap,
    typename VertexColorMap,
    typename GeometryTraits >
void 
progressive_decompression_filter(HalfedgeGraph& g,
                                 PointMap& pm,
                                 VertexColorMap& v_cm,
                                 const GeometryTraits&/*gt*/,
                                 draco::DecoderBuffer& buffer,
                                 bool dequantize) 
{
    ////////////////////////
    // DECODE FILE HEADER //
    ////////////////////////
    // Retrieve decompression info (predictor, kept position, quantization 
    // settings)
    FEVV::Header_handler HH;
    HH.decode_binary_header(buffer);

    // Retrieve the base mesh
    FEVV::Filters::Coarse_mesh_decoder<
      HalfedgeGraph,
      PointMap,
      typename GeometryTraits::Vector,
      typename GeometryTraits::Point,
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
      typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor >
      coarse(g, pm);
    coarse.decode_coarse_mesh(buffer);

    // Create predictor and kept position objects according to header info
    FEVV::Filters::VKEPT_POSITION _vkept = HH.get_vkept();
    FEVV::Filters::PREDICTION_TYPE _pred = HH.get_pred();
    FEVV::Filters::Kept_position<
      HalfedgeGraph,
      PointMap>* KP = nullptr;
    FEVV::Filters::Predictor<
      HalfedgeGraph,
      PointMap>* predictor = nullptr;


    if(_vkept == FEVV::Filters::VKEPT_POSITION::HALFEDGE)
    {
      KP = new FEVV::Filters::Halfedge<
        HalfedgeGraph,
        PointMap >(g, pm);
    }
    else if(_vkept == FEVV::Filters::VKEPT_POSITION::MIDPOINT)
    {
      KP = new FEVV::Filters::Midpoint<
        HalfedgeGraph,
        PointMap >(g, pm);
    }
    else
      throw std::runtime_error("progressive_decompression_filter cannot handle kept position type.");


    if(_pred == FEVV::Filters::PREDICTION_TYPE::BUTTERFLY)
    {
      std::cout << "butterfly" << std::endl;
      predictor = new FEVV::Filters::Butterfly<
        HalfedgeGraph,
        PointMap>(g, KP, pm);
    }
    else if(_pred == FEVV::Filters::PREDICTION_TYPE::DELTA)
    {
      std::cout << "delta" << std::endl;
      predictor = new FEVV::Filters::Delta_predictor<
        HalfedgeGraph,
        PointMap>(g, KP, pm);
    }
    else if(_pred == FEVV::Filters::PREDICTION_TYPE::POSITION)
    {
      std::cout << "position" << std::endl;
      predictor = new FEVV::Filters::Raw_positions< HalfedgeGraph,
        PointMap >(g, KP, pm);
    }
    else
    {
      delete KP;
      throw std::runtime_error("progressive_decompression_filter cannot handle geometric prediction type.");
    }

    ///////////////////////////////
    // DECODE REFINEMENT BATCHES //
    ///////////////////////////////
    // Create a batch decompressor object to decode the refinement information
    // of the following LoDs and update consequently the halfedge graph g and 
    // its point map pm.
    FEVV::Filters::Batch_decompressor<
      HalfedgeGraph,
      PointMap,
      VertexColorMap>
      _batch(g,
        pm,
        predictor,
        KP,
        HH,
        v_cm);

    // Refine mesh until the buffer is fully decoded
    // Implements the while loop of Algorithm 2.
    while(buffer.remaining_size() != 0)
    {
      _batch.decompress_binary_batch(buffer);
    }
    if(dequantize) // if we want, dequantize attributes (it is useful to not
    {              // dequantize for tests)
      // These lines implement line 25 of Algorithm 2.
      FEVV::Filters::Uniform_dequantization< HalfedgeGraph, PointMap > dq(
        g, pm, HH.get_quantization(), HH.get_dimension(), HH.get_init_coord());
      dq.point_dequantization();
    }
  }

/**
 * \brief Takes a binary file as an input, will decode the compression  
 *   settings, the coarse mesh and refine the mesh until the end of file is 
 *   reached. 
 *
 * @param[out] g Empty halfedge graph that is used to reconstruct the
 *                decoded mesh topology.
 * @param[out] pm Empty pointmap associated with g that is used to 
 *                reconstruct the decoded mesh geometry (vertex positions).
 * @param[out] v_cm Empty vertex color map associated with g, used to debug
 *                topological issues during the decoding.
 * @param[in] gt The geometry trait object. Used to automatically find
 *            the GeometryTraits type.
 * @param[in] input_file_path The binary file full name to read and decode.
 * @param[in] dequantize True to dequantize the vertex positions at
 *            the end of the decompression, else vertex positions remain 
 *            quantized, used to test (without dequantization).
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename GeometryTraits >
void
progressive_decompression_filter(HalfedgeGraph &g,
                                 PointMap &pm,
                                 VertexColorMap &v_cm,
                                 const GeometryTraits & gt,
                                 const std::string &input_file_path,
                                 bool dequantize)
{
#ifdef TIMING
  auto time_point_before_decomp = std::chrono::high_resolution_clock::now();
#endif
  std::fstream filebin;
  filebin.open(input_file_path, std::fstream::in | std::fstream::binary);
  draco::DecoderBuffer buffer;
  buffer.set_bitstream_version(DRACO_BITSTREAM_VERSION(2, 2));

  std::vector< char > data;
  char b;
  // Open file, load binary data.
  while(filebin.get(b))
  {
    data.push_back(b);
  }
  filebin.close();

  // Init draco buffer with the content of the binary file.
  buffer.Init(data.data(), data.size());
  
  progressive_decompression_filter(g, pm, v_cm, gt, buffer, dequantize);

#ifdef TIMING
  auto time_point_after_decomp = std::chrono::high_resolution_clock::now();
  auto duration_decomp = std::chrono::duration_cast<std::chrono::milliseconds>(time_point_after_decomp - time_point_before_decomp);
  std::cout << "Decompression time: " << duration_decomp.count() << " milliseconds" << std::endl;
#endif
}

/**
 * \brief Takes a binary file as an input, will decode the compression
 *   settings, the coarse mesh and refine the mesh until the end of file is
 *   reached. The geometry trait object is set automatically (syntactic sugar).
 *
 * @param[out] g Empty halfedge graph that is used to reconstruct the
 *                decoded mesh topology.
 * @param[out] pm Empty pointmap associated with g that is used to
 *                reconstruct the decoded mesh geometry (vertex positions).
 * @param[out] v_cm Empty vertex color map associated with g, used to debug
 *                topological issues during the decoding.
 * @param[in] input_file_path The binary file full name to read and decode.
 * @param[in] dequantize True to dequantize the vertex positions at
 *            the end of the decompression, else vertex positions remain
 *            quantized, used to test (without dequantization).
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
progressive_decompression_filter(HalfedgeGraph &g,
                                 PointMap &pm,
                                 VertexColorMap &v_cm,
                                 const std::string &input_file_path,
                                 bool dequantize)

{
  GeometryTraits gt(g);
  progressive_decompression_filter< HalfedgeGraph, 
                                    PointMap, 
                                    VertexColorMap, 
                                    GeometryTraits>(
      g, pm, v_cm, gt, input_file_path, dequantize);
}

} // namespace Filters
} // namespace FEVV
