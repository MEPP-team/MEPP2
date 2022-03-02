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
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/EdgeLengthMetric.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/VolumePreserving.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/QEM3D.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Decompression/BatchDecompressor.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/RawPositions.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/DeltaPredictor.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Helpers/HeaderHandler.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Decompression/BinaryBatchDecoder.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/quantification.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/dequantization.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Butterfly.h"


namespace FEVV {
namespace Filters {

/**
 * \brief  Takes a binary file as an input, will decode the coarse mesh, the 
 * compression settings and refine the mesh until the end of file is reached. 
 * The final mesh will be saved in progressive_decompressionFilterOutput.obj
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename EdgeColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
progressive_decompression_filter(HalfedgeGraph &g,
                                 PointMap &pm,
                                 VertexColorMap &v_cm,
                                 EdgeColorMap &e_cm,
                                 const GeometryTraits &/*gt*/,
                                 std::string &input_file_path,
                                 bool dequantize)
{
  FEVV::HeaderHandler HH;
  std::fstream filebin;
  filebin.open(input_file_path, std::fstream::in | std::fstream::binary);
  draco::DecoderBuffer buffer;
  buffer.set_bitstream_version(DRACO_BITSTREAM_VERSION(2, 2));

  std::vector< char > data;
  char b;
  //Open file, load binary data
  while(filebin.get(b))
  {
    data.push_back(b);
  }
  filebin.close();

  //init draco buffer with the content of the binary file
  buffer.Init(data.data(), data.size());
  FEVV::Filters::CoarseMeshDecoder<
      HalfedgeGraph,
      PointMap,
      typename GeometryTraits::Vector,
      typename GeometryTraits::Point,
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
      typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor >
      coarse(g, pm);
  //convert draco encoded mesh to a halfedge representation
  coarse.decodeCoarseMesh(buffer);

  // retrieve decompression info (predictor, kept position, quantization 
  // settings)
  HH.DecodeBinaryHeader(buffer);


  using FaceLabelMap =
      typename FEVV::Face_pmap_traits< HalfedgeGraph, int >::pmap_type;
  FaceLabelMap f_lmap;

  f_lmap = FEVV::make_face_property_map< HalfedgeGraph, int >(g);

  //create predictor and kept position objects according to header info
  FEVV::Filters::VKEPT_POSITION _vkept = HH.getVkept();
  FEVV::Filters::PREDICTION_TYPE _pred = HH.getPred();
  FEVV::Filters::KeptPosition<
      HalfedgeGraph,
      PointMap> *KP = nullptr;
  FEVV::Filters::Predictor<
      HalfedgeGraph,
      PointMap> *predictor = nullptr;


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
    predictor = new FEVV::Filters::DeltaPredictor<
        HalfedgeGraph,
        PointMap>(g, KP, pm);
  }
  else if (_pred == FEVV::Filters::PREDICTION_TYPE::POSITION)
  {
    std::cout << "position" << std::endl;
    predictor = new FEVV::Filters::RawPositions< HalfedgeGraph,
      PointMap >(g, KP, pm);
  }
  else
  {
    delete KP;
    throw std::runtime_error("progressive_decompression_filter cannot handle geometric prediction type.");
  }

  FEVV::Filters::BatchDecompressor<
      HalfedgeGraph,
      PointMap,
      EdgeColorMap,
      VertexColorMap>
      _batch(g,
             pm,
             predictor,
             KP,
             HH,
             v_cm,
             e_cm,
             HH.getQuantization());
  //refine mesh until complete
  while(buffer.remaining_size() != 0)
  {
    _batch.decompressBinaryBatch(buffer);
  }
  if(dequantize) // if we want, dequantize attributes (it is useful to not
  {              // dequantize for tests)
    FEVV::Filters::UniformDequantization< HalfedgeGraph, PointMap > dq(
        g, pm, HH.getQuantization(), HH.getDimension(), HH.getInitCoord());

    dq.set_max_length();
    dq.set_quantization_step();
    dq.point_dequantization();
  }
}

// Helper function to simplify (syntactic sugar) the call to
// progressive_decompression_filter
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename EdgeColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
progressive_decompression_filter(HalfedgeGraph &g,
                                 PointMap &pm,
                                 VertexColorMap &v_cm,
                                 EdgeColorMap &e_cm,
                                 std::string &input_file_path,
                                 bool dequantize)

{
  GeometryTraits gt(g);
  progressive_decompression_filter(
      g, pm, v_cm, /*v_nm,*/ e_cm, gt, input_file_path, dequantize);
}
} // namespace Filters
} // namespace FEVV
