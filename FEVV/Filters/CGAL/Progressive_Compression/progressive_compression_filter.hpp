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
#include <CGAL/boost/graph/helpers.h> // is_triangle_mesh
#include "FEVV/Filters/CGAL/Progressive_Compression/Parameters.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/Edge_length_metric.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/Volume_preserving.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/QEM_3D.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Batch_collapser.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Raw_positions.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Delta_predictor.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Helpers/Header_handler.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Binary_batch_encoder.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Geometric_metrics.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Uniform_quantization.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Uniform_dequantization.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Preprocessing.h"
#include "FEVV/Filters/Generic/minmax_map.h"
#include "FEVV/Filters/Generic/color_mesh.h"

//#define TIMING
#ifdef TIMING
#include <chrono>
#endif

namespace FEVV {
namespace Filters {

/**
 * \brief Isolated vertices and isolated edges are removed, then vertex 
 *        positions are quantized, and eventually, vertex duplicates are moved 
 *        until there is no more vertex position duplicate.
 *        This function implements line 1 of Algorithm 1.
 *
 * @param[in,out] g Halfedge graph that is processed.  
 * @param[in,out] pm Pointmap associated with g. The input positions are
 *                not quantized, while the ouput positions are.
 * @param[in,out] pq The Uniform_quantization object to apply a XYZ
 *                uniform quantization to points in pm.
 * @param[in] allow_duplicates The boolean set to false to move vertex 
 *            duplicates until there is no more vertex position duplicate.
 */
template< typename HalfedgeGraph,
          typename PointMap>
void
preprocess_mesh(
    HalfedgeGraph &g,
    PointMap &pm,
    FEVV::Filters::Uniform_quantization< HalfedgeGraph, PointMap > &pq,
    bool allow_duplicates = false)
{
  FEVV::Filters::Preprocessing< HalfedgeGraph, PointMap >
      preprocess(g, pm);
  preprocess.process_mesh_before_quantization();

  // vertex positions part
  pq.point_quantization(true);
 
  if(!allow_duplicates)
    preprocess.process_mesh_after_quantization();
}

/**
 * \brief Takes a mesh g, applies batches of simplification until either the
 * number of max batches or the minimum number of vertices is reached. After
 * that, it encodes the coarse mesh and the refinement data into a buffer.
 *        This function implements Algorithm 1 (Progressive encoder).
 *
 * @param[in,out] g Halfedge graph to encode.
 * @param[in,out] pm Pointmap associated with g (vertex positions).
 * @param[in,out] pq The Uniform_quantization object to apply a XYZ
 *                uniform quantization to points in pm.
 * @param[out] v_cm Empty vertex color map associated with g, used to debug
 *                topological issues during the encoding.
 * @param[out] e_cm Empty edge color map associated with g, used to debug
 *                topological issues during the encoding.
 * @param[in] gt The geometry trait object. Needed to encode the base mesh.
 * @param[in,out] EM The metric type, mutable because it calls non const
 *                methods.
 * @param[in] KP The geometric position type.
 * @param[in,out] predict The geometric predictor, mutable because it calls 
 *                non const methods.
 * @param[out] buffer Output encoded draco buffer. 
 * @param[in] measure_path The measure path to save information such as
 *            distortion.
 * @param[in] nb_max_batches The maximum number of simplification batch to
 *            apply.
 * @param[in] nb_min_vertices The minimum number of vertices to reach in the
 *            base mesh. The obtained coarse mesh will have a number of
 *            vertices less or equal to nb_min_vertices, because the condition
 *            is checked after a simplification batch has been applied.
 *            Note that the nb_max_batches may stop the decimation before
 *            the condition on nb_min_vertices is satisfied. Therefore, to
 *            make sure the the condition on nb_min_vertices is satisfied, one
 *            should set a sufficiently high nb_max_batches.
 * @param[in] batch_condition The batch stopping criterion for its edge
 *            collapses.
 * @param[in] preprocess A boolean set to true to apply preprocesses.
 * @param[in] dequantiz A boolean set to true to dequantize coordinates just
 *            after writting the compressed mesh into a file.
 * @param[in] save_preprocess A boolean set to true to save the preprocessed
 *            mesh (usefull for unit tests). Note that when preprocess is false
 *            the preprocessed mesh and the input mesh are the same.
 * @param[in] output_file_path_save_preprocess A path to save the preprocessed
 *            mesh.
 */
template< typename HalfedgeGraph,
  typename PointMap,
  typename VertexColorMap,
  typename EdgeColorMap,
  typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
  void
  progressive_compression_filter(HalfedgeGraph &g, /// Mesh to encode
                                 PointMap &pm, /// Vertex positions property map
                                 FEVV::Filters::Uniform_quantization< HalfedgeGraph, PointMap > &pq,
                                 VertexColorMap &v_cm, /// Vertex colors property map (debug)
                                 EdgeColorMap &e_cm, /// Edge colors property map (debug)
                                 const GeometryTraits &gt, /// Geometry trait object
                                 FEVV::Filters::Error_metric< HalfedgeGraph, 
                                                              PointMap >* EM, /// Local metric type
                                 const FEVV::Filters::Kept_position< HalfedgeGraph,
                                                                     PointMap>* KP, /// Geometric position type
                                 FEVV::Filters::Predictor< HalfedgeGraph, 
                                                           PointMap >* predict, /// Geometric predictor
                                 int nb_q_bits, /// Number of quantization bits
                                 int nb_max_batches,
                                 int nb_min_vertices,
                                 FEVV::Filters::BATCH_CONDITION batch_condition, /// batch stopping criterion for its edge collapses
                                 draco::EncoderBuffer &buffer, /// Output encoded buffer
                                 const std::string &measure_path,
                                 bool preprocess = true, /// true to remove isolated vertices, isolated edges
                                                  /// and vertex duplicates (several compression steps
                                                  /// cannot handle these configurations)
                                 bool dequantiz = false, /// true to dequantize coordinates just after writting the compressed mesh into a file
                                                 /// mainly for display purpose
                                 bool save_preprocess = false,
                                 const std::string& output_file_path_save_preprocess = "",
                                 bool allow_duplicates=false)
{
  if (!CGAL::is_triangle_mesh(g))
    throw std::runtime_error("progressive_compression_filter cannot handle non-pure triangular mesh, exiting (nothing is done).");
  /////////////////////////////////////////////////////////////////////////////
  // Keep the original mesh (to measure distorsion with the right scale).
  size_t num_vertices = FEVV::size_of_vertices(g);
  FEVV::Filters::Geometric_metrics< HalfedgeGraph, PointMap > g_metric(g, pm);
  /////////////////////////////////////////////////////////////////////////////
  ////////////////
  // PREPROCESS //
  ////////////////

  // Getting compression parameters (informations for the bounding box)
  std::vector< double > bb = pq.get_bb_dimension();
  std::vector< double > init = pq.get_init_coord();

  // Create a header object 
  FEVV::Header_handler HH(bb,
    init,
    KP->get_type(),
    predict->get_type(),
    nb_q_bits);

  // For testing, we choose if we use preprocess or not, and we test if we save
  // the result (so we can compare it to the decompressed mesh).
  if (preprocess)
  {
    preprocess_mesh(g, pm, pq, allow_duplicates); // Apply a preprocess to remove isolated vertices and edges
                                                  // then quantize vertex positions,
                                                  // then move duplicate positions when allow_duplicates is false.
  }
  if (save_preprocess)
  {
    FEVV::PMapsContainer pmaps_bag_empty;

    put_property_map(FEVV::edge_color, g, pmaps_bag_empty, e_cm);
    put_property_map(FEVV::vertex_color, g, pmaps_bag_empty, v_cm);
    FEVV::Filters::write_mesh(output_file_path_save_preprocess, g, pmaps_bag_empty);
  }
  double length = pq.get_diagonal();
  //////////////////////////
  // PREPROCESS DONE HERE //
  //////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  bool measure = false; // Generate csv data file for encoded mesh.

  //////////////////////////////////////////////////
  // GENERATE ENCODED INFO FOR THE DIFFERENT LODs //
  //////////////////////////////////////////////////

  // Object used to collapse edges. Also capable of computing L2 distance 
  // between 2 meshes (distortion).
  FEVV::Filters::Batch_collapser<
    HalfedgeGraph,
    PointMap,
    FEVV::Filters::Error_metric< HalfedgeGraph,
    PointMap >,
    EdgeColorMap,
    VertexColorMap >
    batch(g, pm, EM, predict, v_cm, e_cm, batch_condition);
  auto nb_vertices_current_last = FEVV::size_of_vertices(g);

  // Implements the first for loop of Algorithm 1.
  for (int i = 0; i < nb_max_batches; i++)
  {
    batch.collapse_batch(); // Simplify the mesh (fine to coarse)
                            // and encode refinement information. 
                            // Correspond to the while loop of
                            // Algorithm 1.
    if (measure)
    {
      bool skip = ((i % 5) != 0);
      batch.compute_distortion_l2(g_metric, HH, length, skip);
    }
    auto nb_vertices_current = FEVV::size_of_vertices(g);

    if ((nb_vertices_current_last == nb_vertices_current) ||
        (static_cast<int>(nb_vertices_current) <= nb_min_vertices))
    { // Last batch was not capable of removing any vertex or
      // the min number of vertices is reached.
      break;
    }

    nb_vertices_current_last = nb_vertices_current;
  }
  //////////////////////////////////////////////
  // ENCODED INFORMATION GENERATION DONE HERE //
  //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////
  // Add encoded info into the buffer  //
  ///////////////////////////////////////
  auto refinements = batch.get_refinements(); // get refinements 

  FEVV::Filters::
    Binary_batch_encoder<HalfedgeGraph, PointMap>
    binaryencode; // binary file encoding (lightweight)	

  size_t header_size = 0;
  if (measure)
  {
    binaryencode.init_measure_file(measure_path);
  }
  /////////////////////////
  // Add parameters info //
  /////////////////////////
  header_size += HH.encode_binary_header(buffer); // Implements line 35 of Algorithm 1.

  //////////////////////////////////////////////
  // Add encoded base mesh with its positions //
  //////////////////////////////////////////////
  size_t size = binaryencode.quantize_encode_coarse_mesh(g, pm, gt, buffer); // Implements line 36 of Algorithm 1.

  std::cout << "size coarse mesh: " << size << std::endl;
  header_size += size;
  /////////////////////////////////////////////
  // Add refinement info from coarse to fine //
  /////////////////////////////////////////////
  int num_batch = 0;
  auto dist_RMSE = batch.get_RMSE_distortion();
  auto dist_hausdorff = batch.get_hausdorff_distortion();
  int cumulative_bitmask = 0;
  int cumulative_residuals = 0;
  int cumulative_connectivity = 0;
  int cumulative_otherinfo = 0;
  // Implements lines 37 to 43 of Algorithm 1.
  for (int i = static_cast<int>(refinements.size() - 1); i >= 0;
    --i) // refinement batches have to be written in reverse order: from
         // coarse to fine (during the simplification step we go from fine to
         // coarse)
  {
    if (refinements[i].get_bitmask().size() > 0 &&
      refinements[i].get_connectivity().size() > 0 &&
      refinements[i].get_error_prediction().size() > 0)
    {
      // encode all of our data with draco
      auto data_bitmask =
        binaryencode.encode_bitmask(refinements[i].get_bitmask(), buffer);
      auto data_connectivity =
        binaryencode.encode_bitmask(refinements[i].get_connectivity(), buffer);

      std::pair< int, double > data_other_info = std::make_pair(0, 0.0);
      size_t size_other_info = 0; // reverse information

      if (KP->get_type() == FEVV::Filters::VKEPT_POSITION::HALFEDGE)
      {
        data_other_info = binaryencode.encode_bitmask(
          refinements[i].get_reverse_bool(), buffer);
        size_other_info = refinements[i].get_reverse_bool().size();
      }

      auto data_residuals = binaryencode.encode_residuals(
        refinements[i].get_error_prediction(), buffer, HH.get_quantization());

      cumulative_bitmask += (data_bitmask.first);
      cumulative_residuals += (data_residuals.first);
      cumulative_connectivity += (data_connectivity.first);
      cumulative_otherinfo += (data_other_info.first);

      if (measure)
      {
        binaryencode.add_line_measure_file(
          measure_path,
          num_batch,
          data_bitmask, // in bits
          data_connectivity, // in bits
          data_residuals, // in bits
          data_other_info, // in bits
          dist_RMSE[i],
          dist_hausdorff[i],
          refinements[i].get_bitmask().size(), // in bits
          refinements[i].get_connectivity().size(), // in bits
          size_other_info, // in bits
          header_size, // in bits
          num_vertices);
      }
      num_batch++;
    }
#ifdef _DEBUG
    else
    {
      std::cout << "nothing to encode" << std::endl;
    }
#endif
  }
}

/**
 * \brief Takes a mesh g, applies batches of simplification until either the 
 * number of max batches or the minimum number of vertices is reached. After 
 * that, it encodes the coarse mesh and the refinement data to a binary file.
 *
 * @param[in,out] g Halfedge graph to encode.
 * @param[in,out] pm Pointmap associated with g (vertex positions).
 * @param[out] v_cm Empty vertex color map associated with g, used to debug
 *                topological issues during the encoding.
 * @param[out] e_cm Empty edge color map associated with g, used to debug
 *                topological issues during the encoding.
 * @param[in] gt The geometry trait object. Needed to encode the base mesh.
 * @param[in] params The progressive compression parameters (metric type,
 *            position type, prediction type, nb_q_bits, allow_duplicates).
 * @param[in] nb_max_batches The maximum number of simplification batch to
 *            apply.
 * @param[in] nb_min_vertices The minimum number of vertices to reach in the
 *            base mesh. The obtained coarse mesh will have a number of
 *            vertices less or equal to nb_min_vertices, because the condition
 *            is checked after a simplification batch has been applied.
 *            Note that the nb_max_batches may stop the decimation before
 *            the condition on nb_min_vertices is satisfied. Therefore, to
 *            make sure the the condition on nb_min_vertices is satisfied, one
 *            should set a sufficiently high nb_max_batches.
 * @param[in] batch_condition The batch stopping criterion for its edge
 *            collapses.
 * @param[in] output_path The main path to set the measure_path and the
 *            binary_path (only when set to "" for the latter).
 * @param[in,out] binary_path The path to write the compressed mesh file.
 *                If not set by the user (=""), then binary_path is set to
 *                output_path+predictor+metric+keptposition+quantization.bin
 * @param[in] preprocess A boolean set to true to apply preprocesses. 
 * @param[in] dequantiz A boolean set to true to dequantize coordinates just
 *            after writting the compressed mesh into a file.
 * @param[in] save_preprocess A boolean set to true to save the preprocessed 
 *            mesh (usefull for unit tests). Note that when preprocess is false
 *            the preprocessed mesh and the input mesh are the same. 
 * @param[in] output_file_path_save_preprocess A path to save the preprocessed
 *            mesh.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename EdgeColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
progressive_compression_filter(HalfedgeGraph &g, /// Mesh to encode
                               PointMap &pm, /// Vertex positions property map
                               VertexColorMap &v_cm, /// Vertex colors property map (debug)
                               EdgeColorMap &e_cm, /// Edge colors property map (debug)
                               const GeometryTraits &gt, /// Geometry trait object
                               const FEVV::Filters::Parameters &params, /// Progressive compression parameters 
                               int nb_max_batches,
                               int nb_min_vertices,
                               FEVV::Filters::BATCH_CONDITION batch_condition, /// batch stopping criterion for its edge collapses
                               const std::string &output_path, /// needed by paths automatically set (e.g. measure_path,
                                                               /// binary_path when set to "", etc.)
                               std::string &binary_path, /// If the binary path is not set by the user, 
                                                         /// it will automatically set the path to 
                                                         /// output_path+predictor+metric+keptposition+quantization.bin
                               bool preprocess = true, /// true to remove isolated vertices, isolated edges
                                                /// and vertex duplicates (several compression steps
                                                /// cannot handle these configurations)
                               bool dequantiz = false, /// true to dequantize coordinates just after writting the compressed mesh into a file
                                               /// mainly for display purpose
                               bool save_preprocess = false,
                               const std::string& output_file_path_save_preprocess = "")
{
#ifdef TIMING
  auto time_point_before_comp = std::chrono::high_resolution_clock::now();
#endif
  /////////////////////////////////////////////////////////////////////////////
    // Create quantization object.
  FEVV::Filters::Uniform_quantization< HalfedgeGraph, PointMap > pq(
    g, pm, params.get_quantization());

  // Create dequantization object (will be used in metrics).
  // Its construction is partialy based on the quantization object.
  FEVV::Filters::Uniform_dequantization< HalfedgeGraph, PointMap > dq(
    g,
    pm,
    params.get_quantization(),
    pq.get_bb_dimension(),
    pq.get_init_coord());

    // Create geometric position type object.
  FEVV::Filters::Kept_position< HalfedgeGraph,
    PointMap> *KP = nullptr;

  if (params.get_vkept_position() ==
    FEVV::Filters::VKEPT_POSITION::HALFEDGE) // Kept at the position of the
                                             // halfedge target
  {
#ifdef _DEBUG
    std::cout << "HALFEDGE" << std::endl;
#endif
    KP = new FEVV::Filters::
      Halfedge< HalfedgeGraph, PointMap >(
        g, pm);
  }
  else if (params.get_vkept_position() ==
    FEVV::Filters::VKEPT_POSITION::MIDPOINT) // kept at the middle of the
                                             // halfedge
  {
#ifdef _DEBUG
    std::cout << "MIDPOINT" << std::endl;
#endif
    KP = new FEVV::Filters::
      Midpoint< HalfedgeGraph, PointMap >(
        g, pm);
  }
  else
    throw std::runtime_error("progressive_compression_filter cannot handle kept position type.");

  // Create metric object.
  FEVV::Filters::
    Error_metric< HalfedgeGraph, PointMap >
    *EM = nullptr;

  if (params.get_metric() == FEVV::Filters::METRIC_TYPE::EDGE_LENGTH)
  {
    EM = new FEVV::Filters::Edge_length_metric< HalfedgeGraph,
      PointMap >(
        g, pm, KP, dq);
#ifdef _DEBUG
    std::cout << "EDGE_LENGTH" << std::endl;
#endif
  }
  else if (params.get_metric() == FEVV::Filters::METRIC_TYPE::VOLUME_PRESERVING)
  {
    EM = new FEVV::Filters::Volume_preserving<
      HalfedgeGraph,
      PointMap >( // metric computing the local volume difference between
                      // the locally collapsed neighborhood and the original one
        g,
        pm,
        KP,
        dq);
#ifdef _DEBUG
    std::cout << "VOLUME_PRESERVING" << std::endl;
#endif
  }
  else if (params.get_metric() == FEVV::Filters::METRIC_TYPE::QEM_3D)
  {
    EM = new FEVV::Filters::QEM_3D< HalfedgeGraph,
      PointMap >( // Quadric error metric
        g,
        pm,
        KP,
        dq);
#ifdef _DEBUG
    std::cout << "QEM_3D" << std::endl;
#endif
  }
  else
  {
    delete KP;
    throw std::runtime_error("progressive_compression_filter cannot handle metric type.");
  }

  // Create geometric predictor object.
  FEVV::Filters::Predictor< HalfedgeGraph,
    PointMap > *predict = nullptr;

  if (params.get_prediction() == FEVV::Filters::PREDICTION_TYPE::BUTTERFLY)
  {
    predict = new FEVV::Filters::
      Butterfly< HalfedgeGraph, PointMap >(g, KP, pm);
#ifdef _DEBUG
    std::cout << "BUTTERFLY" << std::endl;
#endif
  }
  else if (params.get_prediction() == FEVV::Filters::PREDICTION_TYPE::DELTA)
  {
    predict = new FEVV::Filters::Delta_predictor< HalfedgeGraph,
      PointMap >(g, KP, pm);
#ifdef _DEBUG
    std::cout << "DELTA" << std::endl;
#endif
  }
  else if (params.get_prediction() == FEVV::Filters::PREDICTION_TYPE::POSITION)
  {
    predict = new FEVV::Filters::Raw_positions< HalfedgeGraph,
      PointMap >(g, KP, pm);
#ifdef _DEBUG
    std::cout << "POSITION" << std::endl;
#endif
  }
  else
  {
    delete KP;
    delete EM;
    throw std::runtime_error("progressive_compression_filter cannot handle geometric prediction type.");
  }

  std::string measure_path = // used when measure is set to true (Generate csv 
                             // data file for encoded mesh)
    output_path + predict->get_as_string() +
    EM->get_as_string() + KP->get_as_string() +
    std::to_string(params.get_quantization()) + ".csv";
  // If the binary path is not set by the user, it will automatically set the
  // path to output_path[predictor][metric][keptposition][quantization].bin
  if (binary_path == "")
  {
    binary_path = output_path + predict->get_as_string() +
      EM->get_as_string() + KP->get_as_string() +
      std::to_string(params.get_quantization()) + ".bin";
  }
  /////////////////////////////////////////////////////////////////////////////
  draco::EncoderBuffer buffer;

  progressive_compression_filter(g,
                                 pm, 
                                 pq,
                                 v_cm,
                                 e_cm, 
                                 gt, 
                                 EM,
                                 KP,
                                 predict,
                                 params.get_quantization(),
                                 nb_max_batches,
                                 nb_min_vertices,
                                 batch_condition,
                                 buffer,
                                 measure_path,
                                 preprocess, 
                                 dequantiz,
                                 save_preprocess,
                                 output_file_path_save_preprocess,
                                 params.get_allow_duplicates());
  /////////////////////////////////////////////////////////////////////////////
  // write this info into binary file
  std::ofstream binary_output_file;

  binary_output_file.open(binary_path, std::fstream::binary | std::ofstream::out | std::ofstream::trunc);
  binary_output_file.write(buffer.data(), buffer.size());
  binary_output_file.close();

  delete KP;
  delete EM;
  delete predict;

  if (dequantiz)
  {
    dq.point_dequantization();
  }

#ifdef TIMING
  auto time_point_after_comp = std::chrono::high_resolution_clock::now();
  auto duration_comp = std::chrono::duration_cast<std::chrono::milliseconds>(time_point_after_comp - time_point_before_comp);
  std::cout << "Compression time: " << duration_comp.count() << " milliseconds" << std::endl;
#endif
}

/**
 * \brief Takes a mesh g, applies batches of simplification until either the
 * number of max batches or the minimum number of vertices is reached. After
 * that, it encodes the coarse mesh and the refinement data to a binary file.
 * The geometry trait object is set automatically (syntactic sugar).
 *
 * @param[in,out] g Halfedge graph to encode.
 * @param[in,out] pm Pointmap associated with g (vertex positions).
 * @param[out] v_cm Empty vertex color map associated with g, used to debug
 *                topological issues during the encoding.
 * @param[out] e_cm Empty edge color map associated with g, used to debug
 *                topological issues during the encoding.
 * @param[in] params The progressive compression parameters.
 * @param[in] output_path The main path to set the measure_path and the
 *            binary_path (only when set to "" for the latter).
 * @param[in,out] binary_path The path to write the compressed mesh file.
 *                If not set by the user (=""), then binary_path is set to
 *                output_path+predictor+metric+keptposition+quantization.bin
 * @param[in] nb_max_batches The maximum number of simplification batch to
 *            apply.
 * @param[in] nb_min_vertices The minimum number of vertices to reach in the
 *            base mesh. The obtained coarse mesh will have a number of
 *            vertices less or equal to nb_min_vertices, because the condition
 *            is checked after a simplification batch has been applied.
 *            Note that the nb_max_batches may stop the decimation before
 *            the condition on nb_min_vertices is satisfied. Therefore, to
 *            make sure the the condition on nb_min_vertices is satisfied, one
 *            should set a sufficiently high nb_max_batches.
 * @param[in] batch_condition The batch stopping criterion for its edge
 *            collapses.
 * @param[in] preprocess A boolean set to true to apply preprocesses.
 * @param[in] dequantiz A boolean set to true to dequantize coordinates just
 *            after writting the compressed mesh into a file.
 * @param[in] save_preprocess A boolean set to true to save the preprocessed
 *            mesh (usefull for unit tests). Note that when preprocess is false
 *            the preprocessed mesh and the input mesh are the same.
 * @param[in] output_file_path_save_preprocess A path to save the preprocessed
 *            mesh.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename EdgeColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
progressive_compression_filter(HalfedgeGraph &g, /// input mesh
                               PointMap &pm, /// its vertex position map
                               VertexColorMap &v_cm, /// its vertex color map (debug)
                               EdgeColorMap &e_cm, /// its edge color map (debug)
                               const FEVV::Filters::Parameters &params, /// input progressive compression param
                               const std::string &output_path,
                               std::string &binary_path,
                               int nb_max_batches,
                               int nb_min_vertices,
                               FEVV::Filters::BATCH_CONDITION batch_condition, /// batch stopping criterion for its edge collapses
                               bool preprocess = true,
                               bool dequantiz = false,
                               bool save_preprocess = false,
                               const std::string& output_file_path_save_preprocess = "")

{
  GeometryTraits gt(g);
  progressive_compression_filter(g,
                                 pm,
                                 v_cm,
                                 e_cm,
                                 gt,
                                 params,
                                 output_path,
                                 binary_path,
                                 nb_max_batches,
                                 nb_min_vertices,
                                 batch_condition,
                                 preprocess,
                                 dequantiz,
                                 save_preprocess, 
                                 output_file_path_save_preprocess);
}

} // namespace Filters
} // namespace FEVV