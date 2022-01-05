#pragma once

#include "FEVV/Wrappings/Geometry_traits.h"
#include <boost/graph/graph_traits.hpp>
#include "FEVV/Wrappings/properties.h"
#include <CGAL/boost/graph/helpers.h> // is_triangle_mesh
#include "FEVV/Filters/CGAL/Progressive_Compression/Parameters.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/EdgeLengthMetric.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/VolumePreserving.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/QEM3D.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/BatchCollapser.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/RawPositions.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/DeltaPredictor.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/HeaderHandler.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/BinaryBatchEncoder.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/geometric_metrics.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/quantification.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/dequantization.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/preprocessingLD.h"
#include "FEVV/Filters/Generic/minmax_map.h"
#include "FEVV/Filters/Generic/color_mesh.h"

template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap >
void
preprocess_mesh(
    HalfedgeGraph &g,
    PointMap &pm,
    VertexColorMap &v_cm,
    FEVV::Filters::UniformQuantization< HalfedgeGraph, PointMap > &pq)
{

  FEVV::Filters::PreprocessingLD< HalfedgeGraph, PointMap, VertexColorMap >
      preprocess(g, pm, v_cm);
  preprocess.process_mesh_before_quantization();

  // vertex positions part
  pq.find_min_and_max();
  pq.set_bounding_box();
  pq.set_quantization_step();
  pq.point_quantization();
 
  preprocess.process_mesh_after_quantization();
}
/**
 * \brief Simplifies a mesh and encodes the necessary information to restore it to its
 * original state.
 *
 * For a given mesh g, applies batches of simplification until either the 
 * number of max batches or the minimum number of vertices is reached.  After 
 * that, encodes the coarse mesh and the refinement data to a binary file.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename EdgeColorMap,
          //typename VertexNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
progressive_compression_filter(HalfedgeGraph &g, /// Mesh to encode
                               PointMap &pm, /// Vertex positions property map
                               VertexColorMap &v_cm, /// Vertex colors property map
                               EdgeColorMap &e_cm, /// Edge colors property map
                               //VertexNormalMap &/*v_nm*/, /// Vertex normals property map
                               const GeometryTraits &gt, /// Geometry trait object
                               const FEVV::Filters::Parameters &params, /// Progressive compression parameters 
                               const std::string &output_path, /// needed by paths automatically set (e.g. measure_path,
                                                               /// binary_path when set to "", etc.)
                               std::string &binary_path, /// If the binary path is not set by the user, 
                                                         /// it will automatically set the path to 
                                                         /// output_path+predictor+metric+keptposition+quantization.bin
                               int nb_max_batches,
                               int nb_min_vertices,
                               FEVV::Filters::BATCH_CONDITION batch_condition, /// batch stopping criterion for its edge collapses
                               bool preprocess, /// true to remove isolated vertices, isolated edges
                                                /// and vertex duplicates (several compression steps
                                                /// cannot handle these configurations)
                               bool dequantif,
                               bool save_preprocess,
                               const std::string& output_file_path_save_preprocess)
{
  if (!CGAL::is_triangle_mesh(g))
    throw std::runtime_error("progressive_compression_filter cannot handle non-pure triangular mesh, exiting (nothing is done).");

  FEVV::Filters::KeptPosition< HalfedgeGraph,
                               PointMap> *KP = nullptr;

  if(params.getVKeptPosition() ==
     FEVV::Filters::VKEPT_POSITION::HALFEDGE) // Kept at the position of the
                                              // halfedge target
  {
#ifdef _DEBUG
    std::cout << "Halfedge" << std::endl;
#endif
    KP = new FEVV::Filters::
        Halfedge< HalfedgeGraph, PointMap >(
            g, pm);
  }
  else if(params.getVKeptPosition() ==
     FEVV::Filters::VKEPT_POSITION::MIDPOINT) // kept at the middle of the
                                              // halfedge
  {
#ifdef _DEBUG
    std::cout << "Midpoint" << std::endl;
#endif
    KP = new FEVV::Filters::
        Midpoint< HalfedgeGraph, PointMap >(
            g, pm);
  }
  else
    throw std::runtime_error("progressive_compression_filter cannot handle kept position type.");

  // keeping the original mesh (to measure distorsion with the right scale)

  auto vertices_pair = vertices(g);
  size_t num_vertices =
      std::distance(vertices_pair.first, vertices_pair.second);
  FEVV::Filters::GeometricMetrics< HalfedgeGraph, PointMap > g_metric(g, pm);
  g_metric.initialize_AABB_tree_for_init_LoD();
  g_metric.subsample_LoD_init();

  // PREPROCESS

  // Create quantization objects
  FEVV::Filters::UniformQuantization< HalfedgeGraph, PointMap > pq(
      g, pm, params.getQuantization());

  // for testing, we choose if we use preprocess or not, and we test if we save
  // the result (so we can compare it to the decompressed mesh)
  if(preprocess)
  {
    preprocess_mesh(g, pm, v_cm, pq); // quantize vertex positions and 
                                                 // UV coordinates
  }
  if(save_preprocess)
  {
    FEVV::PMapsContainer pmaps_bag_empty;

    put_property_map(FEVV::edge_color, g, pmaps_bag_empty, e_cm);
    put_property_map(FEVV::vertex_color, g, pmaps_bag_empty, v_cm);
    FEVV::Filters::write_mesh(output_file_path_save_preprocess, g, pmaps_bag_empty);
  }

  double length = pq.getDiagonal();


  // create dequantization objects (will be used in metrics)
  FEVV::Filters::UniformDequantization< HalfedgeGraph, PointMap > dq(
      g,
      pm,
      params.getQuantization(),
      pq.get_bb_dimension(),
      pq.get_init_coord());

  dq.set_max_length();
  dq.set_quantization_step();

  // PREPROCESS DONE HERE

  bool measure = false; // Generate csv data file for encoded mesh

  // Getting compression parameters (informations for the bounding box)


  std::vector< double > bb = pq.get_bb_dimension();
  std::vector< double > init = pq.get_init_coord();


  // Create a header object (will be written in a file)
  FEVV::HeaderHandler HH(bb,
                         init,
                         params.getVKeptPosition(),
                         params.getPrediction(),
                         params.getQuantization() );

  // Object used to compute L2 distance between 2 meshes (Distortion)


  // Create abstract Metric (in order to compute the weight of each edge) and
  // KeptPosition (Where the vertex resulting from a collapse is placed)
  // objects: they will be "overriden" by an object from a derived class
  FEVV::Filters::
      ErrorMetric< HalfedgeGraph, PointMap >
          *EM = nullptr;


  if(params.getMetric() == FEVV::Filters::METRIC_TYPE::EDGE_LENGTH)
  {
    EM = new FEVV::Filters::EdgeLengthMetric< HalfedgeGraph,
                                              PointMap >(
        g, pm, KP, dq);
#ifdef _DEBUG
    std::cout << "Edge Length" << std::endl;
#endif
  }
  else if(params.getMetric() == FEVV::Filters::METRIC_TYPE::VOLUME_PRESERVING)
  {
    EM = new FEVV::Filters::VolumePreserving<
        HalfedgeGraph,
        PointMap >( // metric computing the local volume difference between
                        // the locally collapsed neighborhood and the original one
        g,
        pm,
        KP,
        dq);
#ifdef _DEBUG
    std::cout << "Volume preserving" << std::endl;
#endif
  }
  else if(params.getMetric() == FEVV::Filters::METRIC_TYPE::QEM)
  {
    EM = new FEVV::Filters::QEM3D< HalfedgeGraph,
                                   PointMap >( // Quadric error metric
        g,
        pm,
        KP,
        dq);
#ifdef _DEBUG
    std::cout << "QEM" << std::endl;
#endif
  }
  else
  {
    delete KP;
    throw std::runtime_error("progressive_compression_filter cannot handle metric type.");
  }

  // create geometric predictor object
  FEVV::Filters::Predictor< HalfedgeGraph,
    PointMap > *predict = nullptr;

  if(params.getPrediction() == FEVV::Filters::PREDICTION_TYPE::BUTTERFLY)
  {
    predict = new FEVV::Filters::
        Butterfly< HalfedgeGraph, PointMap >(
            g, KP, pm);
  }
  else if(params.getPrediction() == FEVV::Filters::PREDICTION_TYPE::DELTA)
  {
    predict = new FEVV::Filters::DeltaPredictor< HalfedgeGraph,
                                                 PointMap >(
        g, KP, pm);
  }
  else
  {
    delete KP;
    delete EM;
    throw std::runtime_error("progressive_compression_filter cannot handle geometric prediction type.");
  }

  typename FEVV::Vertex_pmap_traits< HalfedgeGraph, int >::pmap_type index =
      FEVV::make_vertex_property_map< HalfedgeGraph, int >(
          g); // will be useless, only initiliazed because the older version of
              // the code uses it (will eventually be removed)

  // object used to collapse edges
  FEVV::Filters::BatchCollapser<
      HalfedgeGraph,
      PointMap,
      FEVV::Filters::ErrorMetric< HalfedgeGraph,
                                  PointMap >,
      typename FEVV::Vertex_pmap_traits< HalfedgeGraph, int >::pmap_type,
      EdgeColorMap,
      VertexColorMap >
      batch(g, pm, EM, predict, index, v_cm, e_cm, batch_condition);
  auto current_vertices_pair = vertices(g);
  auto nb_vertices_current_last = std::distance(current_vertices_pair.first, current_vertices_pair.second);
  for(int i = 0; i < nb_max_batches; i++)
  {
    batch.collapse_batch(); // simplify mesh (fine to coarse)
    if(measure)
    {
      bool first_call = (i == 0) ? true : false;
      bool skip = ((i % 5) != 0);
      batch.compute_distorsion_l2(g_metric, HH, length, first_call, skip);
    }
    current_vertices_pair = vertices(g);
    auto nb_vertices_current = std::distance(current_vertices_pair.first,
                                             current_vertices_pair.second);
    //std::cout << nb_vertices_current << " " << nb_min_vertices << std::endl;

    if (nb_vertices_current_last == nb_vertices_current)
    { // last batch was not capable of removing any vertex
      break;
    }

    nb_vertices_current_last = nb_vertices_current;

    if(nb_vertices_current < nb_min_vertices)
    {
      break;
    }
  }
  auto refinements = batch.getRefinements(); // get refinements in reverse order
                                             // (coarse to fine)
  FEVV::Filters::
      BinaryBatchEncoder< HalfedgeGraph, PointMap>
          binaryencode; // binary file encoding (lightweight)	


  std::string measure_path = // used when measure is set to true (Generate csv 
                             // data file for encoded mesh)
                             output_path + predict->getMethodasString() +
                             EM->getMethodasString() + KP->getMethodasString() +
                             std::to_string(HH.getQuantization()) + ".csv";
  // If the binary path is not set by the user, it will automatically set the
  // path to output_path[predictor][metric][keptposition][quantization].bin
  if(binary_path == "")
  {
    binary_path = output_path + predict->getMethodasString() +
                 EM->getMethodasString() + KP->getMethodasString() +
                 std::to_string(HH.getQuantization()) + ".bin";
  }

  size_t header_size = 0;
  if(measure)
  {
    binaryencode.InitMeasureFile(measure_path);
  }

  // write coarse mesh
  draco::EncoderBuffer
      buffercoarsemesh; // This object will contain the encoded coarse mesh
                        // (encoded with draco single-rate compression)

  size_t size = 0; 
  // encode a mesh with positions
  size = binaryencode.QuantizeEncodeCoarseMesh(g, pm, gt, buffercoarsemesh); 

  std::cout << "size coarse mesh: " << size << std::endl;
  header_size += size;

  // write this info into binary file
  //std::ofstream coarse_mesh;
  //coarse_mesh.open(coarse_path, std::fstream::binary);
  //coarse_mesh.write(buffercoarsemesh.data(), buffercoarsemesh.size());
  //coarse_mesh.close();

  std::ofstream binary_output_file;
  binary_output_file.open(binary_path, std::fstream::binary);
  binary_output_file.write(buffercoarsemesh.data(), buffercoarsemesh.size());
  binary_output_file.close();

  // write header info into binary file (opened separately in order to empty the
  // file before filling it up again)
  header_size += HH.EncodeBinaryHeader(binary_path);

  binary_output_file.open(binary_path, std::fstream::binary | std::fstream::app);

  // write refinement info into both binary and text files (binary for real use,
  // text for debug)
  int num_batch = 0;
  auto dist = batch.getDistorsion();
  int cumulative_bitmask = 0;
  int cumulative_residuals = 0;
  int cumulative_connectivity = 0;
  int cumulative_otherinfo = 0;
  for(int i = static_cast< int >(refinements.size() - 1); i >= 0;
      --i) // refinement batches have to be written in reverse order: from
           // coarse to fine (during the simplification step we go from fine to
           // coarse)
  {
    if(refinements[i].get_bitmask().size() > 0 &&
       refinements[i].get_connectivity().size() > 0 &&
       refinements[i].get_error_prediction().size() > 0)
    {

      // encode all of our data wth draco
      draco::EncoderBuffer buffer;
      auto data_bitmask =
          binaryencode.EncodeBitMask(refinements[i].get_bitmask(), buffer);
      auto data_connectivity =
          binaryencode.EncodeBitMask(refinements[i].get_connectivity(), buffer);
      auto data_residuals = binaryencode.EncodeResiduals(
          refinements[i].get_error_prediction(), buffer, HH.getQuantization());

      std::pair< int, double > data_other_info = std::make_pair(0, 0.0);
      size_t size_other_info = 0; // reverse information

      if(KP->getType() == FEVV::Filters::VKEPT_POSITION::HALFEDGE)
      {
        data_other_info = binaryencode.EncodeBitMask(
            refinements[i].get_reverse_bool(), buffer);
        size_other_info = refinements[i].get_reverse_bool().size();
      }


      cumulative_bitmask += (data_bitmask.first) * 8;
      cumulative_residuals += (data_residuals.first) * 8;
      cumulative_connectivity += (data_connectivity.first) * 8;
      cumulative_otherinfo += (data_other_info.first) * 8;	

      for(auto buf_it = buffer.buffer()->begin();
          buf_it != buffer.buffer()->end();
          buf_it++)
      {

        binary_output_file << *buf_it;
      }

      if(measure)
      {
        binaryencode.AddLineMeasureFile(
            measure_path,
            num_batch,
            data_bitmask,
            data_connectivity,
            data_residuals,
            data_other_info,
            dist[i],
            refinements[i].get_bitmask().size(),
            refinements[i].get_connectivity().size(),
            size_other_info,
            header_size,
            num_vertices);
      }
      num_batch++;
    }
    else
    {
      std::cout << "nothing to encode" << std::endl;
    }
  }
  binary_output_file.close();

  if(dequantif)
  {
    dq.point_dequantization();
  }
  delete KP;
  delete EM;
  delete predict;
}

// Helper function to simplify (syntactic sugar) the call to
// progressive_compression_filter
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename EdgeColorMap,
          //typename VertexNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
progressive_compression_filter(HalfedgeGraph &g, /// input mesh
                               PointMap &pm, /// its vertex position map
                               VertexColorMap &v_cm, /// its vertex color map
                               EdgeColorMap &e_cm, /// its edge color map
                               //VertexNormalMap &v_nm, /// its vertex normal map
                               const FEVV::Filters::Parameters &params, /// input progressive compression param
                               const std::string &output_path,
                               std::string &binary_path,
                               int nb_max_batches,
                               int nb_min_vertices,
                               FEVV::Filters::BATCH_CONDITION batch_condition, /// batch stopping criterion for its edge collapses
                               bool preprocess,
                               bool dequantif,
                               bool save_preprocess,
                               const std::string& output_file_path_save_preprocess)

{
  GeometryTraits gt(g);
  progressive_compression_filter(g,
                                 pm,
                                 v_cm,
                                 e_cm,
                                 //v_nm,
                                 gt,
                                 params,
                                 output_path,
                                 binary_path,
                                 nb_max_batches,
                                 nb_min_vertices,
                                 batch_condition,
                                 preprocess,
                                 dequantif,
                                 save_preprocess, 
                                 output_file_path_save_preprocess);
}
