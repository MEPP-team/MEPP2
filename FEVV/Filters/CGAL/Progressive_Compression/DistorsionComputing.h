#pragma once

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/BatchCollapser.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/geometric_metrics.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/quantification.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/dequantization.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/preprocessingLD.h"
#include "boost/filesystem.hpp"
#include <iostream>


namespace fs = boost::filesystem;


std::vector< fs::path >
get_files(fs::path const &root, std::string const &ext)
{
  std::vector< fs::path > paths;

  if(fs::exists(root) && fs::is_directory(root))
  {
    for(auto const &entry : fs::recursive_directory_iterator(root))
    {
      if(fs::is_regular_file(entry) && entry.path().extension() == ext)
        paths.emplace_back(entry.path().filename());
    }
  }
  return paths;
}

// Computes the distorsion between an unquantized mesh g and a set of quantized
// meshes (they will be unquantized) 
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename EdgeColorMap,
          typename VertexNormalMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
ComputeDistorsions(HalfedgeGraph &g, /// original mesh
                   PointMap &pm, /// original point map
                   VertexColorMap &v_cm, /// original attribute maps
                   EdgeColorMap &/*e_cm*/,
                   VertexNormalMap &/*v_nm*/,
                   const GeometryTraits &/*gt*/,
                   fs::path const &root // folder containing .obj files which 
				                        // are used to compute distorsions 
										// between them and g
                    )
{
  //initialize AABB tree for the original mesh
  FEVV::Filters::GeometricMetrics< HalfedgeGraph, PointMap > g_metric(g, pm);
  g_metric.initialize_AABB_tree_for_init_LoD();
  g_metric.subsample_LoD_init();

  // PREPROCESS: get quantization parameters
  FEVV::Filters::PreprocessingLD< HalfedgeGraph, PointMap, VertexColorMap >
      preprocess(g, pm, v_cm);
  preprocess.process_mesh_before_quantization();
  
  FEVV::Filters::UniformQuantization< HalfedgeGraph, PointMap > pq(g, pm, 12);
  pq.find_min_and_max();
  pq.set_bounding_box();
  pq.set_quantization_step();

  std::vector< double > bb = pq.get_bb_dimension();
  std::vector< double > init = pq.get_init_coord();
  auto files = get_files(root, ".obj");
  size_t nbf = files.size();
  std::vector< double > _distorsion_per_batch(nbf, 0.);
  for(size_t i = 0; i < nbf; i++)
  {
    HalfedgeGraph current_mesh;
    FEVV::PMapsContainer current_pmaps_bag;
	
	//read quantized mesh
    FEVV::Filters::read_mesh(root.string()+files[i].string(),current_mesh,current_pmaps_bag);
    auto current_pm = get(boost::vertex_point, current_mesh);
	
	//unquantize mesh
    FEVV::Filters::UniformDequantization< HalfedgeGraph, PointMap > dq(
        current_mesh,
        current_pm,
        12,
        pq.get_bb_dimension(),
        pq.get_init_coord());
    dq.set_max_length();
    dq.set_quantization_step();
    dq.point_dequantization();

	//compute distorsion
    _distorsion_per_batch[i] = 
	    g_metric.compute_symetric_L2(current_mesh, false);
  }
  //create distorsion column in	a new .csv file
  std::ofstream dist_file;
  dist_file.open(root.string() + "dist_measures" + ".csv");
  auto it_dist = _distorsion_per_batch.begin();
  dist_file << "DISTORSION"
             << "\n";
  for(; it_dist != _distorsion_per_batch.end(); it_dist++)
  {
    dist_file << (*it_dist) << "\n";
  }

  dist_file.close();
}
