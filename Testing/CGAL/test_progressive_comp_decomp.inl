#pragma once
#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"

#include "FEVV/Filters/CGAL/Progressive_Compression/Helpers/Initializer_compression.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/progressive_decompression_filter.hpp"

#include "FEVV/Tools/IO/FileUtilities.hpp"

#include "Testing/Utils/utils_are_meshes_identical.hpp"

/** \brief Compresses then decompresses a mesh. 
 * \param  argc number of arguments in argv.
 * \param  argv array of C-like strings
 * \return the test status.
 *         0 when the test success (decompressed version is equal to the original mesh)
 *         -1 otherwise.
 */
template< typename MeshT >
int
test_automatic_progressive_compression_decompression(int argc,
                                                     const char **argv
                                                                       )
{
  if (argc < 2)
  {
    std::cout << "test_automatic_progressive_compression_decompression_ld needs an input mesh filename." << std::endl;
    return -1;
  }
  /////////////////////////////////////////////////////////////////////////////
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  /////////////////////////////////////////////////////////////////////////////
  // input and output files
  std::string input_file_path = argv[1];
  std::string output_file_path_save_preprocess = ""; 
  std::string output_file_path_decomp = "";
  std::string measurepath = typeid(m).name();
  if (measurepath.find("Surface_mesh") != std::string::npos)
  {
    measurepath = "Surface_mesh_";
  }
  else if (measurepath.find("Polyhedron_3") != std::string::npos)
  {
    measurepath = "Polyhedron_3_";
  }
  else if (measurepath.find("Linear_cell_complex") != std::string::npos)
  {
    measurepath = "LCC_";
  }
  else if (measurepath.find("AIF") != std::string::npos)
  {
    measurepath = "AIF_";
  }
  else
  {
    measurepath = "Other_mesh_type_";
  }
  std::string compressed_mesh_binary_file = measurepath + FEVV::FileUtils::get_file_name(input_file_path);
  output_file_path_save_preprocess = "progressive_compression_ld_" + compressed_mesh_binary_file + "_original_mesh_after_preprocess.off"; // default value
  output_file_path_decomp = "progressive_decompression_ld_" + compressed_mesh_binary_file + "_output.off"; // default value
  compressed_mesh_binary_file = compressed_mesh_binary_file + ".bin";
  /////////////////////////////////////////////////////////////////////////////
  if (argc > 2)
  {
    output_file_path_save_preprocess = argv[2];
  }
  if (argc > 3)
  {
    output_file_path_decomp = argv[3];
  }
  if (argc > 4)
  {
    std::cout << "test_automatic_progressive_compression_decompression_ld needs at maximum 3 input parameters: input mesh filename, output filename for the preprocessed mesh, and an output filename for the decompressed mesh." << std::endl;
    return -1;
  }
  /////////////////////////////////////////////////////////////////////////////
  // read mesh from file
  FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

  typename FEVV::PMap_traits< FEVV::vertex_color_t, MeshT >::pmap_type v_cm;
  typename FEVV::PMap_traits< FEVV::edge_color_t, MeshT >::pmap_type e_cm;
  typename FEVV::PMap_traits< FEVV::face_normal_t, MeshT >::pmap_type f_nm;
  typename FEVV::PMap_traits< FEVV::vertex_normal_t, MeshT >::pmap_type v_nm;
  set_mesh_and_properties(
      m, pmaps_bag, v_cm, e_cm, f_nm, v_nm);

  FEVV::Filters::Parameters params(FEVV::Filters::PREDICTION_TYPE::BUTTERFLY,
                                   FEVV::Filters::VKEPT_POSITION::MIDPOINT,
                                   FEVV::Filters::METRIC_TYPE::QEM,
                                   true,
                                   false,
                                   12);
  auto pm = get(boost::vertex_point, m);
  auto gt_ = FEVV::Geometry_traits< MeshT >(m);
  /////////////////////////////////////////////////////////////////////////////
  progressive_compression_filter(m,
                                 pm,
                                 v_cm,
                                 e_cm,
                                 // v_nm,
                                 gt_,
                                 params,
                                 measurepath,
                                 compressed_mesh_binary_file,
                                 200,
                                 0,
                                 FEVV::Filters::BATCH_CONDITION::ALL_EDGES,
                                 true,
                                 true,
                                 true,
                                 output_file_path_save_preprocess);

  MeshT m_decomp;
  FEVV::PMapsContainer pmaps_bag_decomp;
  typename FEVV::PMap_traits< FEVV::vertex_color_t, MeshT >::pmap_type
      v_cm_decomp;
  typename FEVV::PMap_traits< FEVV::edge_color_t, MeshT >::pmap_type
      e_cm_decomp;
  typename FEVV::PMap_traits< FEVV::face_normal_t, MeshT >::pmap_type
      f_nm_decomp;
  typename FEVV::PMap_traits< FEVV::vertex_normal_t, MeshT >::pmap_type
      v_nm_decomp;
  auto pm_decomp = get(boost::vertex_point, m_decomp);
  set_mesh_and_properties(m_decomp,
                          pmaps_bag_decomp,
                          v_cm_decomp,
                          e_cm_decomp,
                          f_nm_decomp,
                          v_nm_decomp);


  // apply filter

  FEVV::Filters::progressive_decompression_filter(m_decomp,
                                                  pm_decomp,
                                                  v_cm_decomp,
                                                  //v_nm_decomp,
                                                  e_cm_decomp,
                                                  compressed_mesh_binary_file,
                                                  false);
  // write mesh to file
  FEVV::PMapsContainer pmaps_bag_decomp_without_material;

  put_property_map(
      FEVV::edge_color, m_decomp, pmaps_bag_decomp_without_material, e_cm);
  put_property_map(
      FEVV::vertex_color, m_decomp, pmaps_bag_decomp_without_material, v_cm);
  FEVV::Filters::write_mesh(
      output_file_path_decomp, m_decomp, pmaps_bag_decomp_without_material);

  bool equal_or_not = are_meshes_equal(output_file_path_save_preprocess,
                                       output_file_path_decomp,
                                       false,
                                       0,
	  0,
                                       true);
  if(equal_or_not)
  {
    std::cout << "meshes are equal!" << std::endl;
    return 0;
  }
  else
  {
    std::cout << "meshes are not equal" << std::endl;
    return -1;
  }

  
}
