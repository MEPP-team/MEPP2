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

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"

#include "FEVV/Filters/CGAL/Progressive_Compression/progressive_compression_filter.hpp"
#include "FEVV/Filters/CGAL/Progressive_Compression/Helpers/distortion_computing.h"

/**
 * \brief Function to extract or create property map associated with a
 *        property map bag (for generic processing).
 * @param[in] m The input mesh associated with the the property map bag. 
 * @param[in,out] pmaps_bag The property map bag. It is updated whenever
 *                a property map creation occurs. 
 * @param[out] v_cm The obtained vertex color map that is also added in the 
 *             pmaps_bag when not already present.
 * @param[out] e_cm The obtained edge color map that is also added in the
 *             pmaps_bag when not already present.
 * @param[out] v_nm The obtained vertex normal map that is also added in the
 *             pmaps_bag when not already present.
 */
template< typename MeshT >
void
set_mesh_and_properties(
    const MeshT &m,
    FEVV::PMapsContainer &pmaps_bag,
    typename FEVV::PMap_traits< FEVV::vertex_color_t, MeshT >::pmap_type &v_cm,
    typename FEVV::PMap_traits< FEVV::edge_color_t, MeshT >::pmap_type &e_cm,
    typename FEVV::PMap_traits< FEVV::vertex_normal_t, MeshT >::pmap_type &v_nm)

{
  // Note: the property maps must be extracted from the
  //       property maps bag, and explicitly passed as
  //       parameters to the filter, in order to make
  //       clear what property is used by the filter.

  // Retrieve or create vertex-color property map.
  if(has_map(pmaps_bag, FEVV::vertex_color))
  {
    std::cout << "use existing vertex-color map" << std::endl;
    v_cm = get_property_map(FEVV::vertex_color, m, pmaps_bag);
  }
  else
  {
    std::cout << "create vertex-color map" << std::endl;
    v_cm = make_property_map(FEVV::vertex_color, m);
    // Store property map in property maps bag.
    put_property_map(FEVV::vertex_color, m, pmaps_bag, v_cm);
  }
  // Retrieve or create edge-color property map.
  if(has_map(pmaps_bag, FEVV::edge_color))
  {
    std::cout << "using existing edge-color map" << std::endl;
    e_cm = get_property_map(FEVV::edge_color, m, pmaps_bag);
  }
  else
  {
    std::cout << "create edge-color map" << std::endl;
    e_cm = make_property_map(FEVV::edge_color, m);
    // Store property map in property maps bag.
    put_property_map(FEVV::edge_color, m, pmaps_bag, e_cm);
  }
  // Retrieve or create vertex-normal property map.
  if(has_map(pmaps_bag, FEVV::vertex_normal))
  {
    std::cout << "use existing vertex-normal map" << std::endl;
    v_nm = get_property_map(FEVV::vertex_normal, m, pmaps_bag);
  }
  else
  {
    std::cout << "create vertex-normal map" << std::endl;
    v_nm = make_property_map(FEVV::vertex_normal, m);
    // Store property map in property maps bag.
    put_property_map(FEVV::vertex_normal, m, pmaps_bag, v_nm);
  }
}

/**
 * \brief A mesh type templated main(argc, argv) function that
 *         - loads a mesh from a file,
 *         - applies the \ref progressive_compression_filter generic filter,
 *         - and write the resulting mesh to a file.
 */
template< typename MeshT >
int
progressive_compression_main(int argc, const char **argv)
{
  int mode = 0;

  int number_batches = 15;
  int min_number_vertices = -1;
  int bits_quantization = 12;
  std::string output_file_path_save_preprocess = 
  "progressive_compression_original_mesh_after_preprocess.obj"; // not used in this file
  if(argc < 2)
  {
    std::cout << "Apply a progressive compression filter to the input mesh." << std::endl;
    std::cout << "Usage:  " << argv[0] << "  input_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0]
              << "  ../Testing/Data/CubeNonTriangleFaces.off" << std::endl;
    return EXIT_FAILURE;
  }

  std::string path_binary = "";
  FEVV::Filters::BATCH_CONDITION batch_stop =
      FEVV::Filters::BATCH_CONDITION::ALL_EDGES;
  if(argc > 2)
    mode = atoi(argv[2]);

  if(argc > 3)
  {
    path_binary = argv[3];
    if (path_binary == ".")
    {
      std::cerr << "The binary path cannot be set to \".\". Therefore, it will be set to default value." << std::endl;
      path_binary = "";
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  /////////////////////////////////////////////////////////////////////////////
  std::string output_path = "";
  if (argc > 4)
  {
    output_path = argv[4];
    if (output_path == ".")
    {
      std::cerr << "The output path path cannot be set to \".\". Therefore, it will be set to \"./\"." << std::endl;
      output_path = "./";
    }
  }
  else
  {
    output_path = typeid(m).name();
    if (output_path.find("Surface_mesh") != std::string::npos)
    {
      output_path = "./Surface_mesh_";
    }
    else if (output_path.find("Polyhedron_3") != std::string::npos)
    {
      output_path = "./Polyhedron_3_";
    }
    else if (output_path.find("Linear_cell_complex") != std::string::npos)
    {
      output_path = "./LCC_";
    }
    else if (output_path.find("AIF") != std::string::npos)
    {
      output_path = "./AIF_";
    }
    else
    {
      output_path = "./Other_mesh_type_";
    }
  }

  FEVV::Filters::PREDICTION_TYPE prediction_type =
	    FEVV::Filters::PREDICTION_TYPE::BUTTERFLY;
  FEVV::Filters::METRIC_TYPE metric_type = FEVV::Filters::METRIC_TYPE::QEM_3D;
  FEVV::Filters::VKEPT_POSITION vkeptpos =
      FEVV::Filters::VKEPT_POSITION::HALFEDGE;

  if (argc > 5)
  {
    int predint = atoi(argv[5]);
    prediction_type = static_cast<FEVV::Filters::PREDICTION_TYPE>(predint);
  }
  if (argc > 6)
  {
    int metrint = atoi(argv[6]);
    metric_type = static_cast<FEVV::Filters::METRIC_TYPE>(metrint);
  }
  if (argc > 7)
  {
    int vkeptint = atoi(argv[7]);
    vkeptpos = static_cast<FEVV::Filters::VKEPT_POSITION>(vkeptint);
  }
  if (argc > 8)
    number_batches = atoi(argv[8]);
  if(argc > 9)
    min_number_vertices = atoi(argv[9]);
  if(argc > 10)
  {
    int batchint = atoi(argv[10]);
    batch_stop = static_cast< FEVV::Filters::BATCH_CONDITION >(batchint);
  }
  if(argc > 11)
  {
    bits_quantization = atoi(argv[11]);
  }
  // Input and output files:
  std::string input_file_path = argv[1];
#if(DEBUG)  
  std::string output_file_path = "progressive_compression_filter_base_mesh.obj"; 
#endif
  // Read mesh from file.
  if (mode != 1) {
      FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

      if (min_number_vertices == -1) {
          // set it to 5% of init mesh vertices as in CPM
          min_number_vertices = std::round(FEVV::size_of_vertices(m) * 0.05);
      }
  }

  typename FEVV::PMap_traits< FEVV::vertex_color_t, MeshT >::pmap_type v_cm;
  typename FEVV::PMap_traits< FEVV::edge_color_t, MeshT >::pmap_type e_cm;
  typename FEVV::PMap_traits< FEVV::vertex_normal_t, MeshT >::pmap_type v_nm;

  // Create or retrieve property maps. 
  set_mesh_and_properties(
      m, pmaps_bag, v_cm, e_cm, v_nm);

  // Retrieve point property map (aka geometry)
  auto pm = get(boost::vertex_point, m);

  auto gt_ = FEVV::Geometry_traits< MeshT >(m);
  // Apply filter
  if(mode == 0)
  {
    std::cout << "Single Compression" << std::endl;
    FEVV::Filters::Parameters params(
        prediction_type, vkeptpos, metric_type, true, false, bits_quantization);
    //////////////////////////////////////////////
    // Compression algorithm call (Algorithm 1) //
    //////////////////////////////////////////////	
    progressive_compression_filter(m,
                                   pm,
                                   v_cm,
                                   e_cm,
                                   // v_nm,
                                   gt_,
                                   params,
                                   number_batches,
                                   min_number_vertices,
                                   batch_stop,
                                   output_path,
                                   path_binary,
                                   true,
                                   true,
                                   false, 
                                   output_file_path_save_preprocess);
#if(DEBUG)
     { // Write mesh to file: only write color maps (debug).
       FEVV::PMapsContainer pmaps_bag_empty;
       put_property_map(FEVV::edge_color, m, pmaps_bag_empty, e_cm);
       put_property_map(FEVV::vertex_color, m, pmaps_bag_empty, v_cm);
       FEVV::Filters::write_mesh(output_file_path, m, pmaps_bag_empty);
     }
#endif	 
  }
  else if(mode == 1)
  {
    std::cout << "Measure mode, all metrics being tested" << std::endl;
    std::vector< FEVV::Filters::PREDICTION_TYPE > available_predictions = {
                                FEVV::Filters::PREDICTION_TYPE::BUTTERFLY,
                                FEVV::Filters::PREDICTION_TYPE::DELTA
                                //,FEVV::Filters::PREDICTION_TYPE::POSITION
                           };
    std::vector< FEVV::Filters::VKEPT_POSITION > available_operators = {
                                FEVV::Filters::VKEPT_POSITION::MIDPOINT,
                                FEVV::Filters::VKEPT_POSITION::HALFEDGE
                           };
    std::vector< FEVV::Filters::METRIC_TYPE > available_metrics = {
                                FEVV::Filters::METRIC_TYPE::EDGE_LENGTH,
                                FEVV::Filters::METRIC_TYPE::VOLUME_PRESERVING,
                                FEVV::Filters::METRIC_TYPE::QEM_3D
                           };

    std::vector< int > tested_quantizations = {10, 12, 16};
    
    bool first = true;
    int cpt = 0;
    for(size_t i = 0; i < available_predictions.size(); i++)
    {
      for(size_t j = 0; j < available_operators.size(); j++)
      {
        for(size_t k = 0; k < available_metrics.size(); k++)
        {
          for(size_t l = 0; l < tested_quantizations.size(); l++)
          {
            MeshT mesh;
            number_batches = 70;
            FEVV::PMapsContainer pmaps_bag_second;
            FEVV::Filters::read_mesh(input_file_path, mesh, pmaps_bag_second);

            if (first && (min_number_vertices == -1)) {
                // set it to 5% of init mesh vertices as in CPM
                min_number_vertices = std::round(FEVV::size_of_vertices(mesh) * 0.05);
                first = false;
            }

            typename FEVV::PMap_traits< FEVV::vertex_color_t, MeshT >::pmap_type
                v_cm2;
            typename FEVV::PMap_traits< FEVV::edge_color_t, MeshT >::pmap_type
                e_cm2;
            typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                        MeshT >::pmap_type v_nm2;

            set_mesh_and_properties(mesh,
                                    pmaps_bag_second,
                                    v_cm2,
                                    e_cm2,
                                    v_nm2);


            // Retrieve point property map (aka geometry).
            auto pm2 = get(boost::vertex_point, mesh);

            FEVV::Filters::Parameters params(available_predictions[i],
                                             available_operators[j],
                                             available_metrics[k],
                                             true,
                                             false,
                                             tested_quantizations[l]);
            
            path_binary = ""; // Must be "empty" in order to be set automatically.
            //////////////////////////////////////////////
            // Compression algorithm call (Algorithm 1) //
            //////////////////////////////////////////////			
            progressive_compression_filter(mesh,
                                           pm2,
                                           v_cm2,
                                           e_cm2,
                                           gt_,
                                           params,
                                           number_batches,
                                           min_number_vertices,
                                           batch_stop,
                                           output_path,
                                           path_binary,
                                           true,
                                           true,
                                           false,
                                           output_file_path_save_preprocess);

            ++cpt;
            std::cout << "finished = " << static_cast<double>(cpt)/(available_predictions.size() * available_operators.size() * available_metrics.size() * tested_quantizations.size()) << std::endl;
          }
        }
      }
    }
  }
  else if(mode == 2)
  {
    compute_distortions(m, pm, argv[3]);
  }

  return 0;
}
