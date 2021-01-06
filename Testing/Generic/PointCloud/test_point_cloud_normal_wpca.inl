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

#include "FEVV/Filters/Generic/PointCloud/WeightedPCANormals/point_cloud_normal_wpca.hpp"

#include "Testing/Utils/utils_identical_text_based_files.hpp"

#include <string>

/**
 * \brief  Generic test of Weighted PCA filter.
 */
template< typename PointCloud >
int test_point_cloud_normal_wpca(int argc, const char **argv)
{
    if (argc != 6 && argc != 7)
    {
        std::cout<<std::endl;
        printf("./normal_estimation [input pointcloud] [output_pointcloud.ply]  [number of neighbors] [noise] [curvature] [OPTIONNAL:reference_file]\n");
        printf("Estimated noise (m),\n \t example : 0.0001 \n");
        printf("Min tolerated curvature radius (m) \n \t example for curved objects = 0.1 \n \t for intersecting planes = 10000000000\n");
        printf("It is recommended to scale the objects (for diagonal to be 1) to better parameterize the min tolerated curvature radius\n");

        std::cout << "Examples: \n"
                  << argv[0] << "  casting.xyz  casting.curvature.ply 10 0 0.01 casting.curvature.ref.ply\n"
                  << std::endl;
        return EXIT_FAILURE;
    }
    else
    {
        // display parameters summary
        std::cout << "\nParameters summary:" << std::endl;
        std::cout<<std::endl;
        for (int i = 1 ; i<argc; i++)
            std::cout<<argv[i]<<"   ";
        std::cout<<std::endl<<std::endl;
    }

    // input and output files
    std::string input_file_path = argv[1];
    std::string output_file_path = argv[2];

    // read mesh from file
    PointCloud cloud;
    FEVV::PMapsContainer pmaps_bag;
    FEVV::Filters::read_mesh(input_file_path, cloud, pmaps_bag);

    // Note: the property maps must be extracted from the
    //       property maps bag, and explicitely passed as
    //       parameters to the filter, in order to make
    //       clear what property is used by the filter

    // retrieve or create vertex-normal property map
    using VertexNormalMap =
        typename FEVV::PMap_traits< FEVV::vertex_normal_t, PointCloud >::pmap_type;
    VertexNormalMap v_nm;
    if(has_map(pmaps_bag, FEVV::vertex_normal))
    {
        std::cout << "use existing vertex-normal map" << std::endl;
        v_nm = get_property_map(FEVV::vertex_normal, cloud, pmaps_bag);
    }
    else
    {
        std::cout << "create vertex-normal map" << std::endl;
        v_nm = make_property_map(FEVV::vertex_normal, cloud);
        // store property map in property maps bag
        put_property_map(FEVV::vertex_normal, cloud, pmaps_bag, v_nm);
    }

    // retrieve point property map (aka geometry)
    auto pm = get(boost::vertex_point, cloud);

    // apply filter
    int n_neigh = static_cast<int>(atoi(argv[3]));
    float noise = static_cast<float>(atof(argv[4])/sqrt(3));
    noise = std::max(noise, 0.000001f);
    float curvature = static_cast<float>(atof(argv[5]));
    FEVV::Filters::compute_weighted_pca(cloud, pm, v_nm, n_neigh, noise, curvature);

    //std::string output_file_path = create_output_name(input_file_path, "weightedPCA", n_neigh, "pcd");

    // save the point cloud
    std::cout << "saving point cloud with normals to " << output_file_path;
    FEVV::Filters::write_mesh(output_file_path, cloud, pmaps_bag);
    
    // check output file
    std::string reference_file_path;
    if(argc == 7)
        reference_file_path = argv[6];

    if(!reference_file_path.empty())
    {
      std::cout << "Comparing output file" << std::endl;
      std::cout << "  '" << output_file_path << "'" << std::endl;
      std::cout << "with reference file" << std::endl;
      std::cout << "  '" << reference_file_path << "'" << std::endl;
      std::cout << "..." << std::endl;

      // use text file comparator
      if(!identical_text_based_files(output_file_path, reference_file_path))
      {
        std::cout << "Files are different!" << std::endl;
        return EXIT_FAILURE;
      }

      std::cout << "Files are identical." << std::endl;
    }

    return 0;
}
