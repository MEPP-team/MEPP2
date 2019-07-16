// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

#if 1

// *** MEPP2 code calling CGAL native filters.         ***
// *** Point and normal stored in a CGAL::Point_set_3. ***

#include "FEVV/DataStructures/DataStructures_cgal_point_set.h"
  // for FEVV::CGALPointSet
#include "FEVV/Wrappings/Graph_traits_cgal_point_set.h"
  // for boost::graph_traits< FEVV::CGALPointSet >
#include "FEVV/Wrappings/Geometry_traits_cgal_point_set.h"
  // for FEVV::RetrieveKernel< FEVV::CGALPointSet >
#include "FEVV/Wrappings/Graph_properties_cgal_point_set.h"
  // for get(FEVV::CGALPointSetPointMap&, ...)
#include "FEVV/Wrappings/properties_cgal_point_set.h"
  // for FEVV::PMap_traits< FEVV::CGALPointSet > and vector-property-maps

#include "FEVV/Filters/CGAL/Point_set/cgal_point_set_reader.hpp"
  // for FEVV::Filters::read_mesh< FEVV::CGALPointSet >
#include "FEVV/Filters/CGAL/Point_set/cgal_point_set_writer.hpp"
  // for FEVV::Filters::write_mesh< FEVV::CGALPointSet >

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>


// main
int main(int argc, char *argv[])
{
  // parse arguments
  if(argc != 3)
  {
    std::cout << "Usage:  " << argv[0]
              << "  input_mesh_filename  output_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0] << "  ../Testing/Data/tetra.xyz  tetra.out.ply"
              << std::endl;
    return EXIT_FAILURE;
  }
  std::string input_file = argv[1];
  std::string output_file = argv[2];

  //----------------------------------

  // create point cloud
  typedef  FEVV::CGALPointSet  PointCloudT;
  PointCloudT pc;

  // load point cloud
  std::cout << "load point cloud" << std::endl;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file, pc, pmaps_bag);

  //----------------------------------

  // retrieve Point Map
  auto pm = get(boost::vertex_point, pc);

  // create vertex-normal property map
  using VertexNormalMap = typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                                      PointCloudT >::pmap_type;
  VertexNormalMap v_nm;
  std::cout << "create vertex-normal map" << std::endl;
  v_nm = make_property_map(FEVV::vertex_normal, pc);

  // store property map in property maps bag
  put_property_map(FEVV::vertex_normal, pc, pmaps_bag, v_nm);

  // compute normals
  std::cout << "compute normals" << std::endl;
  typedef typename boost::property_traits< VertexNormalMap >::value_type Normal;
  // Estimates normals direction.
  // Note: pca_estimate_normals() requiresa range of points
  // as well as property maps to access each point's position and normal.
  const int nb_neighbors = 18; // K-nearest neighbors = 3 rings
  CGAL::pca_estimate_normals< CGAL::Sequential_tag >(
      pc,
      nb_neighbors,
      CGAL::parameters::point_map(pm).normal_map(v_nm));

  // Orients normals.
  // Note: mst_orient_normals() requires a range of points
  // as well as property maps to access each point's position and normal.
  std::cout << "orient normals" << std::endl;
  //auto /*indices_range_iterator*/ unoriented_points_begin =
      CGAL::mst_orient_normals(
          pc,
          nb_neighbors,
          CGAL::parameters::point_map(pm).normal_map(v_nm));

  //----------------------------------

  // save point cloud with normals
  std::cout << "save point cloud with normals" << std::endl;
  FEVV::Filters::write_mesh(output_file, pc, pmaps_bag);
  
  return 0;
}

#else

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <utility> // defines std::pair
#include <list>
#include <fstream>

// *** CGAL Native code.                               ***
// *** Point with normal vector stored in a std::pair. ***

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
// Point with normal vector stored in a std::pair.
typedef std::pair< Point, Vector > PointVectorPair;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

int main(int argc, char*argv[])
{
  // parse arguments
  if(argc != 3)
  {
    std::cout << "Usage:  " << argv[0]
              << "  input_mesh_filename  output_mesh_filename" << std::endl;
    std::cout << "Example:  " << argv[0] << "  ../Testing/Data/tetra.xyz  tetra.xyz.with_normals.ply"
              << std::endl;
    return EXIT_FAILURE;
  }
  std::string input_file = argv[1];
  std::string output_file = argv[2];

  // Reads a .xyz point set file in points[].
  std::cout << "load point cloud" << std::endl;
  std::list< PointVectorPair > points;
  std::ifstream stream(input_file);
  if(!stream || !CGAL::read_xyz_points(
                    stream,
                    std::back_inserter(points),
                    CGAL::parameters::point_map(
                        CGAL::First_of_pair_property_map< PointVectorPair >())))
  {
    std::cerr << "Error: cannot read file " << input_file << std::endl;
    return EXIT_FAILURE;
  }

  // Estimates normals direction.
  // Note: pca_estimate_normals() requiresa range of points
  // as well as property maps to access each point's position and normal.
  std::cout << "compute normals" << std::endl;
  const int nb_neighbors = 18; // K-nearest neighbors = 3 rings
  CGAL::pca_estimate_normals< Concurrency_tag >(
      points,
      nb_neighbors,
      CGAL::parameters::point_map(
          CGAL::First_of_pair_property_map< PointVectorPair >())
          .normal_map(CGAL::Second_of_pair_property_map< PointVectorPair >()));

  // Orients normals.
  // Note: mst_orient_normals() requires a range of points
  // as well as property maps to access each point's position and normal.
  std::cout << "orient normals" << std::endl;
  std::list< PointVectorPair >::iterator unoriented_points_begin =
      CGAL::mst_orient_normals(
          points,
          nb_neighbors,
          CGAL::parameters::point_map(
              CGAL::First_of_pair_property_map< PointVectorPair >())
              .normal_map(
                  CGAL::Second_of_pair_property_map< PointVectorPair >()));

  // Optional: delete points with an unoriented normal
  // if you plan to call a reconstruction algorithm that expects oriented
  // normals.
  //points.erase(unoriented_points_begin, points.end());

  // Saves point set.
  // Note: write_xyz_points() requires property maps to access each
  // point position and normal.
  std::cout << "save point cloud with normals" << std::endl;
  std::ofstream out(output_file);
  out.precision(16);
  if(!out ||
     !CGAL::write_xyz_points(
         out,
         points,
         CGAL::parameters::point_map(
             CGAL::First_of_pair_property_map< PointVectorPair >())
             .normal_map(
                 CGAL::Second_of_pair_property_map< PointVectorPair >())))
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

#endif
