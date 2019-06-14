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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <utility> // defines std::pair
#include <list>
#include <fstream>

#if 1

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
  const char *fname = (argc > 1) ? argv[1] : "casting.xyz";

  // Reads a .xyz point set file in points[].
  std::list< PointVectorPair > points;
  std::ifstream stream(fname);
  if(!stream || !CGAL::read_xyz_points(
                    stream,
                    std::back_inserter(points),
                    CGAL::parameters::point_map(
                        CGAL::First_of_pair_property_map< PointVectorPair >())))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // Estimates normals direction.
  // Note: pca_estimate_normals() requiresa range of points
  // as well as property maps to access each point's position and normal.
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
  std::ofstream out(std::string(fname).append(".with_normals.xyz"));
  out.precision(6);
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

#else

// *** Point and normal vector stored in separate containers. ***

// DON'T COMPILE AT THE MOMENT !!!!!!!!!!!!!!!!!!!!
//TODO-elo-fix-it

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

int main(int argc, char*argv[])
{
  const char *fname = (argc > 1) ? argv[1] : "casting.xyz";

  // Reads a .xyz point set file in points[].
  std::vector< Point > points;
  std::ifstream stream(fname);
  if(!stream ||
     !CGAL::read_xyz_points(
         stream,
         std::back_inserter(points),
         CGAL::parameters::point_map(CGAL::Identity_property_map< Point >())))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // Estimates normals direction.
  // Note: pca_estimate_normals() requiresa range of points
  // as well as property maps to access each point's position and normal.
  //std::vector< Vector > normals;
  auto normal_map =
      boost::make_vector_property_map< Vector >(boost::identity_property_map());
  const int nb_neighbors = 18; // K-nearest neighbors = 3 rings
  CGAL::pca_estimate_normals< Concurrency_tag >(
      points,
      nb_neighbors,
      CGAL::parameters::
          point_map(CGAL::Identity_property_map< Point >())
          //.normal_map(CGAL::Second_of_pair_property_map< PointVectorPair >()));
          //TODO-elo how to pass normal prop map here ???????
          .normal_map(normal_map));

  // Orients normals.
  // Note: mst_orient_normals() requires a range of points
  // as well as property maps to access each point's position and normal.
  //std::list< Point >::iterator unoriented_points_begin =
  //    CGAL::mst_orient_normals(
  //        points,
  //        nb_neighbors,
  //        CGAL::parameters::
  //            point_map(CGAL::Identity_property_map< Point >())
  //            .normal_map(normal_map));

  // Optional: delete points with an unoriented normal
  // if you plan to call a reconstruction algorithm that expects oriented
  // normals.
  //points.erase(unoriented_points_begin, points.end());

  // Saves point set.
  // Note: write_xyz_points() requires property maps to access each
  // point position and normal.
  std::ofstream out(std::string(fname).append(".with_normals.xyz"));
  out.precision(6);
  if(!out ||
     !CGAL::write_xyz_points(
         out,
         points,
         CGAL::parameters::
             point_map(CGAL::Identity_property_map< Point >())
             .normal_map(normal_map)))
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

#endif
