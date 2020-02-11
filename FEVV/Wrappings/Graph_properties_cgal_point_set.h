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
#pragma once

#include "FEVV/Wrappings/Graph_traits_cgal_point_set.h"
#include <CGAL/boost/graph/properties.h> // for boost::vertex_point_t

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <memory> // for std::unique_ptr
#include <utility> // for std::make_pair


namespace FEVV {

using CGALPointSetPointMap = CGALPointSet::Point_map;

} // namespace FEVV


namespace CGAL {

// note:
//   FEVV::CGALPointSet is a typedef to CGAL::Point_set_3<...> ;
//   functions f(FEVV::CGALPointSet) must be declared in the same
//   namespace as CGALPointSet real underlying class for the
//   name lookup mecanism to work properly ;
//   see https://en.cppreference.com/w/cpp/language/adl

/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 */
inline
FEVV::CGALPointSetPointMap
get(const boost::vertex_point_t, FEVV::CGALPointSet &ps)
{
  return ps.point_map();
}

// --------------- Nearest neighbors search ---------------

struct CGALPointSetNNSearch
{
  // k-d tree type
  using vertex_descriptor =
      boost::graph_traits< FEVV::CGALPointSet >::vertex_descriptor;
  using Traits_base = CGAL::Search_traits_3< FEVV::CGALPointSetKernel >;
  using Traits      = CGAL::Search_traits_adapter< vertex_descriptor,
                                              FEVV::CGALPointSetPointMap,
                                              Traits_base >;
  using KdTree      = CGAL::Kd_tree< Traits >;
  using KdTreePtr   = std::unique_ptr< KdTree >;
    // note : CGAL::Kd_Tree can not be returned by a function because its copy
    //        constructor is private ; a compilation error occurs even in case
    //        of copy elision ; so we use a smart pointer as a workaround
};


// MEPP2 PointCloud concept
//!
//! \brief  Initialize a k-d tree with all cloud points.
//!
inline
CGALPointSetNNSearch::KdTreePtr
create_kd_tree(const FEVV::CGALPointSet &pc)
{
  typedef  CGALPointSetNNSearch::Traits     Traits;
  typedef  CGALPointSetNNSearch::KdTree     KdTree;
  typedef  CGALPointSetNNSearch::KdTreePtr  KdTreePtr;

  // retrieve Point map
  auto ppmap = pc.point_map();

  // create kdtree (insert all points in the tree)
  KdTreePtr kd_tree(
      new KdTree(boost::counting_iterator< std::size_t >(0),
                 boost::counting_iterator< std::size_t >(pc.number_of_points()),
                 KdTree::Splitter(),
                 Traits(ppmap)));

  return kd_tree;
}


// MEPP2 PointCloud concept
//!
//! \brief   Search the k nearest neighbors of a geometric point.
//! \return  A pair containing a vector of vertex descriptors
//!          (the ids of the k nearest neighbors) and a vector of distances
//!          (the distance to each nearest neighbor) ; both vectors have size k.
//!
inline
std::pair< std::vector< CGALPointSetNNSearch::vertex_descriptor >,
           std::vector< double > >
kNN_search(
    const CGALPointSetNNSearch::KdTree &kd_tree,
    unsigned int k,
    const FEVV::CGALPointSetPoint &query,
    const FEVV::CGALPointSet& pc)
{
  typedef CGALPointSetNNSearch::vertex_descriptor       vertex_descriptor;
  typedef CGALPointSetNNSearch::Traits                  Traits;
  typedef CGAL::Orthogonal_k_neighbor_search< Traits >  K_neighbor_search;
  typedef K_neighbor_search::Distance                   Distance;

  // retrieve Point map
  auto ppmap = pc.point_map();

  // search K nearest neighbors
  Distance sq_dist(ppmap); // squared distance
  K_neighbor_search search(kd_tree, query, k, 0, true, sq_dist);

  // convert ids and distances
  std::vector< vertex_descriptor > nn_ids;
  std::vector< double > nn_distances;
  for(auto it = search.begin(); it != search.end(); it++)
  {
    nn_ids.push_back(it->first);
    nn_distances.push_back(sq_dist.inverse_of_transformed_distance(it->second));
  }

  // return pair (vector_of_ids, vector_of_distances)
  return make_pair(nn_ids, nn_distances);
}


// MEPP2 PointCloud concept
//!
//! \brief   Search for the nearest neighbors of a geometric point in the
//!          given radius.
//! \return  A pair containing a vector of vertex descriptors
//!          (the ids of the nearest neighbors) and a vector of distances
//!          (the distance to each nearest neighbor) ; both vectors have
//!          the same size.
//!
inline
std::pair< std::vector< CGALPointSetNNSearch::vertex_descriptor >,
           std::vector< double > >
radius_search(
    const CGALPointSetNNSearch::KdTree &kd_tree,
    double radius,
    const FEVV::CGALPointSetPoint &query,
    const FEVV::CGALPointSet& pc)
{
  typedef CGALPointSetNNSearch::vertex_descriptor       vertex_descriptor;
  typedef CGALPointSetNNSearch::Traits                  Traits;
  typedef CGAL::Orthogonal_incremental_neighbor_search< Traits >
                                                        NN_incremental_search;
  typedef NN_incremental_search::Distance               Distance;

  // retrieve Point map
  auto ppmap = pc.point_map();

  // search nearest neighbors
  Distance sq_dist(ppmap); // squared distance
  NN_incremental_search search(kd_tree, query, 0.0, true, sq_dist);

  // get nearest neighbors ids with distance < radius
  std::vector< vertex_descriptor > nn_ids;
  std::vector< double > nn_distances;
  for(auto it = search.begin(); it != search.end(); ++it)
  {
    double distance = sq_dist.inverse_of_transformed_distance(it->second);
    if(distance > radius)
      break; // stop searching

    // add current point to list
    nn_ids.push_back(it->first);
    nn_distances.push_back(distance);
  }

  // return pair (vector_of_ids, vector_of_distances)
  return make_pair(nn_ids, nn_distances);
}

} // namespace CGAL


namespace boost {

template<>
struct property_traits< FEVV::CGALPointSet >
{
  typedef FEVV::CGALPointSetPointMap::value_type    value_type;
  typedef FEVV::CGALPointSetPointMap::reference     reference;
};

} // namespace boost
