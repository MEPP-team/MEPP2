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

#include "FEVV/Wrappings/Graph_traits_pcl_point_cloud.h"
#include <CGAL/boost/graph/properties.h> // for boost::vertex_point_t

#include <stdexcept> // for std::runtime_error


#if defined _MSC_VER
#pragma warning(disable : 4334) // warning C4334 in FLANN
#endif
#include <pcl/kdtree/kdtree_flann.h>
#include <utility> // for std::make_pair
#include <memory>  // for std::unique_ptr
#include <cmath>   // for std::sqrt


namespace FEVV {

// --------------- Point map ---------------

//note: get(boost::vertex_point_t,...) must return a shallow copy of the point
//      map and NOT a reference, so we must create a real PCLPointCloudPointMap
//      class that:
//       - implement shallow copy
//       - implement operator[]
//      see
//        https://github.com/CGAL/cgal/blob/ea20dfd63fcdec0b98258c9c47b0cbb88cdb356c/BGL/include/CGAL/boost/graph/properties_OpenMesh.h#L185
//      and
//        https://www.boost.org/doc/libs/1_57_0/libs/property_map/doc/property_map.html
//        https://www.boost.org/doc/libs/1_57_0/libs/property_map/doc/LvaluePropertyMap.html
class PCLPointCloudPointMap
{
public:
  typedef FEVV::PCLPoint                       value_type;
  typedef FEVV::PCLPoint&                      reference;
  typedef typename boost::graph_traits< FEVV::PCLPointCloud >::vertex_descriptor
      key_type;
  typedef boost::read_write_property_map_tag   category;

  PCLPointCloudPointMap(FEVV::PCLPointCloud &pc) : m_pc(&pc) {}

  PCLPointCloudPointMap(const FEVV::PCLPointCloudPointMap &pm) : m_pc(pm.m_pc)
  {
  }

  value_type operator[](key_type k) const
  {
    // 'points' attribute of pcl::PointCloud<>
    // is a std::vector<FEVV::PCLEnrichedPoint>
    const FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
    return FEVV::PCLPoint(enriched_point.x,
                          enriched_point.y,
                          enriched_point.z);
  }

  void set(key_type k, const FEVV::PCLPoint &point)
  {
    // 'points' attribute of pcl::PointCloud<>
    // is a std::vector<FEVV::PCLEnrichedPoint>
    FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
    enriched_point.x = point[0];
    enriched_point.y = point[1];
    enriched_point.z = point[2];
  }

private:
  FEVV::PCLPointCloud *m_pc; //TODO-elo replace by a reference
};

} // namespace FEVV


namespace pcl {

/**
 * \brief  Returns the points property map (aka the geometry) of the mesh.
 * \note   Must be in 'pcl' namespace because the real type of
 *         FEVV::PCLPointCloud is pcl::PointCloud<...> and the
 *        'get(vertex_point,...)' must be in the same namespace as one
 *         of its parameters.
 */
inline
FEVV::PCLPointCloudPointMap
get(const boost::vertex_point_t, FEVV::PCLPointCloud &pc)
{
  FEVV::PCLPointCloudPointMap pm(pc);
  return pm;
}


} // namespace pcl


namespace FEVV {

//! Specialization of get(point_map, key) for PCLPointCloud
inline
FEVV::PCLPointCloudPointMap::value_type
get(const FEVV::PCLPointCloudPointMap &pm,
    FEVV::PCLPointCloudPointMap::key_type key)
{
  return pm[key];
}

//! Specialization of put(point_map, key, value) for PCLPointCloud
inline
void
put(FEVV::PCLPointCloudPointMap &pm,
    FEVV::PCLPointCloudPointMap::key_type key,
    const FEVV::PCLPointCloudPointMap::value_type &value)
{
  pm.set(key, value);
}

// --------------- Normal map ---------------

// implement a normal map with shallow copy
// normal is stored inside PCLEnrichedPoint structure
class PCLPointCloudNormalMap
{
public:
  typedef FEVV::PCLVector                      value_type;
  typedef FEVV::PCLVector&                     reference;
  typedef typename boost::graph_traits< FEVV::PCLPointCloud >::vertex_descriptor
      key_type;
  typedef boost::read_write_property_map_tag   category;

  PCLPointCloudNormalMap() : m_pc(nullptr) {}

  PCLPointCloudNormalMap(FEVV::PCLPointCloud &pc) : m_pc(&pc) {}

  // shallow copy
  PCLPointCloudNormalMap(const FEVV::PCLPointCloudNormalMap &nm) : m_pc(nm.m_pc)
  {
  }

  value_type operator[](key_type k) const
  {
    // 'points' attribute of pcl::PointCloud<> 
    // is a std::vector<FEVV::PCLEnrichedPoint>
    if(m_pc)
    {
      const FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
      return FEVV::PCLVector(enriched_point.normal_x,
                             enriched_point.normal_y,
                             enriched_point.normal_z);
    }
    else
    {
      throw std::runtime_error(
          "FEVV::PCLPointCloudNormalMap: invalid use of un-initialized map");
      return FEVV::PCLVector(0, 0, 0); // dummy return to avoid a warning
    }
  }

  void set(key_type k, const FEVV::PCLVector &normal)
  {
    // 'points' attribute of pcl::PointCloud<>
    // is a std::vector<FEVV::PCLEnrichedPoint>
    if(m_pc)
    {
      FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
      enriched_point.normal_x = normal[0];
      enriched_point.normal_y = normal[1];
      enriched_point.normal_z = normal[2];
    }
    else
    {
      throw std::runtime_error("Invalid un-initialized PCLPointCloudNormalMap");
    }
  }

private:
  FEVV::PCLPointCloud *m_pc;
};

//! Specialization of get(normal_map, key) for PCLPointCloud
inline
FEVV::PCLPointCloudNormalMap::value_type
get(const FEVV::PCLPointCloudNormalMap &nm,
    FEVV::PCLPointCloudNormalMap::key_type key)
{
  return nm[key];
}

//! Specialization of put(normal_map, key, value) for PCLPointCloud
inline
void
put(FEVV::PCLPointCloudNormalMap &nm,
    FEVV::PCLPointCloudNormalMap::key_type key,
    const FEVV::PCLPointCloudNormalMap::value_type &value)
{
  nm.set(key, value);
}

// --------------- Color map ---------------

// implement a color map with shallow copy
// color is stored inside PCLEnrichedPoint structure
class PCLPointCloudColorMap
{
public:
  typedef FEVV::PCLColor                       value_type;
  typedef FEVV::PCLColor&                      reference;
  typedef typename boost::graph_traits< FEVV::PCLPointCloud >::vertex_descriptor
      key_type;
  typedef boost::read_write_property_map_tag   category;

  PCLPointCloudColorMap() : m_pc(nullptr) {}

  PCLPointCloudColorMap(FEVV::PCLPointCloud &pc) : m_pc(&pc) {}

  // shallow copy
  PCLPointCloudColorMap(const FEVV::PCLPointCloudColorMap &cm) : m_pc(cm.m_pc)
  {
  }

  value_type operator[](key_type k) const
  {
    // 'points' attribute of pcl::PointCloud<> 
    // is a std::vector<FEVV::PCLEnrichedPoint>
    if(m_pc)
    {
      const FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
      return FEVV::PCLColor(enriched_point.r,
                            enriched_point.g,
                            enriched_point.b);
    }
    else
    {
      throw std::runtime_error(
          "FEVV::PCLPointCloudColorMap: invalid use of un-initialized map");
      return FEVV::PCLColor(0, 0, 0); // dummy return to avoid a warning
    }
  }

  void set(key_type k, const FEVV::PCLColor &color)
  {
    // 'points' attribute of pcl::PointCloud<>
    // is a std::vector<FEVV::PCLEnrichedPoint>
    if(m_pc)
    {
      FEVV::PCLEnrichedPoint &enriched_point = m_pc->points[k];
      enriched_point.r = color[0];
      enriched_point.g = color[1];
      enriched_point.b = color[2];
    }
    else
    {
      throw std::runtime_error("Invalid un-initialized PCLPointCloudNormalMap");
    }
  }

private:
  FEVV::PCLPointCloud *m_pc;
};

//! Specialization of get(color_map, key) for PCLPointCloud
inline
FEVV::PCLPointCloudColorMap::value_type
get(const FEVV::PCLPointCloudColorMap &cm,
    FEVV::PCLPointCloudColorMap::key_type key)
{
  return cm[key];
}

//! Specialization of put(color_map, key, value) for PCLPointCloud
inline
void
put(FEVV::PCLPointCloudColorMap &cm,
    FEVV::PCLPointCloudColorMap::key_type key,
    const FEVV::PCLPointCloudColorMap::value_type &value)
{
  cm.set(key, value);
}

} // namespace FEVV


namespace boost {

template<>
struct property_traits< FEVV::PCLPointCloudPointMap >
{
  typedef FEVV::PCLPointCloudPointMap::value_type    value_type;
  typedef FEVV::PCLPointCloudPointMap::reference     reference;
};

} // namespace boost

// --------------- Nearest neighbors search ---------------

namespace pcl {

// NN-search types
struct PCLPointCloudNNSearch
{
  typedef  pcl::KdTreeFLANN< FEVV::PCLEnrichedPoint >  KdTree;
  typedef  std::unique_ptr< KdTree >                   KdTreePtr;
};


// MEPP2 PointCloud concept
//!
//! \brief  Initialize a k-d tree with all cloud points.
//!
inline
PCLPointCloudNNSearch::KdTreePtr
create_kd_tree(const FEVV::PCLPointCloud &pc)
{
  typedef  PCLPointCloudNNSearch::KdTree     KdTree;
  typedef  PCLPointCloudNNSearch::KdTreePtr  KdTreePtr;

  // create k-d tree
  KdTreePtr kd_tree(new KdTree);

  // copy all points in the tree
  kd_tree->setInputCloud(pc.makeShared());
  // note : cloud is passed as a shared pointer here
  //        so all data is copied twice:
  //        1. whole point cloud copied to shared pointer
  //           due to makeShared()
  //        2. all points X,Y,Z copied into kdtree by
  //           setInputCloud()
  
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
std::pair< std::vector< typename boost::graph_traits<
               FEVV::PCLPointCloud >::vertex_descriptor >,
           std::vector< double > >
kNN_search(
    const PCLPointCloudNNSearch::KdTree &kd_tree,
    unsigned int k,
    const FEVV::PCLPoint &query,
    const FEVV::PCLPointCloud &pc) // useless for PCL
{
  typedef typename boost::graph_traits< FEVV::PCLPointCloud >::vertex_descriptor
      vertex_descriptor;

  // prepare query
  FEVV::PCLEnrichedPoint search_point;
  search_point.x = query[0];
  search_point.y = query[1];
  search_point.z = query[2];

  // search K nearest neighbors
  std::vector< int > nn_ids_tmp(k);
  std::vector< float > nn_distances_tmp(k); // squared distance
  int nn_nbr =
      kd_tree.nearestKSearch(search_point, k, nn_ids_tmp, nn_distances_tmp);
  // note : kd_tree.nearestKSearch() wants a std::vector< int > and a
  //        std::vector< float > as outputs, and nothing else ; the more
  //        generic std::vector< vertex_descriptor > and 
  //        std::vector <double > trigger a compilation error as of PCL 1.9.1

  // note : should launch an exception if nn_nbr != k ?

  // convert ids and distances
  std::vector< vertex_descriptor > nn_ids(nn_nbr);
  std::vector< double > nn_distances(nn_nbr);
  for(int i = 0; i < nn_nbr; i++)
  {
    nn_ids[i] = static_cast< vertex_descriptor >(nn_ids_tmp[i]);
    nn_distances[i] = std::sqrt(nn_distances_tmp[i]);
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
std::pair< std::vector< typename boost::graph_traits<
               FEVV::PCLPointCloud >::vertex_descriptor >,
           std::vector< double > >
radius_search(
    const PCLPointCloudNNSearch::KdTree &kd_tree,
    double radius,
    const FEVV::PCLPoint &query,
    const FEVV::PCLPointCloud &pc) // useless for PCL
{
  typedef typename boost::graph_traits< FEVV::PCLPointCloud >::vertex_descriptor
      vertex_descriptor;

  // prepare query
  FEVV::PCLEnrichedPoint search_point;
  search_point.x = query[0];
  search_point.y = query[1];
  search_point.z = query[2];

  // search nearest neighbors
  std::vector< int > nn_ids_tmp;
  std::vector< float > nn_distances_tmp; // squared distance
  int nn_nbr =
      kd_tree.radiusSearch(search_point, radius, nn_ids_tmp, nn_distances_tmp);

  // convert ids and distances
  std::vector< vertex_descriptor > nn_ids(nn_nbr);
  std::vector< double > nn_distances(nn_nbr);
  for(int i = 0; i < nn_nbr; i++)
  {
    nn_ids[i] = static_cast< vertex_descriptor >(nn_ids_tmp[i]);
    nn_distances[i] = std::sqrt(nn_distances_tmp[i]);
  }

  // return pair (vector_of_ids, vector_of_distances)
  return make_pair(nn_ids, nn_distances);
}

} // namespace pcl
