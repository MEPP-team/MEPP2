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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/boost/graph/helpers.h> // is_tetrahedron

#include "FEVV/Wrappings/Geometry_traits.h"

#include <iostream>
#include <set>
#include <list>
#include <vector>
#include <tuple> // std::tie
#include <utility>


namespace FEVV {
namespace Filters {

template<
    typename HalfedgeGraph,
    typename PointMap,
    typename VertexColorMap,
    typename vertex_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor,
    typename halfedge_descriptor =
        typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor,
    typename edge_iterator =
        typename boost::graph_traits< HalfedgeGraph >::edge_iterator,
    typename vertex_iterator =
        typename boost::graph_traits< HalfedgeGraph >::vertex_iterator,
    typename Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point,
    typename Geometry = typename FEVV::Geometry_traits< HalfedgeGraph > >

class PreprocessingLD
{
public:
  PreprocessingLD(HalfedgeGraph &g, PointMap &pm, VertexColorMap &cm)
      : _g(g), _pm(pm), _cm(cm)
  {
  }

  //--------------------------input mesh
  //pre-process--------------------------------------------
  void process_mesh_before_quantization()
  {
    // supress all small connected components (isolated vertices and isolated
    // edges)
    this->erase_small_connnected_components(); /// \todo need to be tested

    // does the mesh has borders?
    bool has_border = this->check_if_mesh_with_borders();
    std::clog << "process_mesh_before_quantization: does the mesh has borders ? " << has_border << std::endl;
  }

  void process_mesh_after_quantization()
  {
    // process dupplicates
    this->merge_dupplicates_after_quantization();

    // check if there are not remainong dupplicates
    bool still_dupplicates = this->are_dupplicates();
    if(still_dupplicates)
      std::clog << "WARNING: process_mesh_after_quantization: there are still dupplicates "
                << std::endl;
  }


  /// erase isolated edges and vertices
  void erase_small_connnected_components()
  {
    auto iterator_pair = vertices(_g);
    vertex_iterator vi = iterator_pair.first;
    vertex_iterator vi_end = iterator_pair.second;
    std::vector< vertex_descriptor > isolated_vertices;
    std::vector< halfedge_descriptor > isolated_hedges;

    for(; vi != vi_end; ++vi)
    {

      typename boost::graph_traits< HalfedgeGraph >::degree_size_type
          value_type;
      boost::degree_property_map< HalfedgeGraph > degree_property_map(_g);
      value_type = degree_property_map[*vi];

      if(value_type == 0)
      {
        // isolated vertex case
        isolated_vertices.push_back(*vi);
      }

      if(value_type == 1)
      {
        // isolated edge case
        isolated_hedges.push_back(halfedge(*vi, _g));
      }
    }
    auto itE = isolated_hedges.begin(), 
         itEe = isolated_hedges.end();
    for(; itE != itEe; ++itE)
      remove_edge(edge(*itE, _g), _g);

    auto itV = isolated_vertices.begin(), 
         itVe = isolated_vertices.end();
    for(; itV != itVe; ++itV)
      remove_vertex(*itV, _g);
  }

  bool check_if_mesh_with_borders()
  {
    auto iterator_pair = edges(_g);
    edge_iterator ei = iterator_pair.first;
    edge_iterator ei_end = iterator_pair.second;

    for(; ei != ei_end; ++ei)
    {
      if(CGAL::is_border_edge(halfedge(*ei, _g), _g))
        return true;
    }
    return false;
  }


  void merge_dupplicates_after_quantization()
  {
    auto iterator_pair = vertices(_g);
    vertex_iterator vi = iterator_pair.first;
    vertex_iterator vi_end = iterator_pair.second;
    std::pair< typename std::map< Point,
                                  vertex_descriptor,
                                  std::less< Point > >::iterator,
               bool >
        ret;

    for(; vi != vi_end; ++vi)
    {
      Point pos = get(_pm, *vi);

      // move point until position is not occupied anymore
      while(!(_position.insert(std::pair< Point, vertex_descriptor >(pos, *vi)))
                 .second)
      {
        pos = Point(pos[0], pos[1], pos[2] + 1.0);
        //std::cerr<< "doublon ";
      }
      put(_pm, *vi, pos);
    }
  }

  // simple test to check if all the dupplicates have been removed
  // return true if there are still dupplicates in the mesh and false otherwise
  bool are_dupplicates() const
  {
    auto iterator_pair = vertices(_g);
    vertex_iterator vi = iterator_pair.first;
    vertex_iterator vi_end = iterator_pair.second;

    std::map< Point, vertex_descriptor > map_pos;
    std::pair< typename std::map< Point, vertex_descriptor >::iterator, bool >
        ret;

    for(; vi != vi_end; ++vi)
    {
      const Point& pos = get(_pm, *vi);
      ret = map_pos.insert(std::pair< Point, vertex_descriptor >(pos, *vi));
	  // The ret.second element in the pair is set to true if a new element 
	  // was inserted or false if an equivalent key already existed.
      if(ret.second == false)
        return true;
    }
    return false;
  }


  Point move_dupplicate(vertex_descriptor v, int cas)
  {
    const Point& pos = get(_pm, v);
    Point new_pos;
    Geometry gt(_g);
    int px = static_cast< int >(gt.get_x(pos));
    int py = static_cast< int >(gt.get_y(pos));
    int pz = static_cast< int >(gt.get_z(pos));

    switch(cas)
    {
    case 0:
      new_pos = Point(px, py, pz - 1);
      break;

    case 1:
      new_pos = Point(px, py, pz + 1);
      break;

    case 2:
      new_pos = Point(px, py - 1, pz);
      break;

    case 3:
      new_pos = Point(px, py + 1, pz);
      break;

    case 4:
      new_pos = Point(px - 1, py, pz);
      break;

    case 5:
      new_pos = Point(px + 1, py, pz);
      break;
    }

    return new_pos;
  }


  //--------------------------coarse mesh
  //process--------------------------------------------
  struct P
  {
    int x;
    int y;
    int z;
    bool operator<(const P &q) const
    {
      // std::tie can be used to introduce lexicographical comparison
      return std::tie(x, y, z) < std::tie(q.x, q.y, q.z);
    }
  };

  /////////////////////////////////////////////////////////////////////////////
  void get_new_pos(const std::vector< vertex_descriptor > &doublons,
                   int size_doublons,
                   std::set< P > &new_pos)
  {
    Geometry gt(_g);
    int x_min = gt.get_x(get(_pm, doublons.front())) + 1;
    int x_max = gt.get_x(get(_pm, doublons.back())) - 1;
    int y_min = gt.get_y(get(_pm, doublons.front()));
    int y_max = gt.get_y(get(_pm, doublons.back()));
    int z_min = gt.get_z(get(_pm, doublons.front()));
    int z_max = gt.get_z(get(_pm, doublons.back()));
    int cas = 0;
    int counter = 0;


    if(new_pos.empty())
    {
      P value{gt.get_x(get(_pm, doublons[1])),
              gt.get_y(get(_pm, doublons[1])),
              gt.get_z(get(_pm, doublons[1]))};
      new_pos.insert(value);
    }


    typename std::set< P >::iterator it_new_pos = new_pos.begin();

    while(counter != size_doublons)
    {
      const P& pos = *it_new_pos;
      int x_pos = pos.x;
      int y_pos = pos.y;
      int z_pos = pos.z;
      P value;

      switch(cas)
      {
      case 0:
        if(x_pos + 1 != x_max)
        {
          value = {x_pos + 1, y_pos, z_pos};
          new_pos.insert(value);
          counter += 1;
        }
        break;

      case 1:
        value = {x_pos, y_pos + 1, z_pos};
        new_pos.insert(value);
        counter += 1;
        break;

      case 2:
        value = {x_pos, y_pos, z_pos + 1};
        new_pos.insert(value);
        counter += 1;
        break;

      case 3:
        if(x_pos - 1 != x_min)
        {
          value = {x_pos - 1, y_pos, z_pos};
          new_pos.insert(value);
          counter += 1;
        }
        break;

      case 4:
        value = {x_pos, y_pos - 1, z_pos};
        new_pos.insert(value);
        counter += 1;
        break;

      case 5:
        value = {x_pos, y_pos, z_pos - 1};
        new_pos.insert(value);
        counter += 1;
        break;
      }


      cas += 1;
      cas %= 6;
      if(cas == 0)
        ++it_new_pos;
    }
  }


private:
  HalfedgeGraph &_g;
  PointMap &_pm;
  VertexColorMap &_cm;
  std::map< Point, vertex_descriptor, std::less< Point > > _position;
};
} // namespace Filters
} // namespace FEVV
