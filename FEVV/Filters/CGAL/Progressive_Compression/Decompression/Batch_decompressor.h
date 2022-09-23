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

#include "FEVV/Filters/Generic/generic_writer.hpp"

#include "CGAL/boost/graph/Euler_operations.h"


#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/Edge_length_metric.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/Volume_preserving.h"

#include "FEVV/Tools/Comparator/Spanning_tree_vertex_edge_comparator.hpp"

#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Memory_comparator.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Raw_positions.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Butterfly.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/apply_color.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Uniform_dequantization.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Helpers/Header_handler.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Geometric_metrics.h" // needed for the compilation

#include "FEVV/Filters/CGAL/Progressive_Compression/Decompression/Binary_batch_decoder.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Decompression/Coarse_mesh_decoder.h"

#include "FEVV/DataStructures/DataStructures_cgal_linear_cell_complex.h" // to test is_isomorphic_to method
#include "FEVV/Wrappings/Geometry_traits_cgal_linear_cell_complex.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"


#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <utility>

namespace FEVV {
namespace Filters {

enum class TOPOLOGY_CASE {

  CASE11 = 0,
  CASE12,
  CASE13,
  CASE211,
  CASE212,
  CASE221,
  CASE222,
  SIMPLE
};



/**
 * \brief Batch_decompressor: Given a draco buffer, will decompress a 
 *        non-textured (or not) mesh.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap >
class Batch_decompressor
{


public:
  using vertex_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor;
  using halfedge_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor;
  using face_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::face_descriptor;
  using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;
  using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
  using Geometry = typename FEVV::Geometry_traits< HalfedgeGraph >;

  /**
   * \brief Batch_decompressor
   * @param[in,out] g Empty halfedge graph that is used to reconstruct the
   *                decoded mesh topology.
   * @param[in,out] pm Empty pointmap associated with g that is used to
   *                reconstruct the decoded mesh geometry (vertex positions).
   * @param[in] predictor Pointer to a predictor object (e.g. butterfly).
   * @param[in] vkept Pointer to a mid midpoint or halfedge object.
   * @param[in] header Header object (contains refinement info: quantization 
   *            for example).
   * @param[in,out] vcm vertex color map associated with g, used to debug
   *                topological issues during the decoding   
   * @param[in] bit_quantization bits of quantization
  */
  Batch_decompressor(
      HalfedgeGraph &g,
      PointMap &pm,
      Predictor< HalfedgeGraph, PointMap > *predictor,
      FEVV::Filters::Kept_position< HalfedgeGraph,
                                   PointMap,
                                   Geometry > *vkept,
      const FEVV::Header_handler &header,
      VertexColorMap &vcm)
      : _g(g), _pm(pm), _gt(Geometry(g)), _predictor(predictor), _vkept(vkept),
        _header(header), _vcm(vcm)
  {
    _batch_id = 0;
  }
  ~Batch_decompressor() {}

  /**
   * \brief predict_positions: Given a list of halfedges, will predict every 
   *        pair of vertices of a split.
   *        This function implements line 13 of Algorithm 2.
   * @param[in] h_extent List of halfedges to expand to faces.
   * @param[out] new_points Empty list of pair<Point,Point>. Will be filled by
   *             the function.
   */
  void predict_positions(const std::list< halfedge_descriptor > &h_extent,
                        std::list< std::pair< Point, Point > > &new_points)
  {

    auto pos_it = _positions.begin();
    auto info_it = _other_info_bits.begin(); // reverse information
    for(typename std::list< halfedge_descriptor >::const_iterator h_it =
            h_extent.begin();
        h_it != h_extent.end();
        h_it++)
    {
      halfedge_descriptor hnew;
      halfedge_descriptor h1 = *h_it;
      h_it++;
      halfedge_descriptor h2 = *h_it;

      vertex_descriptor vkept = target(h1, _g);

      if(_vkept->get_type() == VKEPT_POSITION::HALFEDGE)
      {
        bool b = (*info_it);

        _predictor->set_rev(b);
        info_it++;
      }

      std::pair< Point, Point > resulting_points =
          _predictor->place_points(*pos_it, vkept, h1, h2);
      pos_it++;
      new_points.push_back(std::move(resulting_points));
    }
  }

  /// saves spanning tree (useful for debugging)
  void save_spanning_tree(
      const FEVV::Comparator::Spanning_tree_vertex_edge_comparator< HalfedgeGraph,
                                                          PointMap >
          &spanningtree)
  {
    std::ofstream file_encode;
    file_encode.open("decodingspanning" + std::to_string(_batch_id) + ".txt");
    auto span = spanningtree.get_spanning_tree_vertices();
    typename std::list< vertex_descriptor >::const_iterator it = span.begin();

    for( ; it != span.end(); it++)
    {
      const Point& pt = get(_pm, *it);
      std::list< vertex_descriptor > list_v_around =
          FEVV::Comparator::get_adjacent_vertices(*it, _g);
      file_encode << _gt.get_x(pt) << " " << _gt.get_y(pt) << " "
                  << _gt.get_z(pt) << " " << list_v_around.size() << std::endl;
    }
    file_encode.close();
  }

  /**
   * \brief Dequantizes a mesh (restores integers as original doubles)
   */
  void dequantize_mesh()
  {
    Uniform_dequantization< HalfedgeGraph, PointMap > dq(
        _g,
        _pm,
        _header.get_quantization(),
        _header.get_dimension(),
        _header.get_init_coord());

    dq.point_dequantization();
  }
 
  /** 
   * \brief Decompresses a batch using a draco buffer.
   *  This function implements lines 5 to 19 of Algorithm 2.
   * @param[in,out] buffer A draco decoderbuffer with refinement info.
   *            Not a const param because of call to non-const draco 
   *            methods such as Decode.
   */
  void decompress_binary_batch(draco::DecoderBuffer &buffer)
  {
    // Decode data for a batch.
    Binary_batch_decoder< HalfedgeGraph,
                        PointMap,
                        Vector >
        decoder(buffer, _header.get_quantization());
    decoder.decode_bitmask(_bitmask); // Implements line 5 of Algorithm 2.
    decoder.decode_bitmask(_edge_bitmask); // Implements line 6 of Algorithm 2.

    // Implements lines 7 to 9 of Algorithm 2.
    if (_vkept->get_type() == FEVV::Filters::VKEPT_POSITION::HALFEDGE)
    { // For HALFEDGE/target case, an additional bit information is
      // needed
      decoder.decode_bitmask(_other_info_bits);  // Implements line 8 of Algorithm 2.
    }

    // Implements line 10 of Algorithm 2.
    decoder.decode_residuals(_positions, _predictor->get_nb_residuals());
    ///////////////////////////////////////////////////////////////////////////
    // Predict position and refine mesh from the decoded data

    // Build a spanning tree for halfedge graph _g with vertex positions in _pm
    // Implements line 11 of Algorithm 2.
    FEVV::Comparator::
        Spanning_tree_vertex_edge_comparator< HalfedgeGraph, PointMap, Geometry >
            spanningtree(_g, _pm, true);// "true"(=tie-break detection and avoidance)
                                        // as for Batch_collapser to get
                                        // the same adjacency ordering (and not
                                        // too costly for reasonable 
                                        // quantization)
    ///////////////////////////////////////////////////////////////////////////
    std::list< halfedge_descriptor > h_extent;
    // Find every halfedge we have to expand into faces.
    fill_h_extent_list_no_sort(spanningtree, h_extent); // Implements line 12 of 
                                                        // Algorithm 2.
    // Split corresponding vertices
    split_vertices(h_extent); // Implements lines 13 to 19 of Algorithm 2.
		
    std::cout << "batch " << std::to_string(_batch_id) << " done" << std::endl;
    _batch_id++;
    _bitmask.clear();
    _edge_bitmask.clear();
    _positions.clear();
    _other_info_bits.clear();
  }

private:
  HalfedgeGraph &_g;
  PointMap &_pm;
  const Geometry _gt;
  std::vector< std::vector< bool > > neighbours_bitmasks;
  int _batch_id;

  FEVV::Filters::
      Predictor< HalfedgeGraph, PointMap >
          *_predictor;
  FEVV::Filters::Kept_position< HalfedgeGraph,
                               PointMap,
                               Geometry > *_vkept;

  std::list< std::vector< Vector > > _positions;

  std::list< bool > _bitmask;
  std::list< bool > _edge_bitmask;
  std::list< bool > _other_info_bits;
  const FEVV::Header_handler &_header; // quantization info is obtained via the 
                                       // header
  VertexColorMap &_vcm;

  /// This function implements lines 13 to 19 of Algorithm 2.
  bool split_vertices(const std::list< halfedge_descriptor > &h_extent)
  {
    ///////////////////////////////////////////////////////
    // DECODE VERTEX POSITIONS FOR ALL VERTICES TO SPLIT //
    ///////////////////////////////////////////////////////
    std::list< std::pair< Point, Point > > new_points;
    predict_positions(h_extent, new_points); // Implements line 13 of Algorithm 2.
    auto pos_it = _positions.begin();
    auto info_it = _other_info_bits.begin();
    /////////////////////////
    // APPLY VERTEX SPLITS //
    /////////////////////////
    // Implements lines 14 to 19 of Algorithm 2: 
    auto new_points_it = new_points.begin();
    for(typename std::list< halfedge_descriptor >::const_iterator h_it =
            h_extent.begin();
        h_it != h_extent.end();
        h_it++)
    {
      halfedge_descriptor hnew;
      halfedge_descriptor h1 = *h_it;
      h_it++;
      halfedge_descriptor h2 = *h_it;

      TOPOLOGY_CASE current_case = find_case(h1, h2);

      bool reverse = false;
      if(_vkept->get_type() == VKEPT_POSITION::HALFEDGE)
      {
        bool b = (*info_it);

        _predictor->set_rev(b);
        reverse = b;
        info_it++;
      }

      std::pair< Point, Point > resulting_points = *new_points_it;
      new_points_it++;

      if(current_case == TOPOLOGY_CASE::SIMPLE)
      {
        hnew = simple_split_sase(h1, h2);
      }
      else
      {
        if(current_case == TOPOLOGY_CASE::CASE211 ||
           current_case == TOPOLOGY_CASE::CASE212 ||
           current_case == TOPOLOGY_CASE::CASE221 ||
           current_case == TOPOLOGY_CASE::CASE222)
        {
          hnew = cases_2(h1, h2, current_case);
        }
        else
        {

          hnew = cases_1(h1, current_case);
        }
      }
      put(_pm, target(hnew, _g), resulting_points.first);
      put(_pm, source(hnew, _g), resulting_points.second);
      
      if(pos_it != _positions.end())
      {
        pos_it++;
      }
    }

    return true;
  }

  /// Fills the list of halfedges to expand. For cases without any border, two
  /// different halfedges will be added. Otherwise, one halfedge will be added
  /// twice.
  /// This function implements line 12 of Algorithm 2.
  void fill_h_extent_list_no_sort(
      const FEVV::Comparator::
          Spanning_tree_vertex_edge_comparator< HalfedgeGraph, PointMap, Geometry >
              &spanningtree,
      std::list< halfedge_descriptor > &h_extent)

  {
    const std::list< vertex_descriptor >& spanning_vertices =
        spanningtree.get_spanning_tree_vertices();

    // get iterators on topology bitmasks

    // Bitmask indicating whether we split a vertex or not
    std::list< bool >::iterator bit_it = _bitmask.begin();
    std::list< bool >::iterator bit_it_end = _bitmask.end();


    // Btmask indicating the halfedges to expand around a vertex to be split
    std::list< bool >::iterator edge_bit_it = _edge_bitmask.begin();
    std::list< bool >::iterator edge_bit_it_end = _edge_bitmask.end();

    typename std::list< vertex_descriptor >::const_iterator spanning_it =
        spanning_vertices.begin();
    std::list< vertex_descriptor > vsplit_neighbours_list;
		
	// bit optimization stuff (to not encode predictable zeros)
	bool last_was_1 = false;	
    std::list< vertex_descriptor > remaining_adjacent_vertices_of_last_v, tmp;
    std::map<vertex_descriptor, bool> processed_vertices;
		
    // Iterate through the vertices in the order of the spanning tree
    for(; bit_it != bit_it_end; spanning_it++)
    {
      // bit optimization stuff (to not encode predictable zeros)
      processed_vertices[*spanning_it] = true;
      bool is_adjacent_to_former_1 = false;
      if(last_was_1)
      {
        auto it_to_remove = std::find(remaining_adjacent_vertices_of_last_v.begin(), 
                                      remaining_adjacent_vertices_of_last_v.end(), *spanning_it) ;
        is_adjacent_to_former_1 = (it_to_remove != remaining_adjacent_vertices_of_last_v.end());
        if(is_adjacent_to_former_1)
          remaining_adjacent_vertices_of_last_v.erase(it_to_remove);
        if(remaining_adjacent_vertices_of_last_v.empty())
          last_was_1 = false;
      }	  
		
      // if we have to split the current vertex
      if((*bit_it == true) && !is_adjacent_to_former_1)
      {
        // bit optimization stuff 
        last_was_1 = true; 
        tmp = FEVV::Comparator::get_not_processed_adjacent_vertices(*spanning_it, _g, processed_vertices, spanningtree.get_spanning_tree_min_incident_edge(*spanning_it));
        remaining_adjacent_vertices_of_last_v.insert(remaining_adjacent_vertices_of_last_v.end(), tmp.begin(), tmp.end()); 
		  
        // find min edge
        std::vector< bool > current_edge_bit;
        auto e_min_tree = spanningtree.get_spanning_tree_min_incident_edge(*spanning_it);
        halfedge_descriptor h_min_tree = halfedge(e_min_tree, _g);
        if(target(h_min_tree, _g) != *spanning_it)
          h_min_tree = opposite(h_min_tree, _g);
        halfedge_descriptor current_edge = h_min_tree;
        int nb_edges = 0;
        do
        {
          // find the halfedges to extent (0: don't extend, 1: extend)
          if(nb_edges>1)
            current_edge_bit.push_back(false);
          else
          {
            current_edge_bit.push_back(*edge_bit_it);
            if(*edge_bit_it == true)
            {
              // if 1, we add the halfedge to the list
              h_extent.push_back(current_edge);
              nb_edges += 1;
            }
            edge_bit_it++;
          }
          // run around the current vertex
          current_edge = opposite(next(current_edge, _g), _g);
        } while(current_edge != h_min_tree &&
                edge_bit_it != edge_bit_it_end);
        neighbours_bitmasks.push_back(current_edge_bit);
        if(nb_edges == 1)
        {
          h_extent.push_back(h_extent.back());
        }

        if(nb_edges != 1 && nb_edges != 2)
        {
          std::cerr << "fill_h_extent_list_no_sort: error. No halfedge to extent." << std::endl;
        }

        vsplit_neighbours_list.clear();
		
        bit_it++;
      }
	  else if(!is_adjacent_to_former_1)
		  bit_it++;
    }
  }
 
  FEVV::Filters::TOPOLOGY_CASE find_case(halfedge_descriptor h1,
                                        halfedge_descriptor h2)
  {
    if(h1 == h2)
    {
      halfedge_descriptor h = h1;
      halfedge_descriptor opp_h = opposite(h, _g);
      if(!(CGAL::is_border(h, _g)) &&
         !(!CGAL::is_border(
             opp_h,
             _g))) // both halfedges are not border, but a vertex might be
      {
        return TOPOLOGY_CASE::CASE11;
      }
      else
      {
        if(!CGAL::is_border(h, _g) &&
           CGAL::is_border(opp_h, _g)) // opposite on border
        {
          return TOPOLOGY_CASE::CASE12;
        }
        else
        {
          return TOPOLOGY_CASE::CASE13;
        }
      }
    }
    else
    {
      if((face(h1, _g) != boost::graph_traits< HalfedgeGraph >::null_face()) &&
         (face(h2, _g) !=
          boost::graph_traits<
              HalfedgeGraph >::null_face())) // both edges completely internal
      {
        return TOPOLOGY_CASE::SIMPLE;
      }
      else
      {
        if(face(h1, _g) == boost::graph_traits< HalfedgeGraph >::null_face())
        {
          if(face(opposite(h2, _g), _g) ==
             boost::graph_traits< HalfedgeGraph >::null_face())
          {
            return TOPOLOGY_CASE::CASE211;
          }
          else
          {
            return TOPOLOGY_CASE::CASE221;
          }
        }
        else
        {
          if(face(opposite(h1, _g), _g) ==
             boost::graph_traits< HalfedgeGraph >::null_face())
          {
            return TOPOLOGY_CASE::CASE212;
          }
          else
          {
            return TOPOLOGY_CASE::CASE222;
          }
        }
      }
    }
  }

  halfedge_descriptor simple_split_sase(
      halfedge_descriptor h1,
      halfedge_descriptor h2) // Simple Split Case: both h1 and h2 are internal
                              // to the mesh, and none is linked to a null face
  {
    if((face(h1, _g) != boost::graph_traits< HalfedgeGraph >::null_face()) &&
       (face(h2, _g) != boost::graph_traits< HalfedgeGraph >::null_face()))
    {
      halfedge_descriptor hnew =
          CGAL::Euler::split_vertex< HalfedgeGraph >(h1, h2, _g);
      CGAL::Euler::split_face(
          prev(prev(opposite(hnew, _g), _g), _g), opposite(hnew, _g), _g);
      CGAL::Euler::split_face(hnew, next(next(hnew, _g), _g), _g);

      return hnew;
    }
    return h1;
  }


  std::pair< halfedge_descriptor, bool > find_good_border(halfedge_descriptor h)
  {
    halfedge_descriptor hb = h;
    bool ok = false;
    do
    {
      hb = opposite(next(hb, _g), _g);
	  // check if good border or if it's a hole in the mesh
      if(CGAL::is_border(hb, _g))
      {
        ok = true;
      }
      else
      {

        if(CGAL::is_border(opposite(hb, _g), _g))
        {
          hb = opposite(hb, _g);
          ok = true;
        }
      }
    } while((hb != h) && !ok);

    return std::make_pair(hb, ok);
  }

  std::pair< halfedge_descriptor, halfedge_descriptor >
  find_good_halfedges_for_split(halfedge_descriptor h, TOPOLOGY_CASE cas)
  {
    halfedge_descriptor h1, h2;
    if(cas == TOPOLOGY_CASE::CASE11)
    {
      // case 1
      h2 = h;
      // find good border
      std::pair< halfedge_descriptor, bool > paire = this->find_good_border(h2);
      if(paire.second == false)
      {
        FEVV::PMapsContainer test_pmaps_bag;
        color_vertex(
            _g,
            _pm,
            _vcm,
            target(h2, _g),
            typename boost::property_traits< VertexColorMap >::value_type(
                1, 0, 0));
        put_property_map(FEVV::vertex_color, _g, test_pmaps_bag, _vcm);
        FEVV::Filters::write_mesh("test11.obj", _g, test_pmaps_bag);
        std::cerr << "find_good_halfedges_for_split: error, there is no good border found " << std::endl;
        assert(false);
      }
      else
        h1 = paire.first;
    }

    else if(cas == TOPOLOGY_CASE::CASE12)
    {
      // case 2
      h2 = h;
      // find good border
      std::pair< halfedge_descriptor, bool > paire = this->find_good_border(h2);
      if(paire.second == false)
      {
        FEVV::PMapsContainer test_pmaps_bag;
        color_vertex(
            _g,
            _pm,
            _vcm,
            target(h2, _g),
            typename boost::property_traits< VertexColorMap >::value_type(
                1, 0, 0));
        put_property_map(FEVV::vertex_color, _g, test_pmaps_bag, _vcm);
        FEVV::Filters::write_mesh("test12.obj", _g, test_pmaps_bag);
        std::cerr << "find_good_halfedges_for_split: error, there is no good border found " << std::endl;
        assert(false);
      }
      else
        h1 = paire.first;
    }

    else if(cas == TOPOLOGY_CASE::CASE13)
    {
      // case 3
      h2 = prev(opposite(h, _g), _g);
      // find good border
      std::pair< halfedge_descriptor, bool > paire = this->find_good_border(h2);
      if(paire.second == false)
      {
        FEVV::PMapsContainer test_pmaps_bag;
        color_vertex(
            _g,
            _pm,
            _vcm,
            target(h, _g),
            typename boost::property_traits< VertexColorMap >::value_type(
                1, 0, 0));
        color_vertex(
            _g,
            _pm,
            _vcm,
            target(h2, _g),
            typename boost::property_traits< VertexColorMap >::value_type(
                1, 1, 0));
        put_property_map(FEVV::vertex_color, _g, test_pmaps_bag, _vcm);
        FEVV::Filters::write_mesh("test13.obj", _g, test_pmaps_bag);
        std::cerr << "find_good_halfedges_for_split: error, there is no good border found " << std::endl;
        assert(false);
      }
      else
        h1 = paire.first;
    }

    else
    {
      std::cerr << "find_good_halfedges_for_split: error , wrong case " << std::endl;
      assert(false);
    }

    return std::make_pair(h1, h2);
  }


  int get_face_degree(halfedge_descriptor h)
  {

    face_descriptor f = face(h, _g);
    if(f == boost::graph_traits< HalfedgeGraph >::null_face())
    {
      std::cerr << "get_face_degree: null face " << std::endl;
      return 0;
    }
    else
    {
      boost::iterator_range<
          CGAL::Vertex_around_face_iterator< HalfedgeGraph > >
          iterator_range_v = CGAL::vertices_around_face(h, _g);
      int nb_face_vertices = 0;

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif //__GNUC__
      for(auto v : iterator_range_v)
      {
        nb_face_vertices += 1;
      }
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif //__GNUC__

      if(nb_face_vertices == 1 || nb_face_vertices == 2)
      {
        std::cerr << "get_face_degree: error, there are not enough vertices to make a face "
                  << std::endl;
        assert(false);
        return 0;
      }
      else
        return nb_face_vertices;
    }
  }
  halfedge_descriptor cases_1(halfedge_descriptor h, TOPOLOGY_CASE current_case)
  {
    halfedge_descriptor hnew, h1, h2;
    std::pair< halfedge_descriptor, halfedge_descriptor > pair_of_h;
    switch(current_case)
    {
    case FEVV::Filters::TOPOLOGY_CASE::CASE11:
      pair_of_h = find_good_halfedges_for_split(h, current_case);
      h1 = pair_of_h.first;
      h2 = pair_of_h.second;

      hnew = CGAL::Euler::split_vertex(h1, h2, _g);
      if(this->get_face_degree(hnew) == 4)
        CGAL::Euler::split_face(hnew, next(next(hnew, _g), _g), _g);
      break;

    case FEVV::Filters::TOPOLOGY_CASE::CASE12:
      pair_of_h = find_good_halfedges_for_split(h, current_case);
      h1 = pair_of_h.first;
      h2 = pair_of_h.second;
      hnew = CGAL::Euler::split_vertex(h1, h2, _g);

      CGAL::internal::set_border(opposite(hnew, _g), _g);

      CGAL::Euler::split_face(hnew, next(next(hnew, _g), _g), _g);
      break;

    case FEVV::Filters::TOPOLOGY_CASE::CASE13:

      pair_of_h = find_good_halfedges_for_split(h, current_case);
      h1 = pair_of_h.first;
      h2 = pair_of_h.second;
      if(h1 != h2)
      {

        hnew = CGAL::Euler::split_vertex(h1, h2, _g);

        CGAL::internal::set_border(opposite(hnew, _g), _g);

        CGAL::Euler::split_face(prev(hnew, _g), next(hnew, _g), _g);
      }
      else
      {
        color_vertex(
            _g,
            _pm,
            _vcm,
            target(h, _g),
            typename boost::property_traits< VertexColorMap >::value_type(
                1, 0, 0));
        hnew = h;
      }

      break;

    default:
      break;
    }

    return hnew;
  }


  halfedge_descriptor cases_2(halfedge_descriptor h1,
                             halfedge_descriptor h2,
                             TOPOLOGY_CASE current_case)
  {
    halfedge_descriptor first_to_split, second_to_split, border;
    vertex_descriptor first, second, third;
    std::vector< vertex_descriptor > face_vertices;
    face_vertices.reserve(3);
    halfedge_descriptor hnew = CGAL::Euler::split_vertex(h1, h2, _g);
    switch(current_case)
    {

    case FEVV::Filters::TOPOLOGY_CASE::CASE211:
      first_to_split = hnew;
      second_to_split = prev(prev(hnew, _g), _g);
      border = opposite(hnew, _g);
      first = source(hnew, _g);
      second = source(h1, _g);
      third = target(hnew, _g);
      break;
    case FEVV::Filters::TOPOLOGY_CASE::CASE221:
      first_to_split = hnew;
      second_to_split = prev(prev(hnew, _g), _g);
      border = opposite(hnew, _g);
      first = source(hnew, _g);
      second = source(h1, _g);
      third = target(hnew, _g);

      break;

    case FEVV::Filters::TOPOLOGY_CASE::CASE212:
      first_to_split = opposite(hnew, _g);
      second_to_split = prev(prev(opposite(hnew, _g), _g), _g);
      border = hnew;
      first = target(hnew, _g);
      second = source(h2, _g);
      third = source(hnew, _g);
      break;
    case FEVV::Filters::TOPOLOGY_CASE::CASE222:
      first_to_split = opposite(hnew, _g);
      second_to_split = prev(prev(opposite(hnew, _g), _g), _g);
      border = hnew;
      first = target(hnew, _g);
      second = source(h2, _g);
      third = source(hnew, _g);
      break;
    default:
      break;
    }
    CGAL::Euler::split_face(first_to_split, second_to_split, _g);
    // CGAL::internal::set_border(border, _g);
    face_vertices.push_back(first);
    face_vertices.push_back(second);
    face_vertices.push_back(third);
    CGAL::Euler::add_face(face_vertices, _g);
    return hnew;
  }
};
} // namespace Filters
} // namespace FEVV
