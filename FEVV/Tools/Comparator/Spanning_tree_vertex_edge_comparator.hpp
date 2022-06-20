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

#include "FEVV/Tools/Comparator/Vertex_comparators.hpp"
#include "FEVV/Operators/Geometry/triangles.hpp"

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>

#include <limits>
#include <iterator>
#include <map>
#include <set>

#include <stack>

namespace FEVV { 
namespace Comparator 
{
  template<typename Graph>
  std::list< typename boost::graph_traits<Graph>::vertex_descriptor> 
	  get_adjacent_vertices(typename boost::graph_traits<Graph>::vertex_descriptor v, 
                            const Graph& g)
  {
	  typedef boost::graph_traits<Graph>                        GraphTraits;
	  typedef typename GraphTraits::vertex_descriptor           vertex_descriptor;

	  std::list<vertex_descriptor> adjacent_vertices;
	  auto vertices_range = CGAL::vertices_around_target(v, g); // works with two-manifold data structures
	  for (auto v : vertices_range) {
          adjacent_vertices.push_back(v);
	  }

	  return adjacent_vertices;
  }

  template<typename Graph>
  std::list< typename boost::graph_traits<Graph>::vertex_descriptor >
    get_not_processed_adjacent_vertices(typename boost::graph_traits<Graph>::vertex_descriptor v,
            const Graph& g,
            std::map<typename boost::graph_traits<Graph>::vertex_descriptor, bool>& processed_vertices,
            typename boost::graph_traits<Graph>::edge_descriptor min_e)
  {
    typedef boost::graph_traits<Graph>                        GraphTraits;
    typedef typename GraphTraits::vertex_descriptor           vertex_descriptor;
    typedef typename GraphTraits::halfedge_descriptor        halfedge_descriptor;

    std::list<vertex_descriptor> adjacent_vertices;
    if (min_e != edge(GraphTraits::null_halfedge(), g))
    {
      // find the (second) edge that separates old region and new region
      halfedge_descriptor h = halfedge(min_e, g);
      if (target(h, g) != v)
        h = opposite(h, g);

      halfedge_descriptor h_end = h;
      do
      {
        vertex_descriptor v_adj = source(h, g);
        auto it_res = processed_vertices.find(v_adj);
        if ((it_res == processed_vertices.end()) || !it_res->second)
        { // not yet processed
          if (it_res != processed_vertices.end())
            it_res->second = true;
          else
            processed_vertices[v_adj] = true;
          adjacent_vertices.push_back(source(h, g));
        }

        h = opposite(next(h, g), g);
      } while (!(h == h_end));
    }
    else
    {
      auto vertices_range = CGAL::vertices_around_target(v, g); // works with two-manifold data structures
      for (auto v_adj : vertices_range) {
        auto it_res = processed_vertices.find(v_adj);
        if ((it_res != processed_vertices.end()) && it_res->second)
          continue;
        else {
          if (it_res != processed_vertices.end())
            it_res->second = true;
          else
            processed_vertices[v_adj] = true;
          adjacent_vertices.push_back(v_adj);
        }
      }
    }
	  return adjacent_vertices;
  }

/**
 * \brief   Compute the regularity/priority of a neighbor/adjacent point.
 *          To use to identify pivot points.
 *
 * \tparam GeometryTraits The geometric kernel.
 * \param[in] v_s_pos The first point.
 * \param[in] v_t_pos The second point.
 * \param[in] v_s_t_kept_neighbor_pos The neighbor point.
 * \param gt The geometry trait object.
 * \return The regularity of the triange (v_s_t_kept_neighbor_pos, v_s_pos, v_t_pos).
 */
  template< typename GeometryTraits >
  double neighbor_regularity( const typename GeometryTraits::Point& v_s_pos,
                              const typename GeometryTraits::Point& v_t_pos,
                              const typename GeometryTraits::Point& v_s_t_kept_neighbor_pos,
                              const GeometryTraits& gt)
  {
    double prec = 1e-8;
    double perim = FEVV::Operators::Geometry::triangle_perimeter<GeometryTraits>(v_s_t_kept_neighbor_pos, v_s_pos, v_t_pos, gt);
    double area = 0.;
    if (gt.length(v_s_pos, v_t_pos) < prec)
    {
      area = 0.f;
    }
    else
      area = FEVV::Operators::Geometry::triangle_area<GeometryTraits>(v_s_t_kept_neighbor_pos, v_s_pos, v_t_pos, gt);
    double res = ((perim < prec)? area:area/perim);
    static double min_res = 0., max_res = 1024 * 2;
    size_t nb_bit_quantization = static_cast<size_t>(log2(max_res) + 2);
    if (max_res < res)
    {
      max_res = res;
      std::cout << "max_res = " << max_res << std::endl;
    }
    if (nb_bit_quantization > 32)
      nb_bit_quantization = 32;
    size_t two_power_nb_bit = std::pow(2, nb_bit_quantization);
    double delta = (max_res - min_res) / two_power_nb_bit;

    res = std::floor(res / delta) * delta + 0.5 * delta;

    return res;
  }

  /// stable_sort_adjacent_vertices based on true source and target positions
  template<typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  void
    stable_sort_adjacent_vertices(const Graph& /*g*/, 
                                  const PointMap& pm,
                                  typename boost::graph_traits<Graph>::vertex_descriptor /*v*/,
                                  std::list< typename boost::graph_traits<Graph>::vertex_descriptor>& adjacent2v, /// provides adjacent2v in the canonical order
                                  const typename GeometryTraits::Point& v_s_pos,
                                  const typename GeometryTraits::Point& v_t_pos,
                                  const GeometryTraits& gt)
  {
    typedef boost::graph_traits<Graph>              GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
    if (gt.length(v_s_pos, v_t_pos) > 1e-8)
    { // sorting is done only when v_s_pos and v_t_pos are different from each other
      std::map< vertex_descriptor, double > vertex_priority;

      auto it = adjacent2v.begin(), it_e = adjacent2v.end();
      for (; it != it_e; ++it)
      {
        const auto& pos_n = get(pm, *it);
        vertex_priority[*it] = neighbor_regularity(v_s_pos, v_t_pos, pos_n, gt);
      }
      adjacent2v.sort([&vertex_priority](vertex_descriptor v1, vertex_descriptor v2) {
        return (vertex_priority[v1] > vertex_priority[v2]);
      });
    }
  }
  /// stable_sort_adjacent_vertices based on the distance to the bisector plane cutting 
  /// the segment connecting the two furthest vertices in the one-ring neighborhood of v.
  template<typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  void
    stable_sort_adjacent_vertices(const Graph& g,
      const PointMap& pm,
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      std::list< typename boost::graph_traits<Graph>::vertex_descriptor>& adjacent2v, /// provides adjacent2v in the canonical order
      const GeometryTraits& gt)
  {
    typedef boost::graph_traits<Graph>                GraphTraits;
    typedef typename GraphTraits::vertex_descriptor   vertex_descriptor;
    typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
    // sorting is done only when v_s_pos and v_t_pos are different from each other
    std::map< vertex_descriptor, double > vertex_priority;

    // 1) Search for the two furthest vertices
    vertex_descriptor first = GraphTraits::null_vertex(), second = GraphTraits::null_vertex();
    bool find_a_border = false;
    double d_max = -1.0;
    const auto& pos_v = get(pm, v);
    auto it = adjacent2v.begin(), it_e = adjacent2v.end();
    for (; it != it_e; ++it) // find the first intesting vertex (furthest from v) for the border case
    {
      halfedge_descriptor h = halfedge(*it, v, g).first;
      if (CGAL::is_border(h, g))
        find_a_border = true;
      const auto& pos_n = get(pm, *it);
      double d_tmp = gt.length(pos_v, pos_n);
      if (d_tmp > d_max)
      {
        d_max = d_tmp;
        first = *it;
      }
    }
    
    if (!find_a_border)
    {
      // for interior v search for the furthest 2 vertices
      d_max = -1.0;
      it = adjacent2v.begin();
      for (; it != it_e; ++it) 
      {
        const auto& pos_n1 = get(pm, *it);
        auto it2 = it;
        ++it2;
        for (; it2 != it_e; ++it2)
        {
          const auto& pos_n2 = get(pm, *it2);
          double d_tmp = gt.length(pos_n1, pos_n2);
          if (d_tmp > d_max)
          {
            d_max = d_tmp;
            first = *it;
            second = *it2;
          }
        }
      }
    }
    const auto& pos_first = get(pm, first);

    // 2) Compute the plane bisector a.x+b.y+c.z+d=0
    const auto& pos_second = ((find_a_border)? pos_v :get(pm, second));
    auto normal = gt.sub_p(pos_second, pos_first); 
    auto plane_point = ((find_a_border) ? pos_v:gt.add_pv(pos_first, gt.scalar_mult(normal, 0.5f)));
    normal = gt.normalize(normal); // found (a,b,c)
    //double d = -gt.dot_product(normal, gt.sub_p(plane_point,gt.ORIGIN));

    // 3) Compute the priority and then do a stable sort
    it = adjacent2v.begin();
    for (; it != it_e; ++it)
    {
      const auto& pos_n = get(pm, *it);
      double d_tmp = std::abs( gt.dot_product(normal, gt.sub_p(pos_n, plane_point)) );
      vertex_priority[*it] = -d_tmp; // the shortest distance, the higher the priority
    }
    adjacent2v.sort([&vertex_priority](vertex_descriptor v1, vertex_descriptor v2) {
      return (vertex_priority[v1] > vertex_priority[v2]);
    });

  }

  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  class Spanning_tree_vertex_edge_comparator
  {
  public:
  typedef std::size_t                                         IndexType;
	typedef typename FEVV::Vertex_pmap_traits<Graph, IndexType >::pmap_type
                                                              VertexIndexMapType;    
	typedef typename FEVV::Edge_pmap_traits<Graph, IndexType >::pmap_type
                                                              EdgeIndexMapType;

    typedef boost::graph_traits<Graph>                        GraphTraits;
    typedef typename GraphTraits::vertex_descriptor           vertex_descriptor;
    typedef typename GraphTraits::edge_descriptor             edge_descriptor;
    typedef typename GeometryTraits::Scalar                   Scalar;
    typedef typename GeometryTraits::Point                    Point;
  private:
    const Graph& _g;
    const PointMap& _pm;
    const GeometryTraits _gt;
    VertexIndexMapType _vim_st; // spanning tree vertex indices (different from current mesh vertex indices)
    EdgeIndexMapType _eim_st; // spanning tree edge indices

    std::list<vertex_descriptor> _first_vertex_of_each_component; /// the first vertex of each connected component
    std::list<vertex_descriptor> _st_vertices; /// all vertices in the spanning tree order of traversal
    std::list<edge_descriptor> _st_edges; /// all edges in the spanning tree order of traversal
    std::set<edge_descriptor> _has_been_added; // to speed up requests on added edges to the spanning tree

    FEVV::Comparator::Vertex_comparator<Graph, PointMap, GeometryTraits> _cmp_v; // natural ordering + manage tie-break
                                                                                // with vertex degree, then mean of 
                                                                                // adjacent vertex degree

    // Tie-break management
    bool _vertex_tie_break_detected; // will be a root tie-break for 2-manifold meshes
    bool _min_edge_tie_break_detected;
    std::set<vertex_descriptor> _tie_break_vertices; // use set to make sure each tie-break vertex
                                                     // will be processed only once
    void compute_first_vertex_connected_component(bool tie_break_detection)
    {
      _first_vertex_of_each_component.clear();
      _vertex_tie_break_detected = _min_edge_tie_break_detected = false;
      _tie_break_vertices.clear();
	  /////////////////////////////////////////////////////////////////////////
      auto  vertices_range_pair = vertices(_g);
      std::set<vertex_descriptor> remaining_mesh_vertices(vertices_range_pair.first, vertices_range_pair.second);
	  /////////////////////////////////////////////////////////////////////////
	  // Extraction of connected components 
	  std::list< std::list<vertex_descriptor> > connected_components;
	  if (!remaining_mesh_vertices.empty())
	  {
		  do
		  {
              std::map< vertex_descriptor, bool > already_stacked;
			  std::list<vertex_descriptor> component; // current connected component
			  vertex_descriptor v_current = *(remaining_mesh_vertices.begin());

			  {
				  std::list< vertex_descriptor >  list_v_around = get_adjacent_vertices(v_current, _g); 

				  // check if v_current is an isolated vertex (should not happen)
				  if (list_v_around.empty())
				  {
                      already_stacked[v_current] = true;
					  component.push_back(v_current);
					  remaining_mesh_vertices.erase(v_current);
				  }
				  else
				  {
					  std::stack<vertex_descriptor> component_stack;
					  component_stack.push(v_current);

					  while (!component_stack.empty())
					  {
						  vertex_descriptor v = component_stack.top();
						  component_stack.pop();

                          already_stacked[v] = true;
						  component.push_back(v); // v is added to current connected component

              list_v_around = get_not_processed_adjacent_vertices(v, _g, already_stacked, edge(GraphTraits::null_halfedge(), _g));
						  //list_v_around = get_adjacent_vertices(v, _g);
						  for (auto v_around : list_v_around)
						  {
							  typename std::set<vertex_descriptor>::iterator it = remaining_mesh_vertices.find(v_around);
                if (it != remaining_mesh_vertices.end())  // vertex has not been inserted yet
							  {
                  already_stacked[v_around] = true;
								  component_stack.push(v_around);
							  }
							  else
								  continue;
						  }

						  remaining_mesh_vertices.erase(v); // v is removed from remaining vertices
					  }
				  }
			  }

			  connected_components.push_back(component);
			  component.clear();
		  } while (!remaining_mesh_vertices.empty());
	  }
	  /////////////////////////////////////////////////////////////////////////
	  // for each component, find its vertex with smallest coordinates
	  typename std::list< std::list<vertex_descriptor> >::iterator it_list = connected_components.begin();
	  for (; it_list != connected_components.end(); ++it_list)
	  {
          assert(!it_list->empty());
		  auto iter_v = it_list->begin(), 
               iter_v_e = it_list->end();
		  vertex_descriptor best_root = *iter_v; // Must select the vertex with min coordinates as root to ensure the same MST construction
		                                         // if several root exist (2 or more min coord vertices), write a warning
		  ++iter_v;
		  for (; iter_v != iter_v_e; ++iter_v)
		  {
			  if (_cmp_v(*iter_v, best_root))
			  {
				  best_root = *iter_v;
			  }
		  }
      if (tie_break_detection) {
        iter_v = it_list->begin();
        while (iter_v != iter_v_e)
        {
          if ((*iter_v != best_root) && !_cmp_v(*iter_v, best_root) && !_cmp_v(best_root, *iter_v))
          {
            _vertex_tie_break_detected = true;
            _tie_break_vertices.insert(best_root);

            std::list<vertex_descriptor> root_candidates;
            while (iter_v != iter_v_e)
            {
              if (*iter_v != best_root)
              {
                if (!_cmp_v(*iter_v, best_root) && !_cmp_v(best_root, *iter_v))
                  _tie_break_vertices.insert(*iter_v);
                else
                {
                  auto iter_v2 = iter_v;
                  ++iter_v2;

                  while (iter_v2 != iter_v_e)
                  {
                    //while ((iter_v2 != iter_v_e) && (!_cmp_v(*iter_v2, best_root) && !_cmp_v(best_root, *iter_v2))) ++iter_v2; // not needed

                    if ((iter_v2 != iter_v_e) && (!_cmp_v(*iter_v, *iter_v2) && !_cmp_v(*iter_v2, *iter_v)))
                    {
                      _tie_break_vertices.insert(*iter_v);
                      while (iter_v2 != iter_v_e)
                      {
                        if (!_cmp_v(*iter_v, *iter_v2) && !_cmp_v(*iter_v2, *iter_v))
                          _tie_break_vertices.insert(*iter_v2);
                        ++iter_v2;
                      }
                      break;
                    }
                    ++iter_v2;
                  }
                  if (iter_v2 == iter_v_e)
                    root_candidates.push_back(*iter_v); // no tie-break for this candidate
                }
              }
              ++iter_v;
            }
            root_candidates.sort(_cmp_v);
            if (!root_candidates.empty())
              best_root = *(root_candidates.begin());
            else
            {
              std::cerr << "Spanning_tree_vertex_edge_comparator: compute_first_vertex_connected_component (find the root(s)): tie-break detected!" << std::endl;
            }
            break;
          }
          else
            ++iter_v;
        }
      }
		  /////////////////////////////////////////////////////////////////////////
		  _first_vertex_of_each_component.push_back(best_root);

	  }
	  /////////////////////////////////////////////////////////////////////////
	  // order all the smallest vertices
	  _first_vertex_of_each_component.sort(_cmp_v); // (stable sorting: https://www.cplusplus.com/reference/list/list/sort/)
    }

    void compute_st(bool tie_break_detection)
    {
      compute_first_vertex_connected_component(tie_break_detection);

      auto  vertices_range_pair = vertices(_g);
      auto vi = vertices_range_pair.first;
      auto vi_end = vertices_range_pair.second;
      /////////////////////////////////////////////////////////////////////////
      std::map<vertex_descriptor, bool> processed_vertices;
      for (; vi != vi_end; ++vi)
        processed_vertices.insert(std::make_pair(*vi, false));
      /////////////////////////////////////////////////////////////////////////
	  auto edges_range_pair = edges(_g);
	  auto iter_e = edges_range_pair.first;
    auto dist = std::distance(edges_range_pair.first, edges_range_pair.second);
	  for( ;iter_e != edges_range_pair.second; ++iter_e)
		  put(_eim_st, *iter_e, dist);
	  /////////////////////////////////////////////////////////////////////////
      auto it_root_vertices = _first_vertex_of_each_component.begin(), 
           it_root_vertices_e = _first_vertex_of_each_component.end();
      _st_vertices.clear(); _st_edges.clear(); _has_been_added.clear();
      _st_vertices.push_back(*it_root_vertices);
      processed_vertices.find(*it_root_vertices)->second = true;
      auto current_st_it = _st_vertices.begin();
      IndexType cpt = 0, cpt_e=0;
      bool first = true;
      for (; it_root_vertices != it_root_vertices_e; ++it_root_vertices)
      {
        vertex_descriptor current_v = *it_root_vertices;
        assert(current_v != GraphTraits::null_vertex());
        if (!first)
        {
          _st_vertices.push_back(current_v);
          processed_vertices.find(current_v)->second = true;
        }
        else
          first = false;

        do
        {
          // get one ring adjacent vertices
          // make sure min edge ordering of adjacent vertices is used
          // indeed, min edge index ordering is not sensitive to tie_break in 2-manifold
          std::list<vertex_descriptor> adjacent_vertices =
            get_not_processed_adjacent_vertices(current_v, _g, processed_vertices, 
                                                //this->get_spanning_tree_min_incident_edge(current_v) // can be an issue 
                                                this->get_spanning_tree_min_index_incident_edge(current_v) // avoid tie-break issue for 2-manifold (but not for a root and its adjacent vertices)
                                                );
          if (!adjacent_vertices.empty())
          {
            if (degree(current_v, _g) == adjacent_vertices.size())
            { // case of a new region growing start
              adjacent_vertices.sort(_cmp_v); // can be an issue if 2 adjacent vertices have the same coordinates + 
                                              // same degree + same mean degree of adjacent vertices	
                                              // (stable sorting: https://www.cplusplus.com/reference/list/list/sort/)
              if (tie_break_detection) {
                auto it_v_debug = adjacent_vertices.begin(), it_v_debug_e = adjacent_vertices.end();
                for (; it_v_debug != it_v_debug_e; ++it_v_debug)
                {
                  auto it_v2_debug = it_v_debug;
                  if (it_v2_debug != it_v_debug_e)
                    ++it_v2_debug;
                  else
                    continue;
                  if ((it_v2_debug != it_v_debug_e) && !_cmp_v(*it_v_debug, *it_v2_debug) && !_cmp_v(*it_v2_debug, *it_v_debug))
                  {
                    _min_edge_tie_break_detected = true;
                    _tie_break_vertices.insert(*it_v_debug);
                    for (; it_v2_debug != it_v_debug_e; ++it_v2_debug)
                    {
                      if (!_cmp_v(*it_v_debug, *it_v2_debug) && !_cmp_v(*it_v2_debug, *it_v_debug))
                      {
                        _tie_break_vertices.insert(*it_v2_debug);
                        std::cerr << "Spanning_tree_vertex_edge_comparator: compute_st: new region growing start: tie-break detected!" << std::endl;
                      }
                      else
                        break;
                    }
                  }
                  else
                    break;
                }
              }
            }
            // update spanning tree vertices
            _st_vertices.insert(_st_vertices.end(), adjacent_vertices.begin(), adjacent_vertices.end());
            // update spanning tree edges
            auto it_v = adjacent_vertices.begin(), it_v_e = adjacent_vertices.end();
            //std::cout << "START debug adjacent_vertices" << std::endl;
            for (; it_v != it_v_e; ++it_v)
            {
              //std::cout << "(" << get(_pm, *it_v) << ") ;";
              auto res = edge(current_v, *it_v, _g);
              if (res.second)
              {
                put(_eim_st, res.first, cpt_e++);
                _st_edges.push_back(res.first);
                _has_been_added.insert(res.first);
              }
            }
            //std::cout << "\n END debug adjacent_vertices" << std::endl;
            adjacent_vertices.clear();
          }
          put(_vim_st, *current_st_it, cpt++);
          ++current_st_it;
          if (current_st_it == _st_vertices.end())
          {
            --current_st_it;
            break;
          }

          current_v = *current_st_it;
        } while (true);
      }
    }
  public:
    Spanning_tree_vertex_edge_comparator(const Graph& g, const PointMap& pm, bool tie_break_detection=true) :_g(g),
                                                          _pm(pm), 
                                                          _gt(GeometryTraits(g)), 
                                                          _first_vertex_of_each_component(), 
                                                          _st_vertices(), 
                                                          _st_edges(),
                                                          _has_been_added(),
                                                           _cmp_v(FEVV::Comparator::get_vertex_comparator<Graph, PointMap, GeometryTraits>(g, pm)),
                                                          _vertex_tie_break_detected(false),
                                                          _min_edge_tie_break_detected(false),
                                                          _tie_break_vertices()
                                                           {   _vim_st = FEVV::Vertex_pmap_traits<Graph, IndexType >::create(_g);
														       _eim_st = FEVV::Edge_pmap_traits<Graph, IndexType >::create(_g);
															   compute_st(tie_break_detection); }
    Spanning_tree_vertex_edge_comparator(const Graph& g, const PointMap& pm, bool tie_break_detection, const GeometryTraits& gt) :_g(g),
                                                                                    _pm(pm), 
                                                                                    _gt(gt), 
                                                                                    _first_vertex_of_each_component(),
                                                                                    _st_vertices(),
                                                                                    _st_edges(),
                                                                                    _has_been_added(),
                                                                                    _cmp_v(g, pm, gt),
                                                                                    _vertex_tie_break_detected(false),
                                                                                    _min_edge_tie_break_detected(false),
                                                                                    _tie_break_vertices()
                                                                                    { _vim_st = FEVV::Vertex_pmap_traits<Graph, IndexType >::create(_g); 
                                                                                      _eim_st = FEVV::Edge_pmap_traits<Graph, IndexType >::create(_g);
	                                                                                  compute_st(tie_break_detection); }
    bool operator()(vertex_descriptor v1, vertex_descriptor v2) const
    {
      if (v1 == v2)
        return false;
      return get(_vim_st, v1) < get(_vim_st, v2);
    }

	bool operator()(edge_descriptor e1, edge_descriptor e2) const
	{
		if (e1 == e2)
			return false;

		if(_has_been_added.find(e1) == _has_been_added.end())
			return false;
    if (_has_been_added.find(e2) == _has_been_added.end())
			return false;

		return get(_eim_st, e1) < get(_eim_st, e2);
	}

	
	const std::list<vertex_descriptor>& get_first_vertex_of_each_component() const { return _first_vertex_of_each_component; }
	size_t get_nb_of_connected_components() const { return _first_vertex_of_each_component.size(); }
	const std::list<vertex_descriptor>& get_spanning_tree_vertices() const { return _st_vertices; }
	const std::list<edge_descriptor>& get_spanning_tree_edges() const { return _st_edges; }
  /// returns the min incident edge according to Vertex_comparator
  edge_descriptor get_spanning_tree_min_incident_edge(vertex_descriptor v) const 
  { 
    if (degree(v, _g) == 0)
      return edge(GraphTraits::null_halfedge(), _g); // GraphTraits::null_edge() is not defined for all mesh data structures

    std::list< vertex_descriptor >&&  list_v_around = get_adjacent_vertices(v, _g);
    vertex_descriptor min_v = *std::min_element(list_v_around.begin(), list_v_around.end(), _cmp_v);
    auto pair_e_bool = edge(min_v, v, _g);
    return pair_e_bool.first;
  }
  /// returns the min incident edge according to current spanning tree edge indices
  /// edge not yet associated with an index are not considered.
  /// if no incident edge to v has an index, then get_spanning_tree_min_incident_edge is used instead.
  edge_descriptor get_spanning_tree_min_index_incident_edge(vertex_descriptor v) const
  {
    if (degree(v, _g) == 0)
      return edge(GraphTraits::null_halfedge(), _g); // GraphTraits::null_edge() is not defined for all mesh data structures

    edge_descriptor selected_min_edge = edge(GraphTraits::null_halfedge(), _g);
    std::list< vertex_descriptor >&&  list_v_around = get_adjacent_vertices(v, _g);
    auto it_v = list_v_around.begin(), 
         it_v_e = list_v_around.end();
    for (; it_v != it_v_e; ++it_v)
    {
      auto pair_e_bool = edge(v, *it_v, _g);
      if ( _has_been_added.find( pair_e_bool.first ) != _has_been_added.end() )
      {
        if (selected_min_edge == edge(GraphTraits::null_halfedge(), _g))
          selected_min_edge = pair_e_bool.first;
        else if (get(_eim_st, pair_e_bool.first) < get(_eim_st, selected_min_edge))
          selected_min_edge = pair_e_bool.first;
      }
    }
    if (selected_min_edge == edge(GraphTraits::null_halfedge(), _g))
      return get_spanning_tree_min_incident_edge(v);
    else
      return selected_min_edge;
  }
  const Graph &get_spanning_tree_input_graph() const { return _g; }
  const PointMap &get_spanning_tree_vertex_point_map() const { return _pm; }
  IndexType get_spanning_tree_vertex_index(const vertex_descriptor v) const { return get(_vim_st, v); }
  const VertexIndexMapType& get_spanning_tree_vertex_index_map() const { return _vim_st; }

  bool get_tie_break_detected_after_spanning_tree_computation() const {
    return _vertex_tie_break_detected || _min_edge_tie_break_detected;
  }
  const std::set<vertex_descriptor>& get_tie_break_vertices() const { return _tie_break_vertices; }
  };
  /////////////////////////////////////////////////////////////////////////////
  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  static 
	Spanning_tree_vertex_edge_comparator<Graph, PointMap, GeometryTraits>
    get_spanning_tree_comparator(const Graph& g, const PointMap& pm, bool tie_break_detection=true) { return Spanning_tree_vertex_edge_comparator<Graph, PointMap, GeometryTraits>(g, pm, tie_break_detection); }
} // namespace Comparator
} // namespace FEVV

