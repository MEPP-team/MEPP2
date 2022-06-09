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
#include "FEVV/Filters/CGAL/Progressive_Compression/Metrics/QEM_3D.h"

#include "FEVV/Tools/Comparator/Spanning_tree_vertex_edge_comparator.hpp"

#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Memory_comparator.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Raw_positions.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/Butterfly.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/apply_color.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Predictors/stencil.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/Refinement_info.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Uniform_dequantization.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Helpers/Header_handler.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Geometric_metrics.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <limits> 

namespace FEVV {
namespace Filters {

enum class BATCH_CONDITION { ALL_EDGES = 0, REACH_THRESHOLD };

/** \brief Batch_collapser: 
  * Takes an halfedge graph and collapses its edges. 
  * A Batch_collapser object can simplify as many times a mesh as possible. 
  * No need to create an object for each batch.
  * Input: Original mesh. 
  * Output: a simplified mesh.
  **/
template< typename HalfedgeGraph,
          typename PointMap,
          typename Metric,
          typename EdgeColorMap,
          typename VertexColorMap >
class Batch_collapser
{
public:
  using vertex_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor;
  using halfedge_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor;
  using edge_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::edge_descriptor;
  using face_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::face_descriptor;
  using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;
  using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
  using Geometry = typename FEVV::Geometry_traits< HalfedgeGraph >;
  /**
    * \brief Batch_collapser
	* @param g the mesh to simplify
    * @param pm the point map associated to the mesh to simplify
    * @param metric the address (pointer) of an object used to measure the cost
    *        of each edge
    * @param predictor the address (pointer) of an object used to compute the 
    *        geometric residuals
    * @param vcm a property map containing the color of each vertex (useful for
    *        debugging)
    * @param ecm a property map containing the color of each edge (useful for
    *        debugging)
    * @param batch_stop an Enum to know when to stop a batch (collapses every 
    *        edge we can or stop at a threshold)
   **/
  Batch_collapser(
      HalfedgeGraph &g, /// the mesh to simplify
      PointMap &pm, /// the point map associated to the mesh to simplify
      Metric *metric, 
      Predictor< HalfedgeGraph, PointMap >
          *predictor,
      VertexColorMap &vcm, /// for debugging
      EdgeColorMap &ecm, /// for debugging
      BATCH_CONDITION batch_stop)
      : _g(g), _pm(pm), _metric(metric),
        _vcm(vcm), _ecm(ecm), 
        _gt(Geometry(g)), _predictor(predictor)
  {
    _batch_id = 0;
    _forbidden = 0;
    _link_condition = 0;
    _nb_collapse = 0;
    _number_vertices_original = FEVV::size_of_vertices(_g);
    _number_vertices_last = _number_vertices_original;
    _batch_stop = batch_stop;
  }
  ~Batch_collapser() {}
  
  /// \brief stop() : Depending on the chosen Batch condition, will return
  ///  whether we have to continue collapsing or not
  bool stop()
  {
      switch(_batch_stop)
      {
      case BATCH_CONDITION::ALL_EDGES:
        return _metric->is_queue_empty();
        break;
      case BATCH_CONDITION::REACH_THRESHOLD:
        if(!_metric->is_queue_empty())
        {

          return (_metric->get_weight_top() > _metric->get_threshold());
        }
        else
        {
          return true;
        }
        break;

      default:
        return _metric->is_queue_empty();
        break;
      }

  }

  /// Write mesh vertices in a text file (and the number of incident faces) in
  /// the order of the spanning tree (useful for comparing spanning trees 
  /// between compression and decompression).
  void save_spanning_tree(
      const FEVV::Comparator::Spanning_tree_vertex_edge_comparator< HalfedgeGraph,
                                                          PointMap >
          &spanningtree)
  {
    std::ofstream file_encode;
    file_encode.open("encodingspanning" +
                     std::to_string(_batch_id) + ".txt");
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

  void sort_list_memory(
      const FEVV::Comparator::Spanning_tree_vertex_edge_comparator< HalfedgeGraph,
                                                          PointMap >
          &spanningtree)
  {
    Memory_comparator< HalfedgeGraph, PointMap > mc(spanningtree);
    _list_memory.sort(mc);
  }

   /// Computes the L2 distorsion between the original mesh and the current LoD.
  void compute_distortion_l2(
      FEVV::Filters::Geometric_metrics< HalfedgeGraph, PointMap > &g_metric, /// object containing the original mesh
      Header_handler &header, /// header containing every compression info (used to dequantize the mesh)
      double diag, /// diagonal of the mesh bounding box (to normalize)
      bool skip /// we can choose to skip distorsion computation (the value will be -1). For example, we can choose to only compute L2 distorsion on 1 LOD out of 5 to save time
                             )
  {
    if(!skip)
    {
      std::cout << "computing symmetrical L2 (max L2/Hausdorff and RMSE)" << std::endl;
      HalfedgeGraph current_graph = _g; // needs a copy 
      PointMap current_pm = get(boost::vertex_point, current_graph);
      Uniform_dequantization< HalfedgeGraph, PointMap > dq(
          current_graph,
          current_pm,
          header.get_quantization(),
          header.get_dimension(),
          header.get_init_coord());

      dq.point_dequantization();
      double geometric_error =
          g_metric.compute_symmetric_L2(current_graph, true);
      _RMSE_distortion_per_batch.push_back(geometric_error / diag);
      geometric_error =
          g_metric.compute_symmetric_L2(current_graph, false);
      _hausdorff_distortion_per_batch.push_back(geometric_error / diag);
    }
    else
    {
      _RMSE_distortion_per_batch.push_back(-1);
      _hausdorff_distortion_per_batch.push_back(-1);
    }
  }

  private:
    void setup_prediction(
      FEVV::Comparator::Spanning_tree_vertex_edge_comparator< HalfedgeGraph,
      PointMap >
    /*&spanningtree*/)
    {
      auto it_list = _list_memory.begin();


      for( ; it_list != _list_memory.end(); ++it_list)
      {
        _predictor->compute_residuals((*it_list));
      }
    }

public: 
  /// Collapses all the collapsible edges from the mesh
  /// in Batch_collapser::_g
  /// 1) compute all edge costs for current batch;
  /// 2) greedy collapse of the collapsible edges 
  ///    with lowest cost.
  /// while the batch stopping condition is not met.
  void collapse_batch()
  {
    // count every case: forbidden edges, edges which would make our mesh non 
	  // manifold
    _forbidden = 0;
    _link_condition = 0;
    _metric->compute_error(); // compute all edge costs for current batch

    bool collapse_possible = true;
    _nb_collapse = 0;


    //std::cout << "size queue = " << _metric->get_size_queue_t() << std::endl;

    // Collapse until we meet our stopping condition
    while(!stop())
    {
      collapse_possible = collapse_top_queue(); // Get the edge with lowest cost
                                                // and collapse it (if possible).
      if(collapse_possible)
      {
        _nb_collapse += 1;
      }
    }
    int64_t number_current_vertices = FEVV::size_of_vertices(_g);
#ifdef _DEBUG
    // Compute percentage of removed vertices compared to original mesh
    double percentage_vertices = (1 - (double)number_current_vertices /
                                         (double)_number_vertices_original) * 100;
    std::cout << "Percentage vertices removed from start " << percentage_vertices << std::endl;

    percentage_vertices = (1 - (double)number_current_vertices / (double)_number_vertices_last) * 100;
    std::cout << "Percentage vertices removed from last " << percentage_vertices << std::endl;
#endif
    _number_vertices_last = number_current_vertices;
    ///////////////////////////////////////////////////////////////////////////
    FEVV::Comparator::Spanning_tree_vertex_edge_comparator< HalfedgeGraph,
                                                        PointMap > spanningtree(_g, _pm, true);
    sort_list_memory(spanningtree); /// sort _list_memory according to vertex 
                                    /// st traversal
    ///////////////////////////////////////////////////////////////////////////
    // Store refinement info as bitstreams and residuals
    Refinement_info< HalfedgeGraph,
                    PointMap >
        ref_settings(_g, _list_memory);

    // compute connectivity bitmask (to know which halfedges to expand at the
    // decompression step)
    ref_settings.set_connectivity_topology(spanningtree, 
                                           _pm, 
                                           _list_memory); // set one-ring bitmask
    // Compute Predictions
    setup_prediction(spanningtree);
	
    // Concatenate our memory info to make one array per attribute per batch
    // (one array for each bitmask...)
    ref_settings.set_bitMask(spanningtree); // set vertex to split bitmask
    ref_settings.set_error_prediction(); // set residuals array
    ref_settings.set_reverse_bool();// set reverse bitmask
    ///////////////////////////////////////////////////////////////////////////
    _refinements.push_back(ref_settings);
    ///////////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
    std::cout << "batch nb " << _batch_id << " done" << std::endl;
#endif
    _batch_id++;
    _forbidden_edges.clear();
    _list_memory.clear();
  }

  private:
  /// Returns whether the edge collapse causes a normal flip.
  bool are_halfedge_vertex_positions_change_without_normal_flip(halfedge_descriptor h_collapse,
                           const Point &pos_vkept)
  {
    return (is_vertex_position_change_without_normal_flip(h_collapse, source(h_collapse, _g), pos_vkept) &&
            is_vertex_position_change_without_normal_flip(h_collapse, target(h_collapse, _g), pos_vkept));
  }

  /// Returns whether the change of position for the vertex
  /// vertex with pos_vkept will create an incident face normal flip.
  /// \return True if everything is all right, else false.
  bool is_vertex_position_change_without_normal_flip(halfedge_descriptor h_collapse,
                         vertex_descriptor vertex,
                         const Point& pos_vkept)
  {
    assert(vertex == source(h_collapse, _g) || 
	       vertex == target(h_collapse, _g));

    const Point& pvertex = get(_pm, vertex);

    if(pos_vkept == pvertex)
      return true;
    else
    {
      auto iterator_range = CGAL::halfedges_around_target(vertex, _g);
      for(auto h : iterator_range)
      {
        // case where h is a border
        if(CGAL::is_border(h, _g))
          continue;

        face_descriptor f = face(h, _g); // we know that f!=null_face()

        // when f is one of the two incident faces to the collapsed edge
        if(f == face(h_collapse, _g))
          continue;

        if(f == face(opposite(h_collapse, _g), _g))
          continue;

        const Point& p1 = get(_pm, source(h, _g));
        const Point& p0 = get(_pm, source(prev(h, _g), _g));

        if(p0 == pvertex || p1 == pvertex)
        {
          std::cerr << "WARNING: check_faces_normal_around_vertex: we should "
                       "not have dupplicate at this stage ";
          // std::cerr << "merge two vertices that are in the same
          // position?"
          // << std::endl; // I have added "?" because no merge is done
          // here avoid to write things on the console, it slow down the
          // processing
          // std::cerr << "pvertex : " << pvertex << std::endl;
          // std::cerr << "p0 : " << p0 << std::endl;
          // std::cerr << "p1 : " << p1 << std::endl;
        }

        Vector face_normal, face_normal_after_collapse;
        // get the normal of the current face
        // check if the current face is not null
        bool b = FEVV::Math::Vector::are_collinear< Geometry >(
            _gt.sub_p(p0, p1), _gt.sub_p(pvertex, p1));
        if(!b)
          face_normal =
              _gt.normal(p0, p1, pvertex); // gt.unit_normal(p0, p1, pvertex);
        else                               // when the current face is null
        {
          face_normal =
              Vector(std::numeric_limits< typename FEVV::Geometry_traits<
                         HalfedgeGraph >::Scalar >::quiet_NaN(),
                     std::numeric_limits< typename FEVV::Geometry_traits<
                         HalfedgeGraph >::Scalar >::quiet_NaN(),
                     std::numeric_limits< typename FEVV::Geometry_traits<
                         HalfedgeGraph >::Scalar >::quiet_NaN());

          // return false;
        }
        // simulation of the new face after collapse and check if its
        // normal has the same sign 
        auto cross = _gt.cross_product(p0 - p1, p0 - pos_vkept);
        if(cross != Vector(0, 0, 0))
        {

          face_normal_after_collapse = _gt.unit_normal(
              p0, p1, pos_vkept);

          const typename Geometry::Scalar sign =
              _gt.dot_product(face_normal, face_normal_after_collapse);
          if(sign < 0.0)
            return false;
        }
        else
        {
          return false;
        }
      }

      return true;
    }
  }
 
  bool is_forbidden(halfedge_descriptor h, Point /*pvkept*/)
  {
    // forbid border cases for the first version
    return (
        ((_forbidden_edges.find(h) == _forbidden_edges.end() &&
          _forbidden_edges.find(opposite(h, _g)) == _forbidden_edges.end()) &&
         _metric->is_present_in_map(edge(h, _g))));
  }

  /// Used for debugging: instead of collapsing an edge, we color it.
  bool color_top_queue()
  {
    bool mesh_changed = false;
    if(!_metric->is_queue_empty())
    {
      auto current_edge = _metric->getTopQueue();
      // is the  current edge forbidden? 
      auto current_halfedge = halfedge(std::get< 0 >(current_edge), _g);

      if(is_forbidden(current_halfedge, std::get< 2 >(current_edge)))
      {
        if(CGAL::Euler::does_satisfy_link_condition(std::get< 0 >(current_edge),
                                                    _g))
        {
          std::vector< halfedge_descriptor > edges_to_color;
          vertex_descriptor vs = source(current_halfedge, _g);
          vertex_descriptor vt = target(current_halfedge, _g);
          const Point& pvkept = std::get< 2 >(current_edge); // position vkept

          std::vector< Vector > v =
              _predictor->ComputeResiduals(vt, vs, pvkept);
          //_metric->remove_old_edges(vs, vt);
          Collapse_info< HalfedgeGraph, PointMap > memory(_g, _pm);
          memory.set_source_and_target(vs, vt);
          memory.record_v3_v4(current_halfedge);
          edge_descriptor v4target = std::get< 0 >(edge(memory._v4, vt, _g));

          typename boost::property_traits< VertexColorMap >::value_type yellow(
              1, 1, 0);
          ColorEdge(_g, _pm, _ecm, _vcm, v4target, yellow);
          typename boost::property_traits< VertexColorMap >::value_type green(
              0, 1, 0);

          memory.record_pos_removed(current_halfedge);
          // fill memory with vkept init position
          // Point pvkept = std::get< 2 >(current_edge);
          memory.record_pos_vkept(vt);

          memory.record_error_prediction(v);
          mesh_changed = true;
          // std::cout << "collapse done" << std::endl;
          typename boost::property_traits< VertexColorMap >::value_type purple(
              1, 0, 1);
          edge_descriptor v3target = std::get< 0 >(edge(memory._v3, vt, _g));
          ColorEdge(_g, _pm, _ecm, _vcm, v3target, green);


          FEVV::Filters::compute_simple_stencil< HalfedgeGraph,
                                               vertex_descriptor,
                                               halfedge_descriptor >(
              _g, vt, _forbidden_edges, edges_to_color);
          FEVV::Filters::compute_simple_stencil< HalfedgeGraph,
                                               vertex_descriptor,
                                               halfedge_descriptor >(
              _g, vs, _forbidden_edges, edges_to_color);


          memory.record_vkept(vt);

          _list_memory.push_back(std::move(memory));
          // put(_pm, vkept, pvkept);
          typename boost::property_traits< VertexColorMap >::value_type red(
              1, 0, 0);
          typename boost::property_traits< VertexColorMap >::value_type blue(
              0, 0, 1);
          ColorCenter(_g, _pm, _vcm, vt, blue);
          ColorCenter(_g, _pm, _vcm, vs, red);

          FEVV::PMapsContainer pmaps_bag;
          put_property_map(FEVV::vertex_color, _g, pmaps_bag, _vcm);
          FEVV::Filters::write_mesh(
              "tests3dmodelsfirstcollapse.off", _g, pmaps_bag);
          //_metric->UpdateCosts( vkept);
          //_metric->add_new_edges_to_cost_edges(vt);
        }
      }
      else
      {
        // std::cout << "forbidden"  <<std::endl;
      }
      _metric->popQueue(); // dans tous les cas on pop
    }
    return mesh_changed;
  }

  /**
    * @brief collapse_top_queue
    * collapses the edge with the lowest cost, and records the necessary
    * information in a Collapse_info object.  
   **/
  bool collapse_top_queue()
  {
    bool mesh_changed = false;
    if(!_metric->is_queue_empty())
    {
      // we obtain the edge with the lowest cost
      std::tuple< edge_descriptor, double, Point > current_edge =
          _metric->get_top_queue();
		  
      // is the  current edge forbidden?
      auto current_halfedge = halfedge(std::get< 0 >(current_edge), _g);

      if(is_forbidden(current_halfedge, std::get< 2 >(current_edge)))
      {
        if(CGAL::Euler::does_satisfy_link_condition(std::get< 0 >(current_edge),
                                                    _g))
        {
          std::vector< halfedge_descriptor > edges_to_color;
          vertex_descriptor vs = source(current_halfedge, _g);
          vertex_descriptor vt = target(current_halfedge, _g);

          const Point& pvkept = std::get< 2 >(current_edge); // position vkept

          // Create a collapse info objet to store the info we need to
          // reconstruct the edge
          Collapse_info< HalfedgeGraph, PointMap > memory(_g, _pm);
          memory.record_v3_v4(current_halfedge);
          memory.record_v1_v2(get(_pm, vt), get(_pm, vs));

          // give collapse info a unique ID (useful for debugging)
          memory.set_num_collapse(_nb_collapse);

          memory.record_pos_removed(current_halfedge);

          // Check if collapsing this edge would flip normals, thus adding
          // artifacts. If it does, do not collapse.
          bool normal_flip = are_halfedge_vertex_positions_change_without_normal_flip(current_halfedge, pvkept); // true when ok

          if(normal_flip)
          {
              // forbid the simplification of the neighbourhood of the collapsed
              // edge for the current batch
              forbid_edges(_g, vs, _forbidden_edges, edges_to_color);
              forbid_edges(_g, vt, _forbidden_edges, edges_to_color);
              // collapse edge
              vertex_descriptor vkept =
                  CGAL::Euler::collapse_edge(std::get< 0 >(current_edge), _g);
              mesh_changed = true;

              if(memory.get_v3() !=
                 boost::graph_traits< HalfedgeGraph >::null_vertex())
              {
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4834)
#endif
                //halfedge_descriptor h1 =
                    std::get< 0 >(halfedge(memory.get_v3(), vkept, _g));
#if defined _MSC_VER
#pragma warning(pop)
#endif
              }

              if(memory.get_v4() !=
                 boost::graph_traits< HalfedgeGraph >::null_vertex())
              {
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4834)
#endif
                //halfedge_descriptor h2 =
                    std::get< 0 >(halfedge(memory.get_v4(), vkept, _g));
#if defined _MSC_VER
#pragma warning(pop)
#endif
              }

              // set position of the kept vertex
              put(_pm, vkept, pvkept);

              memory.record_pos_vkept(vkept);

              memory.record_vkept(vkept);
              _list_memory.push_back(std::move(memory));
          }
          else
          {
            //  std::cerr << "collapse_top_queue: case not managed yet." << std::endl; /// \todo manage missing cases
          }
          _link_condition++;
        }
      }
      else
      {
        _forbidden++;
        mesh_changed = false;
      }
      _metric->pop_queue(); // dans tous les cas on pop
    }
    return mesh_changed;
  }

public:
  const std::vector< Refinement_info< HalfedgeGraph,
                               PointMap > >&
  get_refinements() const
  {

    return _refinements;
  }

  const std::vector< double >& get_RMSE_distortion() const { return _RMSE_distortion_per_batch; }
  const std::vector< double >& get_hausdorff_distortion() const { return _hausdorff_distortion_per_batch; }
private:
  HalfedgeGraph &_g; /// the mesh being simplified
  PointMap &_pm; 
  Metric *_metric; /// Measuring collapse cost of each edge
  VertexColorMap &_vcm;
  EdgeColorMap &_ecm;
  std::set< halfedge_descriptor > _forbidden_edges; // edges that cannot be collapsed in the current batch
  const Geometry _gt;
  int _batch_id; // current batch number
  int _nb_collapse;
  std::list< Collapse_info< HalfedgeGraph, PointMap > > _list_memory;
  FEVV::Filters::
      Predictor< HalfedgeGraph, PointMap >
          *_predictor;
  std::vector< Refinement_info< HalfedgeGraph,
                               PointMap > >
      _refinements;


  std::vector< double > _RMSE_distortion_per_batch, _hausdorff_distortion_per_batch;
  std::vector< std::vector< bool > > _list_edge_bitmask;
  int _link_condition;
  int _forbidden;
  int64_t _number_vertices_original, _number_vertices_last;
  BATCH_CONDITION _batch_stop;
};

} // namespace Filters
} // namespace FEVV
