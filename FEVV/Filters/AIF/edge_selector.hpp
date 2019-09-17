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

#include "FEVV/Wrappings/Graph_traits_extension_aif.h" 
#include "FEVV/Wrappings/Geometry_traits.h"

#include "FEVV/Operators/Generic/Manifold/no_normal_flip_for_collapse.hpp"

#include "FEVV/Operators/AIF/topology_predicates.hpp"

#include <vector>

namespace FEVV {
namespace Filters {

/*!
 * \brief Function used for cleaning the topology of mesh g.
 *        This can be seen as a preprocessing step for some
 *        geometry processing algorithms.
 *
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \tparam GeometryTraits The geometric kernel when available. This is defaulted
 *         to FEVV::Geometry_traits<FaceGraph>. 
 * \param[in] g The FaceGraph instance.
 * \param[in] forbid_non_satisfying_link_condition Boolean to avoid edge 
 *            non-satisfaying the link condition.
 * \param[in] forbid_non_manifold_edge Boolean to avoid non-manifold edge.
 * \param[in] forbid_edge_collapse_creating_non_manifold_split Boolean to 
 *            avoid edge whose collapse creates a non-manifold vertex split.
 * \param[in] forbid_border_edge Boolean to avoid border edge.
 * \param[in] forbid_inner_edge Boolean to avoid inner edge.
 * \param[in] forbid_non_triangular_incident_face_to_edge Boolean to avoid 
 *            edge whose incident faces are not triangular.
 * \param[in] forbid_edges_that_are_adjacent_to_collapsed_edges Boolean to 
 *            avoid edge adjacent to already selected edges.
 * \param[in] forbid_one_edges_that_are_one_ring_edges_of_collapsed_edge_vertices
 *            Boolean to avoid edge that is onto the one-ring of selected
 *            edges' vertices.
 * \param[in] forbid_edges_that_are_incident_to_one_ring_of_collapsed_edge_vertices
              Boolean to avoid edge that is adjacent to edges onto the one-ring
 *            of selected edges' vertices.
 * \param[in] gt The Geometry Traits object.
 * \return  an std::vector of edge_descriptor corresponding to edges
 *          selected from g.
 */
  template<typename FaceGraph,
           typename GeometryTraits = FEVV::Geometry_traits< FaceGraph > >
  std::vector<typename boost::graph_traits<FaceGraph>::edge_descriptor>
    edge_selector_for_collapse(const FaceGraph& g,
                               bool forbid_non_satisfying_link_condition, 
                               bool forbid_non_manifold_edge, 
                               bool forbid_edge_collapse_creating_non_manifold_split,
                               bool forbid_border_edge, // border edge processing is usually different from inner edges
                               bool forbid_inner_edge, // should be false, except when testing border edges
                               bool forbid_non_triangular_incident_face_to_edge, // most of algorithms can only handle edge collapses if all incident faces are triangular
                               bool forbid_edges_that_are_adjacent_to_collapsed_edges, // should be true most of the time, since local encoding will be tricky otherwise
                               bool forbid_one_edges_that_are_one_ring_edges_of_collapsed_edge_vertices, // must be true for full independent edges collapses; you should set the same value for next parameter
                               bool forbid_edges_that_are_incident_to_one_ring_of_collapsed_edge_vertices, // must be true for full independent edges collapses; you should set the same value for previous parameter
                               const GeometryTraits &gt							   
      )
  {
    typedef FEVV::DataStructures::AIF::AIFTopologyHelpers Helpers;
    std::vector<typename boost::graph_traits<FaceGraph>::edge_descriptor> batch_of_edges_to_collapse;
    std::set<typename boost::graph_traits<FaceGraph>::edge_descriptor> forbidden_edges_to_collapse;
    /////////////////////////////////////////////////////////////////////////////
	auto pm = get(boost::vertex_point, g);
    auto edge_range_mesh = Helpers::edges(g);
    for (auto edge_it = edge_range_mesh.begin(); edge_it != edge_range_mesh.end(); ++edge_it)
    {
      if (forbidden_edges_to_collapse.find(*edge_it) == forbidden_edges_to_collapse.end())
      {
        bool is_a_2_manifold_collapse = Helpers::is_2_manifold_edge(*edge_it);
        bool edge_ok_for_collapse = !forbid_non_manifold_edge ||
          ( is_a_2_manifold_collapse && // cannot encode neither cut vertices, nor complex edges with triplets (v,left, right)
            Helpers::is_one_ring_2_manifold(source(*edge_it, g)) &&
            Helpers::is_one_ring_2_manifold(target(*edge_it, g)));

        //if (forbid_non_satisfying_link_condition && 
        //    !CGAL::Euler::does_satisfy_link_condition(*edge_it, g) // not enough for non-manifold meshes and makes g non const
        //  )
        //{
        //  forbidden_edges_to_collapse.insert(*edge_it);
        //  continue;
        //}

        if (forbid_edge_collapse_creating_non_manifold_split && edge_ok_for_collapse)
        {
          // case where the merge of the 2 one-rings makes a cut vertex that cannot be processed by the 2-manifold vertex split
          if (Helpers::is_surface_interior_edge(*edge_it) &&
            Helpers::is_surface_border_vertex(source(*edge_it, g)) &&
            Helpers::is_surface_border_vertex(target(*edge_it, g))
            )
            edge_ok_for_collapse = false;
        }

        if (forbid_border_edge && edge_ok_for_collapse)
        {
          // forbid border edge collapse
          if (Helpers::is_surface_border_edge(*edge_it))
            edge_ok_for_collapse = false;
        }
        if (forbid_inner_edge && edge_ok_for_collapse)
        {
          // forbid inner edge collapse
          if (Helpers::is_surface_interior_edge(*edge_it))
            edge_ok_for_collapse = false;
        }
        bool is_a_triangular_collapse = true;

        if (forbid_non_triangular_incident_face_to_edge && edge_ok_for_collapse )
        {
          // forbid polygonal edge collapse
          if ( !FEVV::Operators::has_only_incident_triangular_faces(*edge_it, g) )
          {
            edge_ok_for_collapse = false;
            is_a_triangular_collapse = false;
          }
        }

        if (forbid_non_triangular_incident_face_to_edge && is_a_2_manifold_collapse && is_a_triangular_collapse && edge_ok_for_collapse)
        {		
          // local check to avoid face flipping 
          if (!FEVV::Operators::no_normal_flip_for_collapse(g, pm, *edge_it, source(*edge_it, g), get(pm, target(*edge_it, g)), gt) ||
            !FEVV::Operators::no_normal_flip_for_collapse(g, pm, *edge_it, target(*edge_it, g), get(pm, source(*edge_it, g)), gt)
            )
          {
            edge_ok_for_collapse = false;
          }  
        }

        if (forbid_non_satisfying_link_condition && edge_ok_for_collapse)
        {
          auto incident_v1_range = Helpers::adjacent_vertices(source(*edge_it, g));
          auto incident_v2_range = Helpers::adjacent_vertices(target(*edge_it, g));
          auto iter_v1 = incident_v1_range.begin();
          for (; (iter_v1 != incident_v1_range.end()) && edge_ok_for_collapse; ++iter_v1)
          {
            if (*iter_v1 == target(*edge_it, g))
              continue;
            auto iter_v2 = incident_v2_range.begin();
            for (; iter_v2 != incident_v2_range.end(); ++iter_v2)
            {
              if (*iter_v2 == source(*edge_it, g))
                continue;
              if (*iter_v1 == *iter_v2)
              {
                typename boost::graph_traits<FaceGraph>::face_descriptor f = Helpers::common_face(Helpers::common_edge(source(*edge_it, g), *iter_v1),
                  Helpers::common_edge(target(*edge_it, g), *iter_v2)
                  );
                // When there is a common neighbor, and if no face exist between these 3 vertices, we have two situations:
                // 1) vertex in "the middle" of a triangle or 2) tunnel that can be sewn by the collapse 
                if (f == Helpers::null_face())
                {
                  edge_ok_for_collapse = false;
                  break;
                }
                else if ((degree(Helpers::common_edge(source(*edge_it, g), *iter_v1), g) == 1) &&
                  (degree(Helpers::common_edge(target(*edge_it, g), *iter_v2), g) == 1)
                  )
                { // case where you create a dangling edge [to avoid for the time being since most of mesh readers do not support it, except AIF native reader]
                  edge_ok_for_collapse = false;
                  break;
                }
              }
            }
          }
        }
        if (!edge_ok_for_collapse)
        {
          forbidden_edges_to_collapse.insert(*edge_it);
          continue;
        }

        batch_of_edges_to_collapse.push_back(*edge_it);
        // remove incident edges from edges to collapse
        auto v1_ad_container = Helpers::adjacent_vertices((*edge_it)->get_first_vertex());
        auto v2_ad_container = Helpers::adjacent_vertices((*edge_it)->get_second_vertex());

        if (forbid_edges_that_are_adjacent_to_collapsed_edges)
        {
          for (auto it = v1_ad_container.begin(); it != v1_ad_container.end(); ++it)
          {
            if (*it == (*edge_it)->get_second_vertex())
              continue;
            forbidden_edges_to_collapse.insert(edge(*it, (*edge_it)->get_first_vertex(), g).first);
          }

          for (auto it = v2_ad_container.begin(); it != v2_ad_container.end(); ++it)
          {
            if (*it == (*edge_it)->get_first_vertex())
              continue;
            forbidden_edges_to_collapse.insert(edge(*it, (*edge_it)->get_second_vertex(), g).first);
          }
        }

        if (forbid_one_edges_that_are_one_ring_edges_of_collapsed_edge_vertices)
        {
          // remove one-ring edges from edges to collapse
          auto edge_set_v1 = Helpers::get_unordered_one_ring_edges((*edge_it)->get_first_vertex());
          auto edge_set_v2 = Helpers::get_unordered_one_ring_edges((*edge_it)->get_second_vertex());

          for (auto it = edge_set_v1.begin(); it != edge_set_v1.end(); ++it)
          {
            forbidden_edges_to_collapse.insert(*it);
          }

          for (auto it = edge_set_v2.begin(); it != edge_set_v2.end(); ++it)
          {
            forbidden_edges_to_collapse.insert(*it);
          }
        }

        if (forbid_edges_that_are_incident_to_one_ring_of_collapsed_edge_vertices)
        {
          // remove incident edges to one-ring vertices
          for (auto it = v1_ad_container.begin(); it != v1_ad_container.end(); ++it)
          {
            auto edge_range = Helpers::incident_edges(*it);

            for (auto it2 = edge_range.begin(); it2 != edge_range.end(); ++it2)
            {
              forbidden_edges_to_collapse.insert(*it2);
            }
          }

          for (auto it = v2_ad_container.begin(); it != v2_ad_container.end(); ++it)
          {
            auto edge_range = Helpers::incident_edges(*it);

            for (auto it2 = edge_range.begin(); it2 != edge_range.end(); ++it2)
            {
              forbidden_edges_to_collapse.insert(*it2);
            }
          }
        }
      }
    }
    return batch_of_edges_to_collapse;
  }
} // namespace Filters
} // namespace FEVV

