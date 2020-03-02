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
#include "FEVV/Wrappings/Graph_traits_aif.h"
#include "FEVV/Wrappings/Geometry_traits_aif.h"
#include "FEVV/Wrappings/Graph_properties_aif.h"
#include "FEVV/Wrappings/properties_aif.h"

#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshReader.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshWriter.hpp"

#include "FEVV/Operators/AIF/similarity.hpp"
#include "FEVV/Operators/AIF/collapse_edge.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshHelpers.h"

#include "FEVV/Operators/Generic/calculate_face_normal.hpp"

#include "FEVV/Tools/IO/FileUtilities.hpp"

#include <vector>
#include <set>
#include <map>

#include <algorithm>
#include <cstdlib> // EXIT_SUCCESS and EXIT_FAILURE symbols

using namespace FEVV::DataStructures::AIF;
template< typename IterFaceType, typename PointMap >
AIFMeshT extract_vertex_local_neighborhood(IterFaceType begin, IterFaceType end, PointMap pm)
{
  typedef typename boost::graph_traits< AIFMeshT >::vertex_descriptor vertex_descriptor;
  AIFMeshT res;
  auto pm_new = get(boost::vertex_point, res);
  std::map<vertex_descriptor, vertex_descriptor> from_current_to_new;
  auto iter_f = begin;
  for (; iter_f != end; ++iter_f)
  {
    std::vector<vertex_descriptor> face_vertices; // for current new face to add
    auto vertex_range = AIFHelpers::incident_vertices(*iter_f);
    auto iter_v = vertex_range.begin();
    for (; iter_v != vertex_range.end(); ++iter_v)
    {
      if (from_current_to_new.find(*iter_v) == from_current_to_new.end())
      {
        vertex_descriptor v_n = AIFHelpers::add_vertex(res);
        from_current_to_new[*iter_v] = v_n;
        auto p = get(pm, *iter_v); // to avoid some writing bugs
        put(pm_new, v_n, p);
      }
      face_vertices.push_back(from_current_to_new[*iter_v]);
    }

    CGAL::Euler::add_face(face_vertices, res);
  }
  return res;
}
template<typename PointMap>
AIFMeshT extract_vertex_local_neighborhood(
  typename boost::graph_traits< AIFMeshT >::vertex_descriptor v,
  const AIFMeshT& g,
  PointMap pm)
{
  auto face_range = AIFHelpers::incident_faces(v);

  return extract_vertex_local_neighborhood(face_range.begin(), face_range.end(), pm);
}
template<typename PointMap>
AIFMeshT extract_edge_local_neighborhood(
  typename boost::graph_traits< AIFMeshT >::edge_descriptor e,
  const AIFMeshT& g,
  PointMap pm)
{
  typedef typename boost::graph_traits< AIFMeshT >::face_descriptor face_descriptor;

  auto face_range_s = AIFHelpers::incident_faces(source(e,g));
  auto face_range_t = AIFHelpers::incident_faces(target(e, g));

  std::set<face_descriptor> faces(face_range_s.begin(), face_range_s.end());
  faces.insert(face_range_t.begin(), face_range_t.end());

  return extract_vertex_local_neighborhood(faces.begin(), faces.end(), pm);
}

static bool argument_analysis(std::string& arg, const std::string& arg_name, bool update_arg_tolower_case = true)
{
	if (arg != "y" && arg != "n" &&
		arg != "Y" && arg != "N")
	{
		std::cerr << "Cannot understand input for " + arg_name + ". "
			"Please use either \"y\" or \"n\" "
			<< std::endl;
		return false;
	}
	else if (update_arg_tolower_case)
	{
		if (arg == "Y")
			arg = "y";
		else if (arg == "N")
			arg = "n";
	}
	return true;
}

static void replace(std::string& to_modify, const std::string& substring_to_search, const std::string& substring_to_use)
{
	size_t pos = to_modify.find(substring_to_search);
	while (pos != std::string::npos)
	{
		to_modify.replace(pos, substring_to_search.size(), substring_to_use);
		pos = to_modify.find(substring_to_search, pos + substring_to_use.size());
	}
}
template< typename MutableFaceIncidentGraph >
typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor
create_new_edge_with_its_incident_vertices(MutableFaceIncidentGraph &g)
{
	typedef boost::graph_traits< MutableFaceIncidentGraph > GraphTraits;

	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::edge_descriptor edge_descriptor;

	edge_descriptor e_tmp = AIFHelpers::add_edge(g);
	vertex_descriptor s_tmp = AIFHelpers::add_vertex(g), t_tmp = AIFHelpers::add_vertex(g);
	AIFHelpers::link_vertex_and_edge(s_tmp, e_tmp, AIFHelpers::vertex_pos::FIRST);
	AIFHelpers::link_vertex_and_edge(t_tmp, e_tmp, AIFHelpers::vertex_pos::SECOND);

	return e_tmp;
}

template< typename MutableFaceIncidentGraph >
static void remove_adjacent_edges(MutableFaceIncidentGraph &g,
	std::vector< typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor >& selected_edges_to_make_independent,
	std::vector< typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor>& removed_edges)
{
	typedef boost::graph_traits< MutableFaceIncidentGraph > GraphTraits;

	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::edge_descriptor edge_descriptor;
	typedef typename GraphTraits::face_descriptor face_descriptor;

	removed_edges.clear();
	std::set< edge_descriptor > set_removed_edges;
	auto iter_e = selected_edges_to_make_independent.begin(), iter_e_e = selected_edges_to_make_independent.end();
	while (iter_e != iter_e_e)
	{
		vertex_descriptor s = source(*iter_e, g);
		vertex_descriptor t = target(*iter_e, g);

		auto edges_pair = in_edges(s, g);
		auto iter_e2 = edges_pair.first;
		for (; iter_e2 != edges_pair.second; ++iter_e2)
		{
			if (*iter_e2 == *iter_e)
				continue;

			//auto iter_res = std::find(selected_edges_to_make_independent.begin(), selected_edges_to_make_independent.end(), *iter_e2);
			//if (iter_res != selected_edges_to_make_independent.end())
			//{
			//	removed_edges.push_back(*iter_res);
			//	selected_edges_to_make_independent.erase(iter_res);
			//}

			set_removed_edges.insert(*iter_e2);
		}

		edges_pair = in_edges(t, g);
		iter_e2 = edges_pair.first;
		for (; iter_e2 != edges_pair.second; ++iter_e2)
		{
			if (*iter_e2 == *iter_e)
				continue;

			//auto iter_res = std::find(selected_edges_to_make_independent.begin(), selected_edges_to_make_independent.end(), *iter_e2);
			//if (iter_res != selected_edges_to_make_independent.end())
			//{
			//	removed_edges.push_back(*iter_res);
			//	selected_edges_to_make_independent.erase(iter_res);
			//}

			set_removed_edges.insert(*iter_e2);
		}

		++iter_e;
	}
	auto iter_e_s = set_removed_edges.begin();
	auto iter_e_s_e = set_removed_edges.end();
	for (; iter_e_s != iter_e_s_e; ++iter_e_s)
	{
		auto iter_res = std::find(selected_edges_to_make_independent.begin(), selected_edges_to_make_independent.end(), *iter_e_s);
		if (iter_res != selected_edges_to_make_independent.end())
		{
			removed_edges.push_back(*iter_res);
			selected_edges_to_make_independent.erase(iter_res);
		}
	}
}

template< typename MutableFaceIncidentGraph >
static size_t nb_different_incident_face_segment_indices(const MutableFaceIncidentGraph &g,
	typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor e
)
{
	std::set<int> segment_indices;
	auto src_incident_faces_pair = in_edges(e, g);
	auto iter_f = src_incident_faces_pair.first;
	for (; iter_f != src_incident_faces_pair.second; ++iter_f)
	{
		int current_face_id = g.template GetProperty< AIFMeshT::face_type::ptr, int >(
			"f:2_manifold_component_seg", (*iter_f)->GetIndex());

		segment_indices.insert(current_face_id);
	}
	return segment_indices.size();
}
template< typename MutableFaceIncidentGraph >
static bool has_that_incident_face_segment_index(const MutableFaceIncidentGraph &g,
	typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor e,
	int index
)
{
	auto src_incident_faces_pair = in_edges(e, g);
	auto iter_f = src_incident_faces_pair.first;
	for (; iter_f != src_incident_faces_pair.second; ++iter_f)
	{
		int current_face_id = g.template GetProperty< AIFMeshT::face_type::ptr, int >(
			"f:2_manifold_component_seg", (*iter_f)->GetIndex());

		if (current_face_id == index)
			return true;
	}
	return false;
}
template< typename MutableFaceIncidentGraph >
static bool has_only_that_incident_face_segment_index(const MutableFaceIncidentGraph &g,
	typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor e,
	int index
)
{
	auto src_incident_faces_pair = in_edges(e, g);
	auto iter_f = src_incident_faces_pair.first;
	for (; iter_f != src_incident_faces_pair.second; ++iter_f)
	{
		int current_face_id = g.template GetProperty< AIFMeshT::face_type::ptr, int >(
			"f:2_manifold_component_seg", (*iter_f)->GetIndex());

		if (current_face_id != index)
			return false;
	}
	return true;
}
template< typename MutableFaceIncidentGraph >
static void replace_vertex_in_incident_edges(MutableFaceIncidentGraph &g,
 typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor v, // current vertex to replace
 typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor replace, // vertex used for the replace
 typename boost::graph_traits<MutableFaceIncidentGraph >::face_descriptor current_f,
	const std::vector< typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor >&
	complex_edges // to forbid the update of edges that are incident to faces than are incident to complex edges (including current complex edge)
	)
{
	typedef boost::graph_traits< MutableFaceIncidentGraph > GraphTraits;

	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::edge_descriptor edge_descriptor;
	typedef typename GraphTraits::face_descriptor face_descriptor;

	assert(v!= AIFHelpers::null_vertex());
	assert(replace != AIFHelpers::null_vertex());
	assert(current_f != AIFHelpers::null_face());

	int current_face_id = g.template GetProperty< AIFMeshT::face_type::ptr, int >(
		"f:2_manifold_component_seg", current_f->GetIndex());

	std::set<vertex_descriptor> forbidden_vertices;
	std::set<edge_descriptor> forbidden_edges;
	auto iter_ec = complex_edges.begin(), iter_ec_e = complex_edges.end();
	forbidden_edges.insert(iter_ec, iter_ec_e);

	// some incident edges of v vertex are put onto replace vertex
	auto src_incident_edges_range = AIFHelpers::incident_edges(v);
	std::set<edge_descriptor> copy_incident_edges_for_removal(src_incident_edges_range.begin(), src_incident_edges_range.end());
	auto iter_e = copy_incident_edges_for_removal.begin();
	auto iter_e_e = copy_incident_edges_for_removal.end();
	for (; iter_e != iter_e_e; ++iter_e)
	{
		if (forbidden_edges.find(*iter_e)!= forbidden_edges.end())
			continue;
		// here we know that *iter_e is not a complexe edge
		// we thus update only edges that have one or two incident faces 
		// with segment id = current_face_id
		if (has_only_that_incident_face_segment_index(g, *iter_e, current_face_id))
		{ 
			if (AIFHelpers::are_incident(v, *iter_e))
			{
				AIFHelpers::remove_edge(v, *iter_e);
				if(AIFHelpers::are_incident(*iter_e, v))
					AIFHelpers::add_vertex(*iter_e, replace, AIFHelpers::vertex_position(*iter_e, v));
				if (!AIFHelpers::are_incident(replace, *iter_e))
					AIFHelpers::add_edge(replace, *iter_e);
			}
		}
		else if (has_that_incident_face_segment_index(g, *iter_e, current_face_id))
		{
			if (!AIFHelpers::are_incident(v, *iter_e))
				continue;

			auto res_pair = AIFHelpers::add_edge(	AIFHelpers::opposite_vertex(*iter_e, v),
				                                    replace, // the new edge will be added in its incident edges
													g);
			assert(res_pair.second); // the edge has been created (false when the edge was existing)
			edge_descriptor new_e = res_pair.first;

			auto face_range_pair = in_edges(*iter_e, g);
			std::set<face_descriptor> copy_incident_faces_for_removal(face_range_pair.first, face_range_pair.second);
			auto iter_f = copy_incident_faces_for_removal.begin(), iter_f_e = copy_incident_faces_for_removal.end();
			for (; iter_f != iter_f_e; ++iter_f)
			{
				if (g.template GetProperty< AIFMeshT::face_type::ptr, int >(
					"f:2_manifold_component_seg", (*iter_f)->GetIndex()) == current_face_id)
				{
					assert(AIFHelpers::are_incident(*iter_f, *iter_e));
					edge_descriptor pe = AIFHelpers::get_edge_of_face_before_edge(*iter_f, *iter_e);
					AIFHelpers::remove_edge(*iter_f, *iter_e);
					AIFHelpers::remove_face(*iter_e, *iter_f);
					AIFHelpers::add_edge_to_face_after_edge(*iter_f, pe, new_e);
					AIFHelpers::add_face(new_e, *iter_f);
				}
			}                                              
		}
	}
}
template< typename MutableFaceIncidentGraph >
static void calculate_previous_and_after_vertices(
	MutableFaceIncidentGraph &g, 
	typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor e,
	typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor pe,
	typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor ae,
	typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor& v_pe_old,
	typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor&v_ae_old,
	std::map<typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor,
	typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor>& v_old_to_v_new, // mapping v_old->v_new
	std::map<typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor,
	typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor>& v_new_to_old    // mapping v_new->v_old
)
{
	v_pe_old = AIFHelpers::common_vertex(pe, e);
	if (v_pe_old == AIFHelpers::null_vertex())
	{
		if (v_new_to_old.find(pe->get_first_vertex()) != v_new_to_old.end())
			v_pe_old = v_new_to_old[pe->get_first_vertex()];
		else if (v_new_to_old.find(pe->get_second_vertex()) != v_new_to_old.end())
			v_pe_old = v_new_to_old[pe->get_second_vertex()];
	}
	else if (v_old_to_v_new.find(v_pe_old) == v_old_to_v_new.end())
	{
		if (v_new_to_old.find(v_pe_old) != v_new_to_old.end())
			v_pe_old = v_new_to_old[v_pe_old];
	}
	if (v_pe_old == AIFHelpers::null_vertex())
	{
		std::cerr << "you may have a similar face with different face id" << std::endl;
	}
	v_ae_old = AIFHelpers::common_vertex(ae, e);
	if (v_ae_old == AIFHelpers::null_vertex())
	{
		if (v_new_to_old.find(ae->get_first_vertex()) != v_new_to_old.end())
			v_ae_old = v_new_to_old[ae->get_first_vertex()];
		else if (v_new_to_old.find(ae->get_second_vertex()) != v_new_to_old.end())
			v_ae_old = v_new_to_old[ae->get_second_vertex()];
	}
	else if (v_old_to_v_new.find(v_ae_old) == v_old_to_v_new.end())
	{
		if (v_new_to_old.find(v_ae_old) != v_new_to_old.end())
			v_ae_old = v_new_to_old[v_ae_old];
	}
	if (v_ae_old == AIFHelpers::null_vertex())
	{
		std::cerr << "you may have a similar face with different face id" << std::endl;
	}
}
template< typename MutableFaceIncidentGraph, typename PointMap >
static void replace_edge_by_new_one_and_update_incidency(
	MutableFaceIncidentGraph &g,
	PointMap& pm,
	typename boost::graph_traits<MutableFaceIncidentGraph >::face_descriptor f, // current face in which e is replaced
	typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor e, // current edge to replace
	typename boost::graph_traits<MutableFaceIncidentGraph >::edge_descriptor replace, // edge used to replace e
	std::map<typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor, 
	         typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor>& v_old_to_v_new, // mapping v_old->v_new
	std::map<typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor, 
	         typename boost::graph_traits<MutableFaceIncidentGraph >::vertex_descriptor>& v_new_to_old    // mapping v_new->v_old
	 )
{
	typedef boost::graph_traits< MutableFaceIncidentGraph > GraphTraits;

	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::edge_descriptor edge_descriptor;
	typedef typename GraphTraits::face_descriptor face_descriptor;
	///////////////////////////////////////////////////////////////////////////
	if (!AIFHelpers::are_incident(f, e))
		return;

	assert(!AIFHelpers::are_incident(f, replace));
	///////////////////////////////////////////////////////////////////////////
	int current_face_id = g.template GetProperty< AIFMeshT::face_type::ptr, int >(
		"f:2_manifold_component_seg", f->GetIndex());
	///////////////////////////////////////////////////////////////////////////
	// for other components: update incidence relationships
	edge_descriptor pe = AIFHelpers::get_edge_of_face_before_edge(f, e);
	edge_descriptor ae = AIFHelpers::get_edge_of_face_after_edge(f, e);

	vertex_descriptor v_pe_old, v_ae_old; 
	calculate_previous_and_after_vertices(g, e, pe, ae, v_pe_old, v_ae_old, v_old_to_v_new, v_new_to_old);
	AIFHelpers::remove_edge(f, e);
	AIFHelpers::remove_face(e, f);
	AIFHelpers::add_edge_to_face_after_edge(f, pe, replace);
	AIFHelpers::add_face(replace, f);

	AIFHelpers::remove_edge(v_pe_old, e);
	AIFHelpers::remove_edge(v_ae_old, e);
	AIFHelpers::add_edge(v_pe_old, replace);
	AIFHelpers::add_edge(v_ae_old, replace);

	if (degree(replace, g) == 1)
	{ // first time replace is added to a face
		auto p = get(pm, v_pe_old); // to avoid some writing bugs
		put(pm, source(replace, g), p);
		p = get(pm, v_ae_old);
		put(pm, target(replace, g), p);

		AIFHelpers::add_edge(source(replace, g), pe);
		if (has_only_that_incident_face_segment_index(g, pe, current_face_id))
		{
			AIFHelpers::add_vertex(pe, source(replace, g), AIFHelpers::vertex_position(pe, v_pe_old));
		}
		// else the update will be done in replace_vertex_in_incident_edges
		AIFHelpers::add_edge(target(replace, g), ae);
		if (has_only_that_incident_face_segment_index(g, ae, current_face_id))
		{
			AIFHelpers::add_vertex(ae, target(replace, g), AIFHelpers::vertex_position(ae, v_ae_old));
		}
		// else the update will be done in replace_vertex_in_incident_edges
		v_old_to_v_new[v_pe_old] = source(replace, g);
		v_old_to_v_new[v_ae_old] = target(replace, g);
		v_new_to_old[source(replace, g)] = v_pe_old;
		v_new_to_old[target(replace, g)] = v_ae_old;
	}
	//else
	//{ 
	//	AIFHelpers::add_edge(v_old_to_v_new[v_pe_old], pe);
	//	if (has_only_that_incident_face_segment_index(g, pe, current_face_id))
	//	{
	//		if((v_old_to_v_new.find(v_pe_old) != v_old_to_v_new.end()) && !AIFHelpers::are_incident(pe, v_old_to_v_new[v_pe_old]))
	//			AIFHelpers::add_vertex(pe, v_old_to_v_new[v_pe_old], AIFHelpers::vertex_position(pe, v_pe_old));
	//	}
	//	// else the update will be done in replace_vertex_in_incident_edges
	//	AIFHelpers::add_edge(v_old_to_v_new[v_ae_old], ae);
	//	if (has_only_that_incident_face_segment_index(g, ae, current_face_id))
	//	{
	//		if ((v_old_to_v_new.find(v_ae_old) != v_old_to_v_new.end()) && !AIFHelpers::are_incident(ae, v_old_to_v_new[v_ae_old]))
	//			AIFHelpers::add_vertex(ae, v_old_to_v_new[v_ae_old], AIFHelpers::vertex_position(ae, v_ae_old));
	//	}
	//	// else the update will be done in replace_vertex_in_incident_edges
	//}
}

template< typename MutableFaceIncidentGraph, typename PointIndexMap >
static bool has_different_vertex_indices(MutableFaceIncidentGraph &g, const PointIndexMap& idm,
	typename boost::graph_traits<MutableFaceIncidentGraph >::face_descriptor f)
{
	typedef typename boost::graph_traits<MutableFaceIncidentGraph>::vertex_descriptor vertex_descriptor;
	auto vertex_range = AIFHelpers::incident_vertices(f);
	auto iter_v = vertex_range.begin();
	std::map<vertex_descriptor, int> v_to_index;
	for (; iter_v != vertex_range.end(); ++iter_v)
	{
		v_to_index[*iter_v] = get(idm, *iter_v);
	}

	return (degree(f, g) == v_to_index.size()) ;
}

static int process_one_meshfile(	const std::string& input_file_path, 
									const std::string& colorize_mesh, 
									const std::string& remove_isolated_elements, // except isolated faces for the time being
									const std::string& resolve_vertices_with_similar_incident_edges,
									const std::string& make_2_mani_not_2_mani)
{
	typedef AIFMeshReader reader_type;
	typedef AIFMeshWriter writer_type;
	typedef boost::graph_traits< reader_type::output_type > GraphTraits;
	typedef FEVV::Geometry_traits<AIFMeshT> GeometryTraits;

	typedef typename GraphTraits::vertex_iterator vertex_iterator;
	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::edge_iterator edge_iterator;
	typedef typename GraphTraits::edge_descriptor edge_descriptor;
	typedef typename GraphTraits::face_iterator face_iterator;
	typedef typename GraphTraits::face_descriptor face_descriptor;
	/////////////////////////////////////////////////////////////////////////////
	reader_type my_reader;
	reader_type::ptr_output ptr_mesh;
	try
	{
		ptr_mesh = my_reader.read(input_file_path);
	}
	catch (const std::invalid_argument &e)
	{
		std::cerr << "Invalid Argument Exception catch while reading "
			<< input_file_path << " :" << std::endl
			<< e.what() << std::endl;
		exit(EXIT_FAILURE);
	}
	catch (const std::ios_base::failure &e)
	{
		std::cerr << "IOS Failure Exception catch while reading " << input_file_path
			<< " :" << std::endl
			<< e.what() << std::endl;
		exit(EXIT_FAILURE);
	}
	catch (const std::length_error &le)
	{
		std::cerr << "[AIF] Exception caught while reading input file "
			<< input_file_path << ": " << le.what() << std::endl;
		BOOST_ASSERT_MSG(false, "[AIF] Exception caught while reading input file.");
	}
	/////////////////////////////////////////////////////////////////////////////
	GeometryTraits gt(*ptr_mesh);
	auto pos_pm = get(boost::vertex_point, *ptr_mesh);
	auto vertex_idm = get(boost::vertex_index, *ptr_mesh);
	auto edge_idm = get(boost::edge_index, *ptr_mesh);
	/////////////////////////////////////////////////////////////////////////////
	const AIFMeshT::Vector red(255, 0, 0);
	const AIFMeshT::Vector green(0, 255, 0);
	const AIFMeshT::Vector blue(0, 0, 255);
	const AIFMeshT::Vector white(255, 255, 255);
	if (colorize_mesh == "y")
	{
		if (!ptr_mesh->isPropertyMap< AIFMeshT::vertex_type::ptr >("v:color"))
		{ // create a property map to store vertex colors if not already created
			ptr_mesh->AddPropertyMap< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >("v:color");
			if (!ptr_mesh->isPropertyMap< AIFMeshT::vertex_type::ptr >("v:color"))
				throw std::runtime_error("Failed to create vertex-color property map.");
		}
		if (!ptr_mesh->isPropertyMap< AIFMeshT::edge_type::ptr >("e:color"))
		{ // create a property map to store edge colors if not already created
			ptr_mesh->AddPropertyMap< AIFMeshT::edge_type::ptr, AIFMeshT::Vector >("e:color");
			if (!ptr_mesh->isPropertyMap< AIFMeshT::edge_type::ptr >("e:color"))
				throw std::runtime_error("Failed to create edge-color property map.");
		}
		if (!ptr_mesh->isPropertyMap< AIFMeshT::face_type::ptr >("f:color"))
		{ // create a property map to store face colors if not already created
			ptr_mesh->AddPropertyMap< AIFMeshT::face_type::ptr, AIFMeshT::Vector >("f:color");
			if (!ptr_mesh->isPropertyMap< AIFMeshT::face_type::ptr >("f:color"))
				throw std::runtime_error("Failed to create face-color property map.");
		}
	}
	if (make_2_mani_not_2_mani == "y")
	{
		if (!ptr_mesh->isPropertyMap< AIFMeshT::face_type::ptr >("f:2_manifold_component_seg"))
		{ // create a property map to store face colors if not already created
			ptr_mesh->AddPropertyMap< AIFMeshT::face_type::ptr, int >("f:2_manifold_component_seg");
			if (!ptr_mesh->isPropertyMap< AIFMeshT::face_type::ptr >("f:2_manifold_component_seg"))
				throw std::runtime_error("Failed to create face 2-manifold components property map.");
		}
		// compute the face segmentation
		std::set<face_descriptor> set_current_component,
			all_faces(faces(*ptr_mesh).first, faces(*ptr_mesh).second); // not yet processed
		bool first_phase = true // true when it still remains at least one face with no incident complex edge
			, second_phase_first_subphase=true;
		int current_id = 0, nb_used_id = 0;
		auto iter_f_first = std::find_if(all_faces.begin(), all_faces.end(), AIFHelpers::has_no_incident_complex_edge);
		if (iter_f_first == all_faces.end())
		{
			first_phase = false;
			iter_f_first = all_faces.begin();
		}
		set_current_component.insert(*iter_f_first);
		all_faces.erase(iter_f_first);

		while (!set_current_component.empty())
		{
			face_descriptor current_f = *set_current_component.begin();
			ptr_mesh->SetProperty< AIFMeshT::face_type::ptr, int >(
				"f:2_manifold_component_seg", current_f->GetIndex(), current_id);

			set_current_component.erase(set_current_component.begin());

     		auto vector_faces = AIFHelpers::adjacent_faces(current_f, true); 
			auto iter_f = vector_faces.begin(), iter_f_end = vector_faces.end();
			for (; iter_f != iter_f_end; ++iter_f)
			{
				assert(*iter_f != current_f);
				if ((all_faces.find(*iter_f) != all_faces.end()) // not already processed
					&& (!first_phase || AIFHelpers::has_no_incident_complex_edge(*iter_f))
					&& AIFHelpers::have_consistent_orientation(current_f, *iter_f)
					)
				{
					edge_descriptor e_tmp = AIFHelpers::common_edge(current_f, *iter_f);
					if ((!first_phase || AIFHelpers::is_surface_regular_edge(e_tmp))
						//&& (AIFHelpers::is_regular_vertex(source(e_tmp, *ptr_mesh)) 
						//&&  AIFHelpers::is_regular_vertex(target(e_tmp, *ptr_mesh)))
						)
					{
						set_current_component.insert(*iter_f);
						all_faces.erase(*iter_f);
					}
				}
			}
			if (!first_phase)
			{
				if (second_phase_first_subphase)
				{
					if ((set_current_component.size() > 1))
					{ // second phase (segmentation of faces incident to at least one complex edge)
						std::set<face_descriptor> new_face_set;
						// only keep the face with closest normal to current_f (except current_f itself)
						auto iter_f = set_current_component.begin(), iter_f_end = set_current_component.end();
						auto current_n = FEVV::Operators::calculate_face_normal(current_f, *ptr_mesh, pos_pm, gt);
						bool first = true;
						double best_val = -1.1;
						face_descriptor best_f, other_f1 = AIFHelpers::null_face(), other_f2 = AIFHelpers::null_face();
						for (; iter_f != iter_f_end; ++iter_f)
						{
							if (*iter_f == current_f)
								continue;

							auto neighbor_n = FEVV::Operators::calculate_face_normal(*iter_f, *ptr_mesh, pos_pm, gt);
							double tmp = gt.dot_product(current_n, neighbor_n);
							if (first)
							{
								best_val = tmp;
								best_f = *iter_f;
								first = false;
								other_f1 = *iter_f;
							}
							else
							{
								other_f2 = *iter_f;
								if (tmp > best_val)
								{
									best_val = tmp;
									best_f = *iter_f;
								}
							}
						}

						if (set_current_component.size() == 2)
						{
							current_n = FEVV::Operators::calculate_face_normal(other_f1, *ptr_mesh, pos_pm, gt);
							auto neighbor_n = FEVV::Operators::calculate_face_normal(other_f2, *ptr_mesh, pos_pm, gt);
							double tmp = gt.dot_product(current_n, neighbor_n);
							if (tmp < best_val)
								new_face_set.insert(best_f);
						}
						else
							new_face_set.insert(best_f);
						set_current_component.swap(new_face_set);
					}
					///////////////////////////////////////////////////////////////
					second_phase_first_subphase = false;
					if (set_current_component.size() > 0)
					{
						++nb_used_id;
						current_id = nb_used_id;
					}
				}
				else
				{
					std::set<face_descriptor> new_face_set;
					set_current_component.swap(new_face_set);
				}
			}
			if (set_current_component.empty() && !all_faces.empty())
			{
				iter_f_first = std::find_if(all_faces.begin(), all_faces.end(), AIFHelpers::has_no_incident_complex_edge);
				if (iter_f_first == all_faces.end())
				{ // we estimate the correct next face id by looking at its already processed neighbors
					first_phase = false;
					iter_f_first = all_faces.begin();
					auto vector_faces = AIFHelpers::adjacent_faces(*iter_f_first, true);
					auto iter_f = vector_faces.begin(), iter_f_end = vector_faces.end();
					auto current_n = FEVV::Operators::calculate_face_normal(*iter_f_first, *ptr_mesh, pos_pm, gt);
					int id_best_match = -1; // look for the face with a consistent orientation
					                        // that has the closest normal correcpondance
					double best_val = -1.1;
					for (; iter_f != iter_f_end; ++iter_f)
					{
						if ((all_faces.find(*iter_f) == all_faces.end()) // already processed
							//&& AIFHelpers::has_no_incident_complex_edge(*iter_f)
							&& AIFHelpers::have_consistent_orientation(*iter_f_first, *iter_f)
							)
						{
							int current_f_id = ptr_mesh->template GetProperty< AIFMeshT::face_type::ptr, int >(
								"f:2_manifold_component_seg", (*iter_f)->GetIndex());

							auto neighbor_n = FEVV::Operators::calculate_face_normal(*iter_f, *ptr_mesh, pos_pm, gt);
							double tmp = gt.dot_product(current_n, neighbor_n);
							if (id_best_match == -1)
							{
								id_best_match = current_f_id;
								best_val = tmp;
							}
							else if (FEVV::Operators::are_similar_faces(*iter_f_first, *iter_f, *ptr_mesh) // similar in topology 
								     || (std::abs(tmp+1.)<1e-5) // similar in spirit, not in topology
								    )
							{ // they must be matched together
								id_best_match = current_f_id;
								best_val = tmp;
								break;
							}
							else
							{
								if (tmp > best_val)
								{
									best_val = tmp;
									id_best_match = current_f_id;
								}
							}
						}
					}
					if (id_best_match == -1)
					{ // no valid neighbord found => current_id must take the next available integer value
						++nb_used_id;
     					current_id = nb_used_id ;
					}
					else
					{ // at least one valid neighbor => take the same id
						current_id = id_best_match;
					}

					second_phase_first_subphase = true;
				}
				else
				{
					++current_id;
					++nb_used_id;
				}
				set_current_component.insert(*iter_f_first);
				all_faces.erase(iter_f_first);
				
			}
		}
		std::cout << "segmented into " << (nb_used_id + 1) << " 2-manifold components." << std::endl;
	}
	/////////////////////////////////////////////////////////////////////////////
	// 1 - Count number of isolated/complex/normal elements

	// VERTICES
	auto iterator_pair_v = vertices(*ptr_mesh);
	vertex_iterator vi = iterator_pair_v.first;
	vertex_iterator vi_end = iterator_pair_v.second;
	int nb_isolated_vertices = 0, nb_cut_vertices = 0, nb_t_junction_vertices = 0,
		nb_vertices_with_similar_incident_edges = 0;
	std::vector< vertex_descriptor > v_to_remeove;
	std::vector< vertex_descriptor > cut_vertices;
	std::vector< vertex_descriptor > t_junction_vertices;
	for (; vi != vi_end; ++vi)
	{
		if (AIFHelpers::is_isolated_vertex(*vi))
		{
			++nb_isolated_vertices;

			if (remove_isolated_elements == "y")
				v_to_remeove.push_back(*vi);
			else
			{
				if (colorize_mesh == "y")
					ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
						"v:color", (*vi)->GetIndex(), red);
			}
		}
		else if (FEVV::DataStructures::AIF::AIFMeshHelpers::is_a_T_junction_vertex(*vi, *ptr_mesh, gt))
		{
			++nb_t_junction_vertices;
			//std::cout << "T-junction detected!" << std::endl;
			if (make_2_mani_not_2_mani == "y")
				t_junction_vertices.push_back(*vi);
			else
				if (colorize_mesh == "y")
					ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
						"v:color", (*vi)->GetIndex(), blue);
		}
		else if (AIFHelpers::is_cut_vertex(*vi) &&
			!FEVV::DataStructures::AIF::AIFMeshHelpers::has_adjacent_T_junction_vertex(*vi, *ptr_mesh, gt) && // we must process T-junction vectices first
			!AIFHelpers::has_dangling_or_complex_incident_edge(*vi) // we must process dangling and complex edges first 
			)                                                       // and we can even either create new cut vertices 
		{                                                           // during the complex edge processing or resolve a cut
			++nb_cut_vertices;                                      // vertex...
			if (make_2_mani_not_2_mani == "y")
				cut_vertices.push_back(*vi);
			else
				if (colorize_mesh == "y")
					ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
						"v:color", (*vi)->GetIndex(), blue);
		}
		else if (FEVV::Operators::contains_similar_incident_edges(
			*vi, *ptr_mesh))
		{

			++nb_vertices_with_similar_incident_edges;

			if (resolve_vertices_with_similar_incident_edges == "y")
				FEVV::Operators::resolve_similar_incident_edges(*vi, *ptr_mesh);
			else
				if (colorize_mesh == "y")
					ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
						"v:color", (*vi)->GetIndex(), blue);
		}
		else
		{
			if (colorize_mesh == "y")
				ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
					"v:color", (*vi)->GetIndex(), green);
		}
	}

	// EDGES
	auto iterator_pair_e = edges(*ptr_mesh);
	edge_iterator ei = iterator_pair_e.first;
	edge_iterator ei_end = iterator_pair_e.second;
	int nb_isolated_edges = 0, nb_dangling_edges = 0, nb_complex_edges = 0,
		nb_border_edges = 0;
	std::vector< edge_descriptor > e_to_remeove;
	std::vector< edge_descriptor > dangling_edges, complex_edges;
	for (; ei != ei_end; ++ei)
	{
		if (AIFHelpers::is_isolated_edge(*ei))
		{
			++nb_isolated_edges;
			if (remove_isolated_elements == "y")
				e_to_remeove.push_back(*ei);
			else
			{
				if (colorize_mesh == "y")
					ptr_mesh->SetProperty< AIFMeshT::edge_type::ptr, AIFMeshT::Vector >(
						"e:color", (*ei)->GetIndex(), red);
			}
		}
		else if (AIFHelpers::is_dangling_edge(*ei))
		{
			++nb_dangling_edges;
			if (make_2_mani_not_2_mani == "y")
				dangling_edges.push_back(*ei);
			else
			{
				if (colorize_mesh == "y")
					ptr_mesh->SetProperty< AIFMeshT::edge_type::ptr, AIFMeshT::Vector >(
						"e:color", (*ei)->GetIndex(), red);
			}
		}
		else if (AIFHelpers::is_complex_edge(*ei))
		{
			++nb_complex_edges;
			if (make_2_mani_not_2_mani == "y")
				complex_edges.push_back(*ei);
			else
			{
				if (colorize_mesh == "y")
				{
					/* ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
					"v:color", (*ei)->get_first_vertex()->GetIndex(), blue);
					ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
					"v:color", (*ei)->get_second_vertex()->GetIndex(), blue);*/

					ptr_mesh->SetProperty< AIFMeshT::edge_type::ptr, AIFMeshT::Vector >(
						"e:color", (*ei)->GetIndex(), blue);
				}
			}
		}
		//else if (AIFHelpers::has_inconsistent_incident_face_orientation(*ei))
		//{
		//  // should be treated by resolving complex edges
		//}
		else if (AIFHelpers::is_surface_regular_edge(*ei))
		{
			++nb_border_edges;
			if (colorize_mesh == "y")
				ptr_mesh->SetProperty< AIFMeshT::edge_type::ptr, AIFMeshT::Vector >(
					"e:color", (*ei)->GetIndex(), white);
		}
		else
		{
			if (colorize_mesh == "y")
				ptr_mesh->SetProperty< AIFMeshT::edge_type::ptr, AIFMeshT::Vector >(
					"e:color", (*ei)->GetIndex(), green);
		}
	}

	// FACES
	auto iterator_pair_f = faces(*ptr_mesh);
	face_iterator fi = iterator_pair_f.first;
	face_iterator fi_end = iterator_pair_f.second;
	int nb_isolated_faces = 0, nb_degenerated_faces = 0;
	std::vector< face_descriptor > f_to_remeove;
	for (; fi != fi_end; ++fi)
	{
		if (AIFHelpers::is_isolated_face(*fi))
		{ // note that removing isolated faces may create hole in the surface!
		  // especially after processing cut vertices.
		  // Therefore, it is rather advised to preserved them, and use
		  // a higher level algorithm to keep or remove them.
			++nb_isolated_faces;
			//if (remove_isolated_elements == "y")
			//	f_to_remeove.push_back(*fi);
			//else
			{
				if (colorize_mesh == "y")
					ptr_mesh->SetProperty< AIFMeshT::face_type::ptr, AIFMeshT::Vector >(
						"f:color", (*fi)->GetIndex(), red);
			}
		}
		else
		{
			if (colorize_mesh == "y")
				ptr_mesh->SetProperty< AIFMeshT::face_type::ptr, AIFMeshT::Vector >(
					"f:color", (*fi)->GetIndex(), green);
		}

		if (AIFHelpers::is_degenerated_face(*fi))
			++nb_degenerated_faces;
	}
	/////////////////////////////////////////////////////////////////////////////
	// 2 - Remove isolated elements + output filename generation
	std::string suffix = "";
	if (remove_isolated_elements == "y")
	{
		if (nb_isolated_vertices != 0 
			|| nb_isolated_edges != 0 
			//|| nb_isolated_faces != 0
			)
		{
			AIFHelpers::remove_faces(
				f_to_remeove.begin(), f_to_remeove.end(), ptr_mesh);
			AIFHelpers::remove_edges(
				e_to_remeove.begin(), e_to_remeove.end(), ptr_mesh);
			AIFHelpers::remove_vertices(
				v_to_remeove.begin(), v_to_remeove.end(), ptr_mesh);
			suffix = "_without_isolated_elmt";
			if (resolve_vertices_with_similar_incident_edges == "y")
				suffix += "_and_similar_incident_edges";
		}
		else if (resolve_vertices_with_similar_incident_edges == "y")
			suffix = "_without_similar_incident_edges";

		if (make_2_mani_not_2_mani == "y")
			suffix += "_without_non_manifold_elm";
	}
	else if (resolve_vertices_with_similar_incident_edges == "y")
	{
		suffix = "_without_similar_incident_edges";
		if (make_2_mani_not_2_mani == "y")
			suffix += "_without_non_manifold_elm";
	}
	else if (make_2_mani_not_2_mani == "y")
		suffix = "_without_non_manifold_elm";
	/////////////////////////////////////////////////////////////////////////////
	// 3 - Process non-2-manifold elements

	// remove dangling edges (+ update list of cut vertices)
	auto iter_e = dangling_edges.begin(), iter_e_end = dangling_edges.end();
	while (iter_e != iter_e_end)
	{
		if (degree(source(*iter_e, *ptr_mesh), *ptr_mesh) == 1)
		{
			auto v_tmp = target(*iter_e, *ptr_mesh);
			FEVV::Operators::collapse_edge_keep_target(*ptr_mesh, *iter_e);
			if (AIFHelpers::is_cut_vertex(v_tmp))
			{
				++nb_cut_vertices;                                            // vertex...
				if (make_2_mani_not_2_mani == "y")
					cut_vertices.push_back(v_tmp);
				else
					if (colorize_mesh == "y")
						ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
							"v:color", v_tmp->GetIndex(), blue);
			}
		}
		else
		{
			auto v_tmp = source(*iter_e, *ptr_mesh);
			FEVV::Operators::collapse_edge_keep_source(*ptr_mesh, *iter_e);
			if (AIFHelpers::is_cut_vertex(v_tmp))
			{
				++nb_cut_vertices;                                            // vertex...
				if (make_2_mani_not_2_mani == "y")
					cut_vertices.push_back(v_tmp);
				else
					if (colorize_mesh == "y")
						ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
							"v:color", v_tmp->GetIndex(), blue);
			}
		}
		++iter_e;
	}

	// decomplexify complex edges (+ update list of cut vertices)
	std::vector< edge_descriptor > selected_complex_edges(complex_edges.begin(), complex_edges.end()), 
		                           remaining_complex_edges;
	remove_adjacent_edges(*ptr_mesh, selected_complex_edges, remaining_complex_edges);
	while (!selected_complex_edges.empty())
	{
		std::cout << "New decomplexification step for an independent set of complex edges" << std::endl; 
		iter_e = selected_complex_edges.begin(), iter_e_end = selected_complex_edges.end();
		int i = 0;
		while (iter_e != iter_e_end)
		{
			{ // mesh file writing debug
				auto res_mesh =
					extract_edge_local_neighborhood(*iter_e, *ptr_mesh, pos_pm);
				writer_type my_writer;
				try
				{
					my_writer.write(res_mesh,
						FEVV::FileUtils::get_file_name(input_file_path) + "_edge_" + FEVV::StrUtils::convert((*iter_e)->GetIndex()) + ".off");
				}
				catch (...)
				{
					std::cout << "writing failed";
					exit(EXIT_FAILURE);
				}
			}
			// complex edges
			std::cout << "Information of current complex edge (" << i << "/" << selected_complex_edges.size() << "): " << std::endl; // debug
			std::cout << "\tedge id: " << get(edge_idm, *iter_e) << std::endl; // debug
			std::cout << "\tsource vertex: " << get(vertex_idm, source(*iter_e, *ptr_mesh)) << std::endl; // debug
			std::cout << "\ttarget vertex: " << get(vertex_idm, target(*iter_e, *ptr_mesh)) << std::endl; // debug
			std::cout << "\tdegree: " << degree(*iter_e, *ptr_mesh) << std::endl; // debug
			std::cout << "\tnb incident different segment: " << nb_different_incident_face_segment_indices(*ptr_mesh, *iter_e) << std::endl; // debug
			vertex_descriptor v_pe_old = AIFHelpers::null_vertex();
			vertex_descriptor v_ae_old = AIFHelpers::null_vertex();
			auto faces_range_pair = in_edges(*iter_e, *ptr_mesh); // get incident faces
			std::set<face_descriptor> current_faces(faces_range_pair.first, faces_range_pair.second),
				next_faces;
			std::map<int, edge_descriptor> face_id_to_edge;
			std::map<vertex_descriptor, vertex_descriptor> v_old_to_v_new, v_new_to_old;
			bool first = true;
			while (!current_faces.empty())
			{
				face_descriptor current_f = *current_faces.begin();
				int current_id = ptr_mesh->template GetProperty< AIFMeshT::face_type::ptr, int >(
					"f:2_manifold_component_seg", current_f->GetIndex());
				if (first)
				{ // the first component keep the initial complex edge
					face_id_to_edge.insert(std::make_pair(current_id, *iter_e));
				}
				else if (face_id_to_edge.find(current_id) == face_id_to_edge.end())
				{ // other components use a duplicate
					std::cout << "create a new edge" << std::endl; // debug
					edge_descriptor e_tmp = create_new_edge_with_its_incident_vertices(*ptr_mesh);
					face_id_to_edge.insert(std::make_pair(current_id, e_tmp));
				}

				auto iter_f = current_faces.begin(), iter_f_end = current_faces.end();
				for (; iter_f != iter_f_end; ++iter_f)
				{
					if (ptr_mesh->template GetProperty< AIFMeshT::face_type::ptr, int >(
						"f:2_manifold_component_seg", (*iter_f)->GetIndex()) != current_id)
						next_faces.insert(*iter_f);
					else
					{ // during the current step we process all faces with the same current_id
						if (first)
							continue;

						// for other components: update incidence relationships
						edge_descriptor pe = AIFHelpers::get_edge_of_face_before_edge(*iter_f, *iter_e);
						edge_descriptor ae = AIFHelpers::get_edge_of_face_after_edge(*iter_f, *iter_e);
						calculate_previous_and_after_vertices(*ptr_mesh, *iter_e, pe, ae, v_pe_old, v_ae_old, v_old_to_v_new, v_new_to_old);

						replace_edge_by_new_one_and_update_incidency(
							*ptr_mesh,
							pos_pm,
							*iter_f,
							*iter_e,
							face_id_to_edge[current_id], // edge used to replace e
							v_old_to_v_new,
							v_new_to_old);

						std::cout << "\torientation: " << get(vertex_idm, v_pe_old) << " -> " << get(vertex_idm, v_ae_old); // debug
						std::cout << "\treplaced by: " << get(vertex_idm, v_old_to_v_new[v_pe_old]) << " -> " << get(vertex_idm, v_old_to_v_new[v_ae_old]) << " with edge id: " << get(edge_idm, face_id_to_edge[current_id]); // debug
						std::cout << std::endl; // debug

						assert(!AIFHelpers::is_degenerated_face(*iter_f));

						bool is_the_last = true;
						auto iter_f2 = iter_f;
						++iter_f2;
						for (; iter_f2 != iter_f_end; ++iter_f2)
							if (ptr_mesh->template GetProperty< AIFMeshT::face_type::ptr, int >(
								"f:2_manifold_component_seg", (*iter_f2)->GetIndex()) == current_id)
							{
								is_the_last = false;
								break;
							}
						if (is_the_last)
						{
							// do this replacement for the last current face with the same id
							replace_vertex_in_incident_edges(*ptr_mesh, v_pe_old, v_old_to_v_new[v_pe_old], *iter_f, complex_edges);
							replace_vertex_in_incident_edges(*ptr_mesh, v_ae_old, v_old_to_v_new[v_ae_old], *iter_f, complex_edges);

							assert(!AIFHelpers::are_incident(v_pe_old, face_id_to_edge[current_id]));
							assert(!AIFHelpers::are_incident(v_ae_old, face_id_to_edge[current_id]));
							assert(has_different_vertex_indices(*ptr_mesh, vertex_idm, *iter_f));
						}
					}
				}
				first = false;
				current_faces.swap(next_faces);
				next_faces.clear();
			}

			// update cut vertices
			auto iter_map = face_id_to_edge.begin(), iter_map_e = face_id_to_edge.end();
			for (; iter_map != iter_map_e; ++iter_map)
			{
				auto v_tmp = source(iter_map->second, *ptr_mesh);
				if (AIFHelpers::is_cut_vertex(v_tmp))
				{
					++nb_cut_vertices;                                            // vertex...
					if (make_2_mani_not_2_mani == "y")
						cut_vertices.push_back(v_tmp);
					else
						if (colorize_mesh == "y")
							ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
								"v:color", v_tmp->GetIndex(), blue);
				}
				v_tmp = target(iter_map->second, *ptr_mesh);
				if (AIFHelpers::is_cut_vertex(v_tmp))
				{
					++nb_cut_vertices;                                            // vertex...
					if (make_2_mani_not_2_mani == "y")
						cut_vertices.push_back(v_tmp);
					else
						if (colorize_mesh == "y")
							ptr_mesh->SetProperty< AIFMeshT::vertex_type::ptr, AIFMeshT::Vector >(
								"v:color", v_tmp->GetIndex(), blue);
				}
			}

			++iter_e;
			++i;
		}
		selected_complex_edges = remaining_complex_edges;
		remove_adjacent_edges(*ptr_mesh, selected_complex_edges, remaining_complex_edges);
	}
	// duplicate cut vertices [local partition and cutting]
	auto iter_v = cut_vertices.begin(), iter_v_end = cut_vertices.end();
	while (iter_v != iter_v_end)
	{
		//{ // mesh file writing debug
		//	auto res_mesh =
		//		extract_vertex_local_neighborhood(*iter_v, *ptr_mesh, pos_pm);
		//	writer_type my_writer;
		//	try
		//	{
		//		my_writer.write(res_mesh,
		//			FEVV::FileUtils::get_file_name(input_file_path) + "_vertex_" + FEVV::StrUtils::convert((*iter_v)->GetIndex()) + ".off");
		//	}
		//	catch (...)
		//	{
		//		std::cout << "writing failed";
		//		exit(EXIT_FAILURE);
		//	}
		//}
		// cut vertices
		if (AIFHelpers::is_cut_vertex(*iter_v))
			std::cout << "degree of current cut vertex: " << degree(*iter_v, *ptr_mesh) << std::endl; // debug
		else
		{
			std::cout << "current vertex is not more a cut vertex!" << std::endl;
			continue;
		}
		// pieces of surface need to replace *iter_v by a new vertex
		// except for the first piece of surface which keeps *iter_v
		auto face_range = AIFHelpers::incident_faces(*iter_v);

		std::set<face_descriptor> current_faces(face_range.begin(), face_range.end()),
			next_faces;
		std::map<int, vertex_descriptor> face_id_to_vertex;
		bool first = true;
		while (!current_faces.empty())
		{
			face_descriptor current_f = *current_faces.begin();
			int current_id = ptr_mesh->template GetProperty< AIFMeshT::face_type::ptr, int >(
				"f:2_manifold_component_seg", current_f->GetIndex());

			if (first)
			{ // the first component keep the initial cut vertex
				face_id_to_vertex.insert(std::make_pair(current_id, *iter_v));

				// also continue for the same piece of surface else add to next_faces
				// remove edge-adjacent face around *iter_v
				edge_descriptor e1 = AIFHelpers::common_edge(*iter_v, current_f);
				edge_descriptor e2 = AIFHelpers::common_edge(*iter_v, current_f, e1);
				face_descriptor current_f_tmp = current_f;
				while ((degree(e1, *ptr_mesh) == 2) &&
					(current_faces.find(current_f_tmp) != current_faces.end()))
				{
					auto face_range2 = AIFHelpers::incident_faces(e1);
					auto it_f2 = face_range2.begin();
					for (; it_f2 != face_range2.end(); ++it_f2)
					{
						if (*it_f2 != current_f_tmp)
						{
							current_faces.erase(current_f_tmp);
							current_f_tmp = *it_f2;
							e1 = AIFHelpers::common_edge(*iter_v, current_f_tmp, e1);
							break;
						}
					}
				}
				if (current_f_tmp != current_f)
					current_faces.erase(current_f_tmp);
				if (degree(e2, *ptr_mesh) == 2)
					current_faces.insert(current_f);
				while ((degree(e2, *ptr_mesh) == 2) &&
					(current_faces.find(current_f) != current_faces.end()))
				{
					auto face_range2 = AIFHelpers::incident_faces(e2);
					auto it_f2 = face_range2.begin();
					for (; it_f2 != face_range2.end(); ++it_f2)
					{
						if (*it_f2 != current_f)
						{
							current_faces.erase(current_f);
							current_f = *it_f2;
							e2 = AIFHelpers::common_edge(*iter_v, current_f, e2);
							break;
						}
					}
				}

			}
			else
			{
				std::cout << "create a new vertex" << std::endl; // debug
				vertex_descriptor v_tmp = AIFHelpers::add_vertex(*ptr_mesh);
				auto p = get(pos_pm, *iter_v); // to avoid some writing bugs
				put(pos_pm, v_tmp, p);

				face_id_to_vertex.insert(std::make_pair(current_id, v_tmp));
			}

			auto iter_f = current_faces.begin(), iter_f_end = current_faces.end();
			std::map<face_descriptor, bool> face_to_next;
			for (; iter_f != iter_f_end; ++iter_f)
				if (current_f == *iter_f)
					face_to_next[*iter_f] = false;
				else if (AIFHelpers::are_incident_to_vertex_and_edge_connected(*iter_v, current_f, *iter_f))
					face_to_next[*iter_f] = false;
				else
					face_to_next[*iter_f] = true;
			iter_f = current_faces.begin();
			for (; iter_f != iter_f_end; ++iter_f)
			{
				if (face_to_next[*iter_f])
					next_faces.insert(*iter_f);
				else
				{
					if (first)
					{
						continue;
					}

					auto edge_range = AIFHelpers::incident_edges(*iter_f);
					auto it_e = edge_range.begin();
					for (; it_e != edge_range.end(); ++it_e)
					{
						if (AIFHelpers::are_incident(*it_e, *iter_v))
						{
							AIFHelpers::remove_edge(*iter_v, *it_e);
							AIFHelpers::add_edge(face_id_to_vertex[current_id], *it_e);
							AIFHelpers::add_vertex(*it_e, face_id_to_vertex[current_id], AIFHelpers::vertex_position(*it_e, *iter_v));
						}
					}
				}
			}

			first = false;
			current_faces.swap(next_faces);
			next_faces.clear();
		}

		if (AIFHelpers::is_cut_vertex(*iter_v))
		{
			std::cout << "failed to correct current vertex!" << std::endl; // debug
			{ // mesh file writing debug
				auto res_mesh =
					extract_vertex_local_neighborhood(*iter_v, *ptr_mesh, pos_pm);
				writer_type my_writer;
				try
				{
					my_writer.write(res_mesh,
						FEVV::FileUtils::get_file_name(input_file_path) + "_vertex_after_" + FEVV::StrUtils::convert((*iter_v)->GetIndex()) + ".off");
				}
				catch (...)
				{
					std::cout << "writing failed";
					exit(EXIT_FAILURE);
				}
			}
		}
		++iter_v;
	}
	/////////////////////////////////////////////////////////////////////////////
	assert(AIFHelpers::check_mesh_validity(ptr_mesh));
	/////////////////////////////////////////////////////////////////////////////
	// 4 - Save corrected mesh
	writer_type my_writer;
	try
	{
		if (FEVV::FileUtils::get_file_extension(input_file_path) == ".ply")
			my_writer.write(ptr_mesh,
				FEVV::FileUtils::get_file_name(input_file_path) + suffix +
				".off");
		else
			my_writer.write(ptr_mesh,
				FEVV::FileUtils::get_file_name(input_file_path) + suffix +
				FEVV::FileUtils::get_file_extension(input_file_path));
	}
	catch (...)
	{
		std::cout << "writing failed";
		exit(EXIT_FAILURE);
	}
	/////////////////////////////////////////////////////////////////////////////
	std::cout << "The mesh file named "
		<< FEVV::FileUtils::get_file_name(input_file_path) +
		FEVV::FileUtils::get_file_extension(input_file_path);
	if (nb_isolated_vertices > 0 || nb_t_junction_vertices > 0 || nb_cut_vertices > 0 || nb_isolated_edges > 0 ||
		nb_dangling_edges > 0 || nb_complex_edges > 0)
		std::cout << " is not 2-manifold " << std::endl;
	else
		std::cout << " is 2-manifold" << std::endl;

	std::string prefix =
		(remove_isolated_elements == "y")
		? "Number of removed isolated"
		: "Number of isolated";
	std::cout << prefix << " vertices: " << nb_isolated_vertices << std::endl;
	std::cout << prefix << " edges: " << nb_isolated_edges << std::endl;
	std::cout << prefix << " faces (not removed yet): " << nb_isolated_faces << std::endl;
	/////////////////////////////////////////////////////////////////////////////
	prefix = (resolve_vertices_with_similar_incident_edges == "y")
		? "Number of resolved vertices with similar incident edges: "
		: "Number of vertices with similar incident edges: ";
	std::cout << prefix << nb_vertices_with_similar_incident_edges << std::endl;
	/////////////////////////////////////////////////////////////////////////////
	std::cout << "Number of T-junction vertices: " << nb_t_junction_vertices << std::endl;
	std::cout << "Number of cut vertices: " << nb_cut_vertices << std::endl;
	std::cout << "Number of dangling edges: " << nb_dangling_edges << std::endl;
	std::cout << "Number of complex edges: " << nb_complex_edges;
	if(remaining_complex_edges.empty())
		std::cout << std::endl;
	else
		std::cout << " (" << remaining_complex_edges.size() << " not resolved due to local dependency). " << std::endl;
	std::cout << "Number of surface border edges: " << nb_border_edges << std::endl;
	std::cout << "Number of degenerated faces: " << nb_degenerated_faces << std::endl;
	/////////////////////////////////////////////////////////////////////////////
	return EXIT_SUCCESS;
}

int
main(int narg, char **argv)
{
  if(narg < 2 || narg > 6)
  {
	  std::cerr << "Cannot proceed arguments. Please use " << argv[0]
		  << " meshfilename_or_meshfolder colorize_mesh [remove_isolated_elements "
      "[resolve_vertices_with_similar_incident_edges [make_2_mani_not_2_mani]]]"
              << std::endl;
    return EXIT_FAILURE;
  }
  std::string input_path(argv[1]);
  if (!boost::filesystem::exists(input_path) && !boost::filesystem::is_directory(input_path))
  {
	  std::cerr << input_path << " does not exist. Exit. " << std::endl;
	  return EXIT_FAILURE;
  }
  std::string colorize_mesh = (argv[2]);
  if(!argument_analysis(colorize_mesh, "colorizing output mesh", true))
	  return EXIT_FAILURE;
  std::string remove_isolated_elements(((narg > 3) ? argv[3] : "n"));
  if (!argument_analysis(remove_isolated_elements, "removing of not isolated element", true))
	  return EXIT_FAILURE;
  std::string resolve_vertices_with_similar_incident_edges(
      ((narg > 4) ? argv[4] : "n"));
  if (!argument_analysis(resolve_vertices_with_similar_incident_edges, "resolving similar/duplicate incident edges", true))
	  return EXIT_FAILURE;
  std::string make_2_mani_not_2_mani(((narg > 5) ? argv[5] : "n"));
  if (!argument_analysis(make_2_mani_not_2_mani, "resolving not 2 manifold elements", true))
	  return EXIT_FAILURE;
  /////////////////////////////////////////////////////////////////////////////
  if(!boost::filesystem::is_directory(input_path))
    process_one_meshfile(input_path,
                         colorize_mesh,
                         remove_isolated_elements,
                         resolve_vertices_with_similar_incident_edges,
                         make_2_mani_not_2_mani);
  else
  {
	  ///////////////////////////////////////////////////////////////////////
	  boost::filesystem::directory_iterator end_iter;
	  for (boost::filesystem::directory_iterator dir_itr(input_path);
                                                dir_itr != end_iter;
                                                          ++dir_itr)
	  {
		  try
		  {
			  if (boost::filesystem::is_directory(dir_itr->status()))
			  {
				  std::string dir_name_without_base = dir_itr->path().string();

				  std::string command = "\"";
				  for (int i = 0; i < narg; ++i)
				  {
					  switch (i)
					  {
					  case 0:
					  {
						  std::string cmd_tmp = argv[i];
						  //replace(cmd_tmp, " ", "\\ ");
						  command += "\"" + cmd_tmp + "\"";
					  }
						  break;
					  case 1:
						  //replace(dir_name_without_base, " ", "\\ ");
						  command += "\"" + dir_name_without_base + "\"";
						  break;
					  default:
						  command += std::string(argv[i]);
					  }

					  if (i < narg - 1)
						  command += std::string(" ");
				  }
				  command += "\"";
				  int r = system(command.c_str());
                  if (r == -1)
                  {
                    std::cerr << " the following command is run: " << command << " failed (child process generation issue)." << std::endl;
                  }

				  continue; // ignore it for the current process
			  }
			  else if (!boost::filesystem::is_regular_file(dir_itr->status()))
			  {
				  continue; // ignore it
			  }

		  }
		  catch (const std::exception & ex)
		  {
			  std::cout << dir_itr->path().filename() << " " << ex.what() << std::endl;
		  }

		  if (dir_itr->path().extension().string() == ".msh" || dir_itr->path().extension().string() == ".off" || dir_itr->path().extension().string() == ".obj"
#ifdef FEVV_USE_VTK
			  || dir_itr->path().extension().string() == ".vtk" || dir_itr->path().extension().string() == ".vtp" || dir_itr->path().extension().string() == ".vtu"
#endif
			  )
		  {
			  process_one_meshfile(dir_itr->path().string(),
                                   colorize_mesh,
                                   remove_isolated_elements,
                                   resolve_vertices_with_similar_incident_edges,
                                   make_2_mani_not_2_mani);
		  }
	  }
  }
  /////////////////////////////////////////////////////////////////////////////
  system("pause");
  return EXIT_SUCCESS;
}
