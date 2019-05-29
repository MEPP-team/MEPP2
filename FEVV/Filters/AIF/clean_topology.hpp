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

#include <boost/graph/graph_traits.hpp>

#include "FEVV/Filters/AIF/resolve_similar_faces.hpp"
#include "FEVV/Filters/AIF/resolve_similar_edges.hpp"

#include <vector>

namespace FEVV {
namespace Filters {

/*!
 * \brief Function use for cleaning the topology of mesh g.
 *        This can be seen as a preprocessing step for some
 *        geometry processing algorithms.
 *
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param[in,out] g The FaceGraph instance.
 * \param[in] remove_isolated_vertices Boolean to remove isolated
 *            vertices if true. 
 * \param[in] remove_isolated_edges Boolean to remove isolated
 *            edges if true. 
 * \param[in] remove_isolated_faces Boolean to remove isolated
 *            faces if true. 
 * \param[in] remove_similar_edges Boolean to remove similar/parallel
 *            edges if true.  
 * \param[in] remove_similar_faces Boolean to remove similar/parallel
 *            faces if true.   
 * \return  true if a topology change has been done, false othewise.
 */
 
  template<typename FaceGraph>
  bool clean_topology(FaceGraph& g,
                      bool remove_isolated_vertices=true,
                      bool remove_isolated_edges=true,
                      bool remove_isolated_faces=false,
                      bool remove_similar_edges=true,
                      bool remove_similar_faces=true)
  {
    bool ret = false;
	if(remove_isolated_faces)
	{	
      std::vector< typename boost::graph_traits<FaceGraph>::face_descriptor > f_to_remove;
      auto iterator_pair_f = faces(g);
      auto fi = iterator_pair_f.first;
      auto fi_end = iterator_pair_f.second;
      for (; fi != fi_end; ++fi)
      {
        if (FEVV::DataStructures::AIF::AIFTopologyHelpers::is_isolated_face(*fi))
          f_to_remove.push_back(*fi);
      }
      if (!f_to_remove.empty())
        ret = true;
      FEVV::DataStructures::AIF::AIFTopologyHelpers::remove_faces(
        f_to_remove.begin(), f_to_remove.end(), g);
    }
	if(remove_isolated_edges)
	{
      std::vector< typename boost::graph_traits<FaceGraph>::edge_descriptor > e_to_remove;
      auto iterator_pair_e = edges(g);
      auto ei = iterator_pair_e.first;
      auto ei_end = iterator_pair_e.second;
      for (; ei != ei_end; ++ei)
      {
        if (FEVV::DataStructures::AIF::AIFTopologyHelpers::is_isolated_edge(*ei))
          e_to_remove.push_back(*ei);
      }
      if (!e_to_remove.empty())
        ret = true;
      FEVV::DataStructures::AIF::AIFTopologyHelpers::remove_edges(
        e_to_remove.begin(), e_to_remove.end(), g);
    }
	if(remove_isolated_vertices)
	{
      std::vector< typename boost::graph_traits<FaceGraph>::vertex_descriptor > v_to_remove;
      auto iterator_pair_v = vertices(g);
      auto vi = iterator_pair_v.first;
      auto vi_end = iterator_pair_v.second;
      for (; vi != vi_end; ++vi)
      {
        if (FEVV::DataStructures::AIF::AIFTopologyHelpers::is_isolated_vertex(*vi))
          v_to_remove.push_back(*vi);
      }
      std::sort(v_to_remove.begin(),
        v_to_remove.end(),
        [](const typename boost::graph_traits<FaceGraph>::vertex_descriptor & v1, 
           const typename boost::graph_traits<FaceGraph>::vertex_descriptor & v2) -> bool
      {
        return v1->GetIndex() > v2->GetIndex();
      }); // to remove first vertices with higher index value
      if (!v_to_remove.empty())
        ret = true;
      FEVV::DataStructures::AIF::AIFTopologyHelpers::remove_vertices(
        v_to_remove.begin(), v_to_remove.end(), g);
	}
    // remove similar faces and edges
    size_t nb_e = size_of_edges(g), nb_f = size_of_faces(g);
    if(remove_similar_faces)
	  FEVV::Filters::resolve_similar_faces(g, false); 
	if(remove_similar_edges)
      FEVV::Filters::resolve_similar_edges(g);
    if ((nb_e != size_of_edges(g)) || (nb_f != size_of_faces(g)))
      ret = true;
    return ret;
  }

} // namespace Filters
} // namespace FEVV

