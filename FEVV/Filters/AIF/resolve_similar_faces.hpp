#pragma once

#include <boost/graph/graph_traits.hpp>
#include "FEVV/Operators/AIF/topology_predicates.hpp"

namespace FEVV {
namespace Filters {

/**
 * \brief Remove/resolve similar faces for the given mesh g.
 *
 * \tparam  MutableFaceIncidentGraph a Mesh type that provides a Model of the
 *          MutableFaceIncidentGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param   g The MutableFaceIncidentGraph instance.
 * \param   take_into_account_face_rientation A boolean to tell if the similar
 *          faces must have the same orientation.
 */
template< typename MutableFaceIncidentGraph >
static void
resolve_similar_faces(MutableFaceIncidentGraph &g,
                      bool take_into_account_face_rientation = false)
{
  typedef
      typename boost::graph_traits< MutableFaceIncidentGraph >::face_descriptor
          face_descriptor;
  std::vector< face_descriptor > faces_to_remove;
  /////////////////////////////////////////////////////////////////////////////
  auto faces_range_pair = faces(g);
  auto iter_f = faces_range_pair.first;
  for(; iter_f != faces_range_pair.second; ++iter_f)
  {
    if(std::find(faces_to_remove.begin(), faces_to_remove.end(), *iter_f) !=
       faces_to_remove.end())
      continue;
    auto iter_f_bis = iter_f;
    ++iter_f_bis;
    for(; iter_f_bis != faces_range_pair.second; ++iter_f_bis)
    {
      if(Operators::are_similar_faces(
             *iter_f, *iter_f_bis, g, take_into_account_face_rientation))
      {
        faces_to_remove.push_back(*iter_f_bis);
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  typename std::vector< face_descriptor >::iterator it_f(
      faces_to_remove.begin()),
      it_fe(faces_to_remove.end());
  for(; it_f != it_fe; ++it_f)
    remove_face(*it_f, g); // REDUNDANT FACES ARE REMOVED
  faces_to_remove.clear();
}

} // namespace Filters
} // namespace FEVV
