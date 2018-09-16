#pragma once

#include <boost/graph/graph_traits.hpp>

namespace FEVV {
namespace Filters {

/*!
 * \brief Function determining if mesh g is a pure triangular mesh.
 *
 * \tparam  FaceGraph a Mesh type that provides a Model of the
 *          FaceGraph Concept through a boost::graph_traits<>
 *          specialization.
 * \param   g The FaceGraph instance.
 * \return  true if all faces of g are triangles.
 */
template< typename FaceGraph >
inline bool
is_triangle_mesh(const FaceGraph &g)
{
  auto face_range_pair = faces(g);
  auto iter = face_range_pair.first;
  for(; iter != face_range_pair.second; ++iter)
  {
    if(degree(*iter, g) != 3)
      return false;
  }
  return true;
}

} // namespace Filters
} // namespace FEVV

