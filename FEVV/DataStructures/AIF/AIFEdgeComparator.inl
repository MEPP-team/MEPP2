
namespace FEVV {
namespace DataStructures {
namespace AIF {

inline
bool
AIFEdgeComparator::operator()(const AIFEdge &e1, const AIFEdge &e2)
{
  // compare the min vertex of each edge
  AIFEdge::ptr_vertex minVertexE1 = std::min< AIFEdge::ptr_vertex >(
      e1.m_incident_vertices[0], e1.m_incident_vertices[1]);
  AIFEdge::ptr_vertex minVertexE2 = std::min< AIFEdge::ptr_vertex >(
      e2.m_incident_vertices[0], e2.m_incident_vertices[1]);
  if(minVertexE1 != minVertexE2)
    return minVertexE1 < minVertexE2;

  // if min vertices are equal, then compare the max vertex of each edge
  AIFEdge::ptr_vertex maxVertexE1 = std::max< AIFEdge::ptr_vertex >(
      e1.m_incident_vertices[0], e1.m_incident_vertices[1]);
  AIFEdge::ptr_vertex maxVertexE2 = std::max< AIFEdge::ptr_vertex >(
      e2.m_incident_vertices[0], e2.m_incident_vertices[1]);
  return maxVertexE1 < maxVertexE2;
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV
