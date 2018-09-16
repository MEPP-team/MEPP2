#pragma once

// Warning: do NOT include this file outside of AIFEdge.hpp


namespace FEVV {
namespace DataStructures {
namespace AIF {

/**
 * \class	AIFEdgeComparator
 * \brief	This class represents an edge comparator used to compare AIFEdge
 * objects. \see		AIFEdge
 */
class AIFEdgeComparator
{
public:
  bool operator()(const AIFEdge &e1, const AIFEdge &e2);

  bool operator()(const AIFEdge::ptr_edge &PtrE1,
                  const AIFEdge::ptr_edge &PtrE2)
  {
    return operator()(*PtrE1, *PtrE2);
  }
};

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV


#include "FEVV/DataStructures/AIF/AIFEdgeComparator.inl"
