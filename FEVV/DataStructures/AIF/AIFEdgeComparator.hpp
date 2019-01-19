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
