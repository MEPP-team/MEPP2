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
#ifndef __AIFEdge_cxx
#define __AIFEdge_cxx

#include "FEVV/DataStructures/AIF/AIFEdgeComparator.hpp"

namespace FEVV {
namespace DataStructures {
namespace AIF {

inline
AIFEdge::ptr_edge
AIFEdge::New()
{
  ptr_edge ptr(new self());
  return ptr;
}

inline
AIFEdge::ptr_edge
AIFEdge::New(const self &other)
{
  ptr_edge ptr(new self(other));
  return ptr;
}

inline
bool
AIFEdge::operator<(const self &other) const
{
  return AIFEdgeComparator().operator()(
      *this, other); // Use natural ordering defined in AIFEdgeComparator class
}

inline
bool
AIFEdge::operator==(const self &other) const
{
  return !(*this < other) && !(other < *this);
}

inline
void
AIFEdge::Print() const
{
  std::cout << "edge " << this << " at index " << GetIndex() << " ["
            << &(*get_first_vertex()) << "," << &(*get_second_vertex()) << "]"
            << std::endl;
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV

#endif
