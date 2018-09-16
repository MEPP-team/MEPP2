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
