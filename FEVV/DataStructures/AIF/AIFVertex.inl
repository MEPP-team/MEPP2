#pragma once

namespace FEVV {
namespace DataStructures {
namespace AIF {

inline
AIFVertex::ptr_vertex
AIFVertex::New()
{
  ptr_vertex ptr(new self());
  return ptr;
}

inline
AIFVertex::ptr_vertex
AIFVertex::New(const self &other)
{
  ptr_vertex ptr(new self(other));
  return ptr;
}

inline
void
AIFVertex::Print() const
{
  std::cout << "vertex " << this << " at index " << GetIndex() << std::endl;
}

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV
