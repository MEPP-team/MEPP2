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
