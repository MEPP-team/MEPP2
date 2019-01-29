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

#include <vector>
#include <list>
#include <set>

namespace FEVV {
namespace Container {

/*
 * Insert a new element in container.
 * General case for STL containers.
 */
template< typename ContainerType, typename ContainedType >
void
insert(ContainerType &container, const ContainedType &element)
{
  container.insert(container.end(), element);
}

template< typename ContainerType, typename ContainedType >
void
erase(ContainerType &container, const ContainedType &element)
{
  typename ContainerType::iterator it;

  if((it = std::find(container.begin(), container.end(), element)) !=
     container.end())
    container.erase(it);
  else
    throw std::invalid_argument(
        "Helpers::erase element in container -> given element not found.");
}

template< typename ContainedType >
bool
is_sequential_container(const std::vector< ContainedType > &container)
{
  return true;
}

template< typename ContainedType >
bool
is_sequential_container(const std::list< ContainedType > &container)
{
  return true;
}

template< typename ContainedType >
bool
is_sequential_container(const std::set< ContainedType > &container)
{
  return false;
}

} // namespace Container
} // namespace FEVV

