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

namespace FEVV {
namespace DataStructures {
namespace AIF {

/*
 * A class to store pointers to cells.
 * The index of a cell may vary in time according to the
 * remove operations that occur on the container.
 * When a cell is removed from the container, it is
 * not destroyed.
 */
template< typename T >
class AIFCellContainer
{
public:
  typedef typename std::vector< T >::iterator iterator;
  typedef typename std::vector< T >::const_iterator const_iterator;

  /*!
   *			Add element c to the end of container.
   */
  void add(T c)
  {
    c->SetIndex(container.size());
    container.push_back(c);
  }
  /*!
   *			Remove the element located at the i-th position in the
   *container.
   */
  std::size_t remove(std::size_t i)
  {
    std::size_t size = container.size();
    std::size_t lastId = size - 1;

    if(size == 0 || i > lastId)
      throw std::runtime_error(
          "In CellContainer::remove(): index out of range.");

    // Remove element #i.
    // If #i is not the last element, then replace #i by the last
    // element, then free the last slot.
    std::size_t tmp = container[i]->GetIndex();
    bool ascending_order = (tmp==i) && (container[lastId]->GetIndex() == lastId); // not exact, but avoid linear access for meshes without vertex reordering
    //for(std::size_t n=0; (n < lastId) && ascending_order; ++n)
    //  if (container[n]->GetIndex() > container[n + 1]->GetIndex())
    //  {
    //    ascending_order = false;
    //    break;
    //  }
    if (ascending_order)
    {
      container[i] = container[lastId];
      container[i]->SetIndex(i);
    }
    else
    {
      if (tmp == lastId)
      {
        container[i] = container[lastId];
      }
      else
      {
        std::size_t index_max = 0;
        for(std::size_t n=1; n <= lastId; ++n)
          if (container[n]->GetIndex() > container[index_max]->GetIndex())
          {
            index_max = n;
          }

        container[i] = container[lastId];

        container[index_max]->SetIndex(tmp);
      }
    }

    container.pop_back();

    return lastId;
  }

  /*!
   *			Random access to element located at the i-th position in the
   *container.
   */
  T &operator[](std::size_t idx) { return container[idx]; }
  /*!
   *			Random access to element located at the i-th position in the
   *container.
   */
  const T &operator[](std::size_t idx) const { return container[idx]; }
  /*!
   *			Return an iterator pointing to the beginning of the
   *container.
   */
  iterator begin(void) { return container.begin(); }
  /*!
   *			Return an iterator pointing to the end of the container.
   */
  iterator end(void) { return container.end(); }
  /*!
   *			Return a const iterator pointing to the beginning of the
   *container.
   */
  const_iterator begin(void) const { return container.begin(); }
  /*!
   *			Return a const iterator pointing to the end of the
   *container.
   */
  const_iterator end(void) const { return container.end(); }
  /*!
   *			Return a const iterator pointing to the beginning of the
   *container.
   */
  const_iterator cbegin(void) const { return container.cbegin(); }
  /*!
   *			Return a const iterator pointing to the end of the
   *container.
   */
  const_iterator cend(void) const { return container.cend(); }
  /*!
   *			Return the size of the container.
   */
  std::size_t size(void) const { return container.size(); }
  /*!
   *			Reserve a given size for the container.
   */
  void reserve(std::size_t n) { container.reserve(n); }

  /*!
   *			Remove the element pointed by the iterator in the
   *container.
   */
  void erase(const iterator iter) { container.erase(iter); }
  /*!
   *			Display debug information of the container.
   */
  void displayDebugInfos(const char *title)
  {
    std::cout << "Debug, " << title << ":\n";
    std::cout << " size=" << container.size() << '\n';
    for(std::size_t i = 0; i < container.size(); i++)
      std::cout << " [" << i << "]: " << container[i]->GetIndex() << ", "
                << container[i]->GetDegree() << '\n';
  }

private:
  //! the real cell container
  std::vector< T > container;
};

} // namespace AIF
} // namespace DataStructures
} // namespace FEVV

