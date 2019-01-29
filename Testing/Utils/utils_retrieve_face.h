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

#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/include/next.hpp>

/*
 * \brief Utility functions that retrieves the face out
 *        of its index.
 * \tparam MutableFaceGraph Any mesh type having a boost::graph_traits
 *                    specialisation.
 * \param g           The graph which which to search for.
 * \param target_index The index of the target face.
 * \return            When found, the face with target index.
 */   
template<typename MutableFaceGraph>
typename boost::graph_traits<MutableFaceGraph>::face_descriptor
retrieve_face( MutableFaceGraph& g,
               int target_index )
{
  typedef boost::graph_traits<MutableFaceGraph>     GraphTraits;
  typedef typename GraphTraits::face_iterator       face_iterator;
  typedef typename GraphTraits::face_descriptor     face_descriptor;

  face_iterator f = faces(g).first;

  face_descriptor dst = *boost::next(f, target_index);

  return dst;
}
