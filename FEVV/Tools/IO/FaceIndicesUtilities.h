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
#include <set>
#include <stdexcept> // for std::invalid_argument
#include <algorithm> // for std::find

namespace FEVV {
namespace FaceIndicesUtils {
template< typename IndexType >
bool
are_all_face_indices_unique(const std::vector< IndexType > &in_face_indices)
{
  std::set< IndexType > s(in_face_indices.begin(), in_face_indices.end());
  return (in_face_indices.size() == s.size());
}

template< typename IndexType >
void
lines_indices_to_segments_indices(
    std::vector< std::vector< IndexType > > &in_out_lines_indices)
{
  size_t nb_poly_lines = in_out_lines_indices.size();
  for(size_t i = 0; i < nb_poly_lines; ++i)
  {
    size_t current_line_size = in_out_lines_indices[i].size();
    if(current_line_size <= 2)
    {
      if(current_line_size <= 1)
        throw std::invalid_argument(
            "lines_indices_to_segments_indices -> input lines indices not "
            "valid. Need at least 2 indices per line");

      continue;
    }

    std::vector< IndexType > line_segment;
    for(size_t j = 1; j < current_line_size;
        ++j) // the first segment is kept at that level
    {
      switch(line_segment.size())
      {
      case 0:
        line_segment.push_back(in_out_lines_indices[i][j]);
        break;
      case 2:
        line_segment.erase(line_segment.begin());
      case 1:
        line_segment.push_back(in_out_lines_indices[i][j]);
        in_out_lines_indices.push_back(line_segment);
      }
    }
    while(in_out_lines_indices[i].size() > 2)
      in_out_lines_indices[i].pop_back();
  }
}

template< typename IndexType >
void
face_indices_to_face_and_lines_indices(
    const std::vector< IndexType >
        &in_face_indices, // in order to work properly, each polyline length
                          // must be strictly less than the face degree minus 1
    std::vector< IndexType >
        &out_face_indices, // for the time being, we assume there is only one
                           // face described with the indices (only one loop)
    std::vector< std::vector< IndexType > >
        &out_lines_indices) // we can extract several dangling polylines
{
  size_t nb_e = in_face_indices.size();
  out_face_indices.clear();
  out_lines_indices.clear();

  std::vector< IndexType > new_line;
  bool found_elm_i = false, // true when an index is present 2 times and
                            // in-between these 2 occurences
      preserve_next = true, // keep current face index
      only_first_is_cleared = true;
  IndexType where_found_i = -1;
  for(IndexType i = 0; i < nb_e; ++i)
  {
    if(found_elm_i)
    {
      if(i == where_found_i)
      {
        found_elm_i = false;
        preserve_next = !preserve_next;
        where_found_i = -1;
        if(new_line.size() > 1)
        {
          out_lines_indices.push_back(new_line);
          new_line.clear();
        }
        continue; // to ensure that the same vertex is not added twice
      }
    }

    if(!preserve_next)
    { // what is not preserved for face indices is used for polylines...
      if(std::find(new_line.begin(), new_line.end(), in_face_indices[i]) ==
         new_line.end())
        new_line.push_back(in_face_indices[i]);
      if(i == nb_e - 1)
        out_lines_indices.push_back(new_line);
      continue;
    }

    for(IndexType j = i + 1; j < nb_e; ++j)
    {
      if(in_face_indices[i] == in_face_indices[j])
      {
        found_elm_i = true;
        where_found_i = j;
        if(j - i <= nb_e / 2)
        {
          preserve_next = false;
          new_line.push_back(in_face_indices[i]);
        }
        else
        {
          if(only_first_is_cleared && (i > 0))
          {
            out_face_indices.clear();
            only_first_is_cleared = false;
            if(out_lines_indices.size() > 0)
              out_lines_indices[out_lines_indices.size() - 1].insert(
                  out_lines_indices[out_lines_indices.size() - 1].begin(),
                  in_face_indices[i]);
          }
          else
            new_line.push_back(in_face_indices[i]);
        }
        break;
      }
    }

    if(std::find(out_face_indices.begin(),
                 out_face_indices.end(),
                 in_face_indices[i]) == out_face_indices.end())
      out_face_indices.push_back(in_face_indices[i]);
  }
}
} // namespace FaceIndicesUtils
} // namespace FEVV

