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

#include <string>
#include <vector>

#include "FEVV/Types/Material.h"


namespace FEVV {
namespace Types {


/**
 * The representation of a mesh using std::vector
 */
template< typename coordP_type,
          typename coordN_type,
          typename coordT_type,
          typename coordC_type,
          typename index_type >
struct MVR /* Mesh Vector Representation */
{
  std::vector< std::vector< coordP_type > > points_coords;
  std::vector< std::vector< coordN_type > > normals_coords;
  std::vector< std::vector< coordT_type > > texture_coords;
  std::vector< std::vector< coordC_type > > vertex_color_coords;
  std::vector< std::vector< coordC_type > > face_color_coords;
  std::vector< std::vector< index_type > >  lines_indices, faces_indices;
  std::vector< std::vector< index_type > >  texture_face_indices;
  std::vector< std::vector< index_type > >  normal_face_indices;
  std::vector< std::vector< coordC_type > > points_colors, faces_colors,
                                            lines_colors;
  std::vector< index_type >            face_material;
  std::vector< FEVV::Types::Material > materials;
  std::vector< std::vector< std::vector< double > > > field_attributes;
  std::vector< std::string >           field_names;
};


} // namespace Types
} // namespace FEVV
