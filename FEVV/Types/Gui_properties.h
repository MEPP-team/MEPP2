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


namespace FEVV {
namespace Types {


enum class DisplayMode { NO_COLOR, VERTEX_COLOR, FACE_COLOR, TEXTURE };

/**
 * This is the datastructure used to store the GUI properties for a mesh.
 * It is used to trigger a change in GUI from within a plugin.
 */
struct GuiProperties
{
  std::string mesh_name;
  DisplayMode display_mode = DisplayMode::NO_COLOR;
  bool is_visible = true;
  bool is_selected = false;
};


} // namespace Types
} // namespace FEVV

