// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include <QMetaType>
#include "Base/Color.hpp"
#include "Visualization/BaseViewerOSG.h"
#include <osg/Node>


namespace FEVV {

class BaseViewerOSG;

namespace Helpers {

enum class DataType : unsigned char { MODEL = 1, GROUP = 2, DEBUG_OBJECT = 3 };

template< typename NodeType >
struct Model
{
  NodeType *node;
  std::string name;
  DataType type;
  unsigned int position;
  BaseViewerOSG *viewer;
};

} // namespace Helpers

} // namespace FEVV

Q_DECLARE_METATYPE(FEVV::Helpers::Model< osg::Node >)
