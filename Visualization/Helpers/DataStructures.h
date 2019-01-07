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
